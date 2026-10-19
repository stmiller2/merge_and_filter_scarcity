#!/bin/bash
set -uo pipefail
# extglob must be enabled before the functions that use +(...)/?(...) are parsed
shopt -s extglob
# NOTE: -e is intentionally NOT set here. Failures are checked explicitly at
# each step so that one bad sample can be logged and skipped (via `continue`)
# without aborting the whole batch run.

# -----------------------------------------------------------------------------
# File layout (relative to data_filepath):
#   Fastq/                 raw inputs; read-only to the pipeline
#   work/<sample>/         per-sample scratch (PEAR output); deleted if space_saver=TRUE
#   csvs/good_reads/       good_reads_<sample>.csv
#   csvs/poor_reads/       poor_reads_<sample>.csv (.csv.gz if space_saver=TRUE)
#   csvs/info/             info_<sample>.csv
#   logs/                  run_progress.log, condor logs, params_used.env
# -----------------------------------------------------------------------------

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
  echo "Usage: merge_reads.sh <params_file>"
  exit 1
fi

# Load parameters
params_file="$1"
source "$params_file"

# Start clock
start_time=$(date --utc +%s)

# Setup
cd "${data_filepath}" || exit 1
fastq_dir="${data_filepath}/Fastq"
work_root="${data_filepath}/work"
csv_dir="${data_filepath}/csvs"
mkdir -p "${data_filepath}/logs"
exec &> "${data_filepath}/logs/run_progress.log"
mkdir -p "${csv_dir}"/{good_reads,poor_reads,info} "${work_root}"

qf3="${data_filepath}/pipeline/QF3.out"
if [ ! -x "$qf3" ]; then
    echo "Error: ${qf3} not found or not executable. Run process_ngs.sh to compile it."
    exit 1
fi
if [ "$merge" != "TRUE" ] && [ ! -x "${perl_filepath:-}" ]; then
    echo "Error: perl_filepath ('${perl_filepath:-}') not found or not executable on this execute node (needed to reverse-complement R2 when merge=FALSE)."
    exit 1
fi

# -----------------------------------------------------------------------------
# Derive sample names from Fastq/ filenames
# -----------------------------------------------------------------------------
derive_sample_names() {
    local dir="$1"
    local f base
    local -a names=()
    shopt -s nullglob
    # Match R1 or R2 so single-end runs on R2-only data are detected too.
    for f in "${dir}"/*_S[0-9]*_R[12]_001.fastq.gz \
             "${dir}"/*_S[0-9]*_R[12]_001.fastq; do
        base=$(basename "$f")
        base="${base%_S[0-9]*_R[12]_001.fastq.gz}"
        base="${base%_S[0-9]*_R[12]_001.fastq}"
        names+=("$base")
    done
    shopt -u nullglob
    [[ "${#names[@]}" -gt 0 ]] || return 0
    printf '%s\n' "${names[@]}" | awk 'NF && !seen[$0]++'
}

mapfile -t sample_names_array < <(derive_sample_names "${fastq_dir}")
if [ "${#sample_names_array[@]}" -eq 0 ]; then
    echo "Error: no fastq files found in ${fastq_dir} matching {sample}_S{#}_R[12]_001.fastq[.gz] (with or without a _L001_-style lane segment)."
    exit 1
fi
echo "Detected samples: ${sample_names_array[*]}"

# Log helper
log() {
    local type="$1"
    local msg="$2"
    if [ "$type" == "start" ]; then
        echo -e "\n$(date '+%I:%M%p') -- $msg"
    else
        echo -e "$(date '+%I:%M%p') -- $msg"
    fi
}

# Collect sample $1's fastqs in Fastq/ into the array sample_fastqs.
# Strict pattern: {sample}_S{digits}[_L{digits}]_R{1|2}_001.fastq[.gz], so
# sample 'A' does not pick up files belonging to sample 'A_2'.
find_sample_fastqs() {
    local s="$1"
    shopt -s nullglob
    sample_fastqs=( "${fastq_dir}/${s}"_S+([0-9])?(_L+([0-9]))_R[12]_001.fastq?(.gz) )
    shopt -u nullglob
}

# -----------------------------------------------------------------------------
# Read streamers. Both write QF3's input format to stdout, one read per line:
#   SEQUENCE QUALITY
# Nothing is written to disk; the output is piped straight into QF3.
# -----------------------------------------------------------------------------

# Reads from one or more fastq[.gz] files as-is (PEAR output, single-end R1).
fastq_to_qf() {
    zcat -f -- "$@" | paste -d '\t' - - - - | cut -f2,4 --output-delimiter=' '
}

# Reverse-complement R2 (quality string reversed too), optionally prepending R1.
#   revcomp_to_qf <n_r1> [R1 files...] <R2 files...>
#   n_r1 = 0  -> revcomp(R2)            (single-end R2)
#   n_r1 > 0  -> R1 + revcomp(R2)       (paired-end, merge=FALSE)
# Each input is decompressed once; R1/R2 read-count mismatches are an error.
revcomp_perl='
use strict; use warnings;
my $n1 = shift @ARGV;
my @r1 = splice(@ARGV, 0, $n1);
my @r2 = @ARGV;
die "revcomp_to_qf: no R2 files given\n" unless @r2;
sub zopen { open(my $fh, "-|", "zcat", "-f", "--", @_) or die "cannot start zcat: $!\n"; return $fh; }
sub rec {
    my $fh = shift;
    my $h = <$fh>;
    return unless defined $h;
    my $s = <$fh>; my $p = <$fh>; my $q = <$fh>;
    die "truncated fastq record\n" unless defined $q;
    chomp($s, $q);
    return ($s, $q);
}
my $fh2 = zopen(@r2);
my $fh1 = @r1 ? zopen(@r1) : undef;
while (1) {
    my ($s2, $q2) = rec($fh2);
    my ($s1, $q1) = $fh1 ? rec($fh1) : ();
    if (!defined $s2) {
        die "R1 has more reads than R2\n" if defined $s1;
        last;
    }
    die "R2 has more reads than R1\n" if $fh1 && !defined $s1;
    my $rs = reverse $s2;
    $rs =~ tr/ACGTNacgtn/TGCANTGCAN/;
    my $rq = reverse $q2;
    print $fh1 ? "$s1$rs $q1$rq\n" : "$rs $rq\n";
}
close($fh2) or die "zcat failed reading R2 input\n";
if ($fh1) { close($fh1) or die "zcat failed reading R1 input\n"; }
'
revcomp_to_qf() {
    "${perl_filepath}" -e "$revcomp_perl" "$@"
}

# Count reads and write info.csv in the current (work) directory.
write_stats() {
    local f good poor total unassembled discarded pct
    local -a needed=(good_reads.csv poor_reads.csv)
    [ "$merge" == "TRUE" ] && needed+=(combined.unassembled.forward.fastq combined.discarded.fastq)
    for f in "${needed[@]}"; do
        [ -f "$f" ] || { echo "ERROR: ${f} not found, cannot compute stats."; return 1; }
    done
    good=$(wc -l < good_reads.csv)
    poor=$(wc -l < poor_reads.csv)
    if [ "$merge" == "TRUE" ]; then
        unassembled=$(( $(wc -l < combined.unassembled.forward.fastq) / 4 ))
        discarded=$(( $(wc -l < combined.discarded.fastq) / 4 ))
        total=$((unassembled + discarded + poor + good))
        echo "total_reads,unassembled_reads,discarded_reads,poor_reads,good_reads" > info.csv
        echo "${total},${unassembled},${discarded},${poor},${good}" >> info.csv
    else
        total=$((poor + good))
        echo "total_reads,poor_reads,good_reads" > info.csv
        echo "${total},${poor},${good}" >> info.csv
    fi
    pct=$(awk -v g="$good" -v t="$total" 'BEGIN {printf "%.1f", (t > 0 ? 100 * g / t : 0)}')
    echo "           $(printf "%'d" "$good") / $(printf "%'d" "$total") reads passing all filters (${pct}%)"
}

failed_samples=()

# Iterate through samples
for i in "${sample_names_array[@]}"; do
	echo -e "\033[1m$(printf %80s |tr " " "=")\033[0m\n"
	echo -e "\033[1m$(printf %$(((80-(18+${#i}))/2))s |tr " " " ")PROCESSING SAMPLE $i\033[0m\n"
	echo -e "\033[1m$(printf %80s |tr " " "=")\033[0m"

    # Fresh per-sample work directory (clears leftovers from a failed earlier run)
    work_dir="${work_root}/${i}"
    rm -rf -- "${work_dir}"
    if ! mkdir -p "${work_dir}" || ! cd "${work_dir}"; then
        echo "ERROR: cannot create work/${i}/, skipping sample."
        failed_samples+=("$i"); continue
    fi

    find_sample_fastqs "$i"
    if [ "${#sample_fastqs[@]}" -eq 0 ]; then
        echo "ERROR: no input files found for sample ${i} in Fastq/, skipping."
        failed_samples+=("$i"); continue
    fi

    # Split by read (multiple lanes of the same read are concatenated downstream).
    r1_files=(); r2_files=()
    for f in "${sample_fastqs[@]}"; do
        if [[ "$f" == *_R2_001.fastq* ]]; then r2_files+=("$f"); else r1_files+=("$f"); fi
    done

    # Choose how reads reach QF3.
    # qf_ends tells QF3 which side(s) of the final read are sequencing-read 5'
    # ends, where end_slop applies (L = left, R = right, LR = both).
    if [ "$merge" == "TRUE" ]; then
        qf_ends="LR"    # PEAR output: R1 5' end on left, revcomp'd R2 5' end on right
        if [ "${#r1_files[@]}" -ne 1 ] || [ "${#r2_files[@]}" -ne 1 ]; then
            echo "ERROR: merge=TRUE needs exactly one R1 and one R2 file for sample ${i} (found ${#r1_files[@]} R1, ${#r2_files[@]} R2), skipping."
            failed_samples+=("$i"); continue
        fi
        log start "MERGING PAIRED-END READS"
        pear_output=$("${pear_filepath}" -f "${r1_files[0]}" -r "${r2_files[0]}" -o combined \
            -y "${memory}" -j "${cpus}" -v "${pear_overlap}" -g "${pear_stattest}" -p "${pear_pvalue}" 2>&1)
        pear_status=$?
        echo "$pear_output" | grep -E -m 3 "Assembled reads|Discarded reads|Not assembled reads" | awk 'NF' | sed 's/^/           /'
        if [ $pear_status -ne 0 ]; then
            echo "ERROR: PEAR failed for sample ${i} (exit code ${pear_status}):"
            echo "$pear_output"
            failed_samples+=("$i"); continue
        fi
        log stop "READS MERGED"
        qf_input=(fastq_to_qf combined.assembled.fastq)
    elif [ "$singleend" == "FALSE" ]; then
        qf_ends="LR"    # R1 + revcomp(R2): 5' ends on both outer edges, 3' ends at the junction
        if [ "${#r1_files[@]}" -eq 0 ] || [ "${#r2_files[@]}" -eq 0 ]; then
            echo "ERROR: paired-end run but sample ${i} is missing R1 or R2, skipping."
            failed_samples+=("$i"); continue
        fi
        log start "CONCATENATING R1 + REVERSE COMPLEMENT OF R2 (streamed into quality filter)"
        qf_input=(revcomp_to_qf "${#r1_files[@]}" "${r1_files[@]}" "${r2_files[@]}")
    else
        # Single-end
        if [ "${#r1_files[@]}" -gt 0 ] && [ "${#r2_files[@]}" -gt 0 ]; then
            echo "ERROR: sample ${i} has both R1 and R2 fastq files in single-end mode, skipping."
            failed_samples+=("$i"); continue
        elif [ "${#r2_files[@]}" -gt 0 ]; then
            qf_ends="R"     # revcomp'd R2: its 5' end is now on the right
            log start "SINGLE-END READ 2: REVERSE COMPLEMENTING (streamed into quality filter)"
            qf_input=(revcomp_to_qf 0 "${r2_files[@]}")
        else
            qf_ends="L"     # R1 as sequenced: 5' end on the left
            qf_input=(fastq_to_qf "${r1_files[@]}")
        fi
    fi

    # Quality filtering (reads are streamed in; writes good_reads.csv / poor_reads.csv here)
    log start "FILTERING FOR HIGH QUALITY READS"
    if ! "${qf_input[@]}" | "$qf3" "$qf_ends"; then
        echo "ERROR: formatting or quality filtering failed for sample ${i}, skipping."
        failed_samples+=("$i"); continue
    fi
    log stop "QUALITY FILTERING COMPLETE"

    # Stats calculation
    log start "GETTING ASSEMBLY AND FILTERING STATISTICS"
    if ! write_stats; then
        echo "ERROR: stats calculation failed for sample ${i}, skipping."
        failed_samples+=("$i"); continue
    fi
    log stop "READ FATES WRITTEN TO INFO.CSV"

    # Move results into csvs/
    log start "WRITING RESULTS TO csvs/"
    rm -f -- "${csv_dir}/poor_reads/poor_reads_${i}.csv" "${csv_dir}/poor_reads/poor_reads_${i}.csv.gz"
    poor_out="poor_reads.csv"
    if [ "$space_saver" == "TRUE" ]; then
        # Level 4: on sequence text, the default (-6) is ~7x slower for ~5% smaller files
        if command -v pigz > /dev/null 2>&1; then
            pigz -4 -f -p "${cpus}" poor_reads.csv
        else
            gzip -4 -f poor_reads.csv
        fi || { echo "ERROR: failed to gzip poor_reads.csv for sample ${i}."; failed_samples+=("$i"); continue; }
        poor_out="poor_reads.csv.gz"
    fi
    if mv -- good_reads.csv "${csv_dir}/good_reads/good_reads_${i}.csv" &&
       mv -- "${poor_out}" "${csv_dir}/poor_reads/poor_reads_${i}.${poor_out#poor_reads.}" &&
       mv -- info.csv "${csv_dir}/info/info_${i}.csv"; then
        log stop "RESULTS WRITTEN"
    else
        echo "ERROR: moving output files to csvs/ failed for sample ${i}."
        failed_samples+=("$i"); continue
    fi

    # Cleanup. Failed samples never reach this point, so their work/ dir is kept.
    cd "${data_filepath}" || exit 1
    if [ "$space_saver" == "TRUE" ]; then
        rm -rf -- "${work_dir}"
        log stop "SPACE SAVER: REMOVED work/${i}/"
    else
        rmdir -- "${work_dir}" 2> /dev/null || true    # only removed if empty (merge=FALSE modes)
    fi

    echo -e "\nSAMPLE $i COMPLETE\n\n"
done

# Remove work/ if nothing is left in it
cd "${data_filepath}" || exit 1
rmdir -- "${work_root}" 2> /dev/null || true

# Finish
end_time=$(date --utc +%s)
runtime=$((end_time - start_time))
echo -e "\n$(printf %80s |tr " " "=")"
if [ "${#failed_samples[@]}" -gt 0 ]; then
    echo "Failed samples: ${failed_samples[*]} (intermediate files kept in work/ for debugging)"
    log stop "Done WITH ERRORS. Total runtime: $(date -u -d "@${runtime}" +%H:%M:%S)"
    exit 1
else
    log stop "Done. All samples completed successfully. Total runtime: $(date -u -d "@${runtime}" +%H:%M:%S)"
fi
