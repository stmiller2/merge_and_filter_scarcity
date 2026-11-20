#!/bin/bash
set -uo pipefail
# extglob must be enabled before the functions that use +(...)/?(...) are parsed
shopt -s extglob
# NOTE: -e is intentionally NOT set here. Failures are checked explicitly at
# each step so that one bad sample can be logged and skipped (via `continue`)
# without aborting the whole batch run.

# -----------------------------------------------------------------------------
# File layout (relative to pipeline_dir):
#   Fastq/                 raw inputs; read-only to the pipeline
#   work/<sample>/         per-sample scratch (PEAR output); deleted if space_saver=TRUE
#   csvs/good_reads/       good_reads_<sample>.csv
#   csvs/poor_reads/       poor_reads_<sample>.csv (.csv.gz if space_saver=TRUE)
#   csvs/info/             info_<sample>.csv
#   logs/                  run_progress.log, condor logs, params_used.env
# -----------------------------------------------------------------------------

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: merge_reads.sh <params_file> <pipeline_dir>"
  exit 1
fi

# Load parameters. pipeline_dir is passed in by process_ngs.sh (the folder it
# was run from), not read from the params file.
params_file="$1"
pipeline_dir="$2"
source "$params_file"

# Start clock
start_time=$(date --utc +%s)

# Setup
cd "${pipeline_dir}" || exit 1
fastq_dir="${pipeline_dir}/Fastq"
work_root="${pipeline_dir}/work"
csv_dir="${pipeline_dir}/csvs"
mkdir -p "${pipeline_dir}/logs"
exec &> "${pipeline_dir}/logs/run_progress.log"
mkdir -p "${csv_dir}"/{good_reads,poor_reads,info} "${work_root}"

qf3="${pipeline_dir}/pipeline/QF3.out"
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
# Also sets the globals n_total, n_unassembled, n_discarded, n_e1, n_e2, n_e3,
# n_good for the end-of-run summary.
write_stats() {
    local f good poor total unassembled discarded pct poor_counts
    local -a needed=(good_reads.csv poor_reads.csv)
    [ "$merge" == "TRUE" ] && needed+=(combined.unassembled.forward.fastq combined.discarded.fastq)
    for f in "${needed[@]}"; do
        [ -f "$f" ] || { echo "ERROR: ${f} not found, cannot compute stats."; return 1; }
    done
    good=$(wc -l < good_reads.csv)
    # One pass over poor_reads.csv: total lines plus counts per failure code (column 2)
    poor_counts=$(awk -F, '{c[$2]++} END {printf "%d %d %d %d\n", NR, c["E1"], c["E2"], c["E3"]}' poor_reads.csv) ||
        { echo "ERROR: could not count poor_reads.csv."; return 1; }
    read -r poor n_e1 n_e2 n_e3 <<< "$poor_counts"
    if [ "$merge" == "TRUE" ]; then
        unassembled=$(( $(wc -l < combined.unassembled.forward.fastq) / 4 ))
        discarded=$(( $(wc -l < combined.discarded.fastq) / 4 ))
        total=$((unassembled + discarded + poor + good))
        echo "total_reads,unassembled_reads,discarded_reads,poor_reads,good_reads" > info.csv
        echo "${total},${unassembled},${discarded},${poor},${good}" >> info.csv
    else
        unassembled=0
        discarded=0
        total=$((poor + good))
        echo "total_reads,poor_reads,good_reads" > info.csv
        echo "${total},${poor},${good}" >> info.csv
    fi
    if [ "$poor" -ne $((n_e1 + n_e2 + n_e3)) ]; then
        echo "WARNING: poor_reads.csv has $((poor - n_e1 - n_e2 - n_e3)) line(s) without an E1/E2/E3 code."
    fi
    n_total=$total; n_unassembled=$unassembled; n_discarded=$discarded; n_good=$good
    pct=$(awk -v g="$good" -v t="$total" 'BEGIN {printf "%.1f", (t > 0 ? 100 * g / t : 0)}')
    echo "           $(printf "%'d" "$good") / $(printf "%'d" "$total") reads passing all filters (${pct}%)"
}

# -----------------------------------------------------------------------------
# End-of-run summary: where every input read ended up, per sample and overall
# -----------------------------------------------------------------------------
pct_of() {
    awk -v n="$1" -v t="$2" 'BEGIN {printf "%.2f", (t > 0 ? 100 * n / t : 0)}'
}

# summary_line <label> <count> <total>
summary_line() {
    printf "  %-46s %15s  %6s%%\n" "$1" "$(printf "%'d" "$2")" "$(pct_of "$2" "$3")"
}

# summary_block <title> <total> <unassembled> <discarded> <e1> <e2> <e3> <good>
summary_block() {
    local title="$1" total="$2" unas="$3" disc="$4" e1="$5" e2="$6" e3="$7" good="$8"
    echo ""
    echo "${title}: $(printf "%'d" "$total") input reads"
    summary_line "Good (passed all filters)" "$good" "$total"
    if [ "$merge" == "TRUE" ]; then
        summary_line "PEAR: not assembled" "$unas" "$total"
        summary_line "PEAR: discarded" "$disc" "$total"
    fi
    summary_line "Poor E1: >=1 base below Q${q_floor}" "$e1" "$total"
    summary_line "Poor E2: <=${cutoff_pct} of bases at >=Q${q_cutoff}" "$e2" "$total"
    if [ "$slop_on" == "TRUE" ] || [ "$e3" -gt 0 ]; then
        summary_line "Poor E3: too short after removing end slop" "$e3" "$total"
    fi
}

print_run_summary() {
    local row s total unas disc e1 e2 e3 good poor_ext
    local sum_total=0 sum_unas=0 sum_disc=0 sum_e1=0 sum_e2=0 sum_e3=0 sum_good=0
    slop_on=FALSE
    if [[ "${end_slop:-FALSE}" =~ ^[0-9]+$ ]] && (( 10#${end_slop} > 0 )); then slop_on=TRUE; fi

    echo -e "\n\033[1m$(printf %80s |tr " " "=")\033[0m"
    echo -e "\033[1m$(printf %33s)RUN SUMMARY\033[0m"
    echo -e "\033[1m$(printf %80s |tr " " "=")\033[0m"
    echo "Percentages are of each sample's input reads (overall: of all input reads)."

    if [ "${#summary_rows[@]}" -gt 0 ]; then
        for row in "${summary_rows[@]}"; do
            IFS=$'\t' read -r s total unas disc e1 e2 e3 good <<< "$row"
            summary_block "Sample ${s}" "$total" "$unas" "$disc" "$e1" "$e2" "$e3" "$good"
            sum_total=$((sum_total + total)); sum_unas=$((sum_unas + unas)); sum_disc=$((sum_disc + disc))
            sum_e1=$((sum_e1 + e1)); sum_e2=$((sum_e2 + e2)); sum_e3=$((sum_e3 + e3)); sum_good=$((sum_good + good))
        done
        summary_block "OVERALL (${#summary_rows[@]} of ${#sample_names_array[@]} samples)" \
            "$sum_total" "$sum_unas" "$sum_disc" "$sum_e1" "$sum_e2" "$sum_e3" "$sum_good"
    else
        echo -e "\nNo samples completed; nothing to summarize."
    fi
    if [ "${#failed_samples[@]}" -gt 0 ]; then
        echo -e "\nFailed samples (not included above): ${failed_samples[*]}"
    fi

    [ "$space_saver" == "TRUE" ] && poor_ext="csv.gz" || poor_ext="csv"
    echo -e "\nWhere these reads are:"
    echo "  Good                 csvs/good_reads/good_reads_<sample>.csv"
    echo "  Poor (E1/E2/E3)      csvs/poor_reads/poor_reads_<sample>.${poor_ext} (code in column 2)"
    if [ "$merge" == "TRUE" ]; then
        if [ "$space_saver" == "TRUE" ]; then
            echo "  PEAR not assembled   deleted with work/<sample>/ (space_saver=TRUE)"
            echo "  PEAR discarded       deleted with work/<sample>/ (space_saver=TRUE)"
        else
            echo "  PEAR not assembled   work/<sample>/combined.unassembled.{forward,reverse}.fastq"
            echo "  PEAR discarded       work/<sample>/combined.discarded.fastq"
        fi
    fi
    echo "  Per-sample counts    csvs/info/info_<sample>.csv"
}

failed_samples=()
summary_rows=()

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
    cd "${pipeline_dir}" || exit 1
    if [ "$space_saver" == "TRUE" ]; then
        rm -rf -- "${work_dir}"
        log stop "SPACE SAVER: REMOVED work/${i}/"
    else
        rmdir -- "${work_dir}" 2> /dev/null || true    # only removed if empty (merge=FALSE modes)
    fi

    summary_rows+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "$i" "$n_total" "$n_unassembled" "$n_discarded" "$n_e1" "$n_e2" "$n_e3" "$n_good")")

    echo -e "\nSAMPLE $i COMPLETE\n\n"
done

# Remove work/ if nothing is left in it
cd "${pipeline_dir}" || exit 1
rmdir -- "${work_root}" 2> /dev/null || true

print_run_summary

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
