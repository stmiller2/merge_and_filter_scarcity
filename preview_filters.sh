#!/bin/bash
set -euo pipefail
# extglob must be enabled before the functions that use +(...)/?(...) are parsed
shopt -s extglob

usage() {
    cat <<'EOF'
Usage: ./preview_filters.sh <params_file.env> [options]

Estimates the % of reads kept / lost (E1, E2, E3, and PEAR when merge=TRUE)
for a fixed grid of quality-filter settings:
  q_floor    0-20 (step 5)
  q_cutoff   0-30 (step 5)
  cutoff_pct 0.50-0.95 (step 0.05)
Other settings (merge, singleend, end_slop, PEAR params) come from the params file.

Options:
  -n N        reads (pairs, if paired-end) to sample per sample   [default 50000]
  -m M        only scan the first M reads of each sample. Faster for very large
              files, but slightly less representative (the first reads in a
              fastq come from the first tiles imaged)              [default: whole file]
  -s SAMPLE   preview only this sample
  -a          also print the full scan for every sample (default: all samples
              pooled, plus a per-sample table at the current settings)
  --no-color  plain text output
  -h          show this help

Output is printed to the terminal and saved to logs/:
  filter_preview.txt          full report, every sample (plain text)
  filter_preview_grid.csv     counts for every sample x setting (long format)
  filter_preview_samples.csv  per-sample read counts and Q values present
Tip: page through long output with   ./preview_filters.sh params.env | less -R
EOF
}

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
params_file=""
n_sample=50000
max_scan=0
only_sample=""
show_all=0
color_opt=1
while [ "$#" -gt 0 ]; do
    case "$1" in
        -n) n_sample="${2:-}"; shift 2 ;;
        -m) max_scan="${2:-}"; shift 2 ;;
        -s) only_sample="${2:-}"; shift 2 ;;
        -a) show_all=1; shift ;;
        --no-color) color_opt=0; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
        *)
            if [ -n "$params_file" ]; then echo "Unexpected argument: $1" >&2; usage >&2; exit 1; fi
            params_file="$1"; shift ;;
    esac
done

if [ -z "$params_file" ]; then usage >&2; exit 1; fi
if [ ! -f "$params_file" ]; then
    echo "Error: Params file '$params_file' not found." >&2
    exit 1
fi
source "$params_file"

pipeline_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
fastq_dir="${pipeline_dir}/Fastq"
log_dir="${pipeline_dir}/logs"

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------
errors=()
[[ "$n_sample" =~ ^[1-9][0-9]*$ ]] || errors+=("-n must be a positive integer (got '${n_sample}').")
[[ "$max_scan" =~ ^[0-9]+$ ]] || errors+=("-m must be a positive integer (got '${max_scan}').")
if [[ "$max_scan" =~ ^[0-9]+$ ]] && [ "$max_scan" -gt 0 ] && [[ "$n_sample" =~ ^[0-9]+$ ]] && [ "$max_scan" -lt "$n_sample" ]; then
    errors+=("-m (${max_scan}) is smaller than -n (${n_sample}).")
fi
[ -d "$fastq_dir" ] || errors+=("Fastq directory not found at '${fastq_dir}'.")
for v in merge singleend; do
    [[ "${!v:-}" == "TRUE" || "${!v:-}" == "FALSE" ]] || errors+=("Parameter '${v}' must be TRUE or FALSE (got '${!v:-}').")
done
if [[ "${merge:-}" == "TRUE" && "${singleend:-}" == "TRUE" ]]; then
    errors+=("merge=TRUE implies paired-end reads, but singleend is also TRUE. Set one of them to FALSE.")
fi
for v in q_floor q_cutoff; do
    [[ "${!v:-}" =~ ^[0-9]+$ ]] || errors+=("Parameter '${v}' must be a positive integer (got '${!v:-}').")
done
[[ "${cutoff_pct:-}" =~ ^[0-1](\.[0-9]+)?$ ]] || errors+=("cutoff_pct should be a decimal between 0 and 1 (got '${cutoff_pct:-}').")
if [[ "${end_slop:-}" == "FALSE" ]]; then
    end_slop_value=0
elif [[ "${end_slop:-}" =~ ^[0-9]+$ ]]; then
    end_slop_value=$((10#${end_slop}))
else
    errors+=("Parameter 'end_slop' must be FALSE or a non-negative integer (got '${end_slop:-}').")
fi
[ -x "${compiler_filepath:-}" ] || errors+=("compiler_filepath ('${compiler_filepath:-}') does not exist or is not executable.")
[ -x "${perl_filepath:-}" ] || errors+=("perl_filepath ('${perl_filepath:-}') does not exist or is not executable (the preview uses perl for subsampling).")
if [[ "${merge:-}" == "TRUE" ]]; then
    [ -x "${pear_filepath:-}" ] || errors+=("pear_filepath ('${pear_filepath:-}') does not exist or is not executable (required when merge=TRUE).")
    [[ "${memory:-}" =~ ^[0-9]+[KMGkmg]$ ]] || errors+=("memory should be digits followed by K/M/G for PEAR's -y flag, e.g. '4G' (got '${memory:-}').")
fi

# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------
derive_sample_names() {
    local dir="$1"
    local f base
    local -a names=()
    shopt -s nullglob
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

find_sample_fastqs() {
    local s="$1"
    shopt -s nullglob
    sample_fastqs=( "${fastq_dir}/${s}"_S+([0-9])?(_L+([0-9]))_R[12]_001.fastq?(.gz) )
    shopt -u nullglob
}

sample_names_array=()
if [ -d "$fastq_dir" ]; then
    mapfile -t sample_names_array < <(derive_sample_names "$fastq_dir")
fi
if [ "${#sample_names_array[@]}" -eq 0 ]; then
    errors+=("No sample fastq files found in '${fastq_dir}'.")
elif [ -n "$only_sample" ]; then
    found=0
    for s in "${sample_names_array[@]}"; do
        if [ "$s" == "$only_sample" ]; then found=1; fi
    done
    if [ "$found" -eq 1 ]; then
        sample_names_array=("$only_sample")
    else
        errors+=("Sample '${only_sample}' not found in Fastq/. Detected: ${sample_names_array[*]}")
    fi
fi

if [ "${#errors[@]}" -gt 0 ]; then
    echo "================ Preview Checks Failed ================" >&2
    for err in "${errors[@]}"; do echo " - $err" >&2; done
    echo "=======================================================" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/qf_preview.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

qfscan="${tmp_dir}/QFscan.out"
if ! "${compiler_filepath}" -O2 -o "$qfscan" "${pipeline_dir}/pipeline/QFscan.cpp"; then
    echo "Error: compilation of pipeline/QFscan.cpp failed." >&2
    exit 1
fi

mkdir -p "$log_dir"
grid_csv="${log_dir}/filter_preview_grid.csv"
samples_csv="${log_dir}/filter_preview_samples.csv"
report_txt="${log_dir}/filter_preview.txt"
rm -f -- "$grid_csv" "$samples_csv" "$report_txt"

# Reservoir sampler (fixed seed). Modes: fastq (R1/R2 for PEAR), concat (R1 + revcomp(R2)),
# rc (revcomp(R2)), plain (R1). Prints "<reads scanned> <reads kept>".
sampler_perl='
use strict; use warnings;
my ($n, $max, $mode, $out, $n1, @files) = @ARGV;
my @r1 = splice(@files, 0, $n1);
my @r2 = @files;
srand(20261120);
sub zopen { open(my $fh, "-|", "zcat", "-f", "--", @_) or die "cannot start zcat: $!\n"; return $fh; }
sub rec {
    my $fh = shift;
    my $h = <$fh>;
    return unless defined $h;
    my $s = <$fh>; my $p = <$fh>; my $q = <$fh>;
    die "truncated fastq record\n" unless defined $q;
    chomp($h, $s, $q);
    return ($h, $s, $q);
}
sub rc { my ($s, $q) = @_; my $rs = reverse $s; $rs =~ tr/ACGTNacgtn/TGCANTGCAN/; return ($rs, scalar reverse $q); }
my $paired = ($mode eq "fastq" || $mode eq "concat");
my $fa = zopen($mode eq "rc" ? @r2 : @r1);
my $fb = $paired ? zopen(@r2) : undef;
my (@ka, @kb);
my $i = 0;
my $early = 0;
while (1) {
    if ($max > 0 && $i >= $max) { $early = 1; last; }
    my @a = rec($fa);
    my @b = $paired ? rec($fb) : ();
    if (!@a) { die "R2 has more reads than R1\n" if @b; last; }
    die "R1 has more reads than R2\n" if $paired && !@b;
    my $slot = $i < $n ? $i : int(rand($i + 1));
    $i++;
    next if $slot >= $n;
    if ($mode eq "fastq") {
        $ka[$slot] = "$a[0]\n$a[1]\n+\n$a[2]\n";
        $kb[$slot] = "$b[0]\n$b[1]\n+\n$b[2]\n";
    } elsif ($mode eq "concat") {
        my ($rs, $rq) = rc($b[1], $b[2]);
        $ka[$slot] = "$a[1]$rs $a[2]$rq\n";
    } elsif ($mode eq "rc") {
        my ($rs, $rq) = rc($a[1], $a[2]);
        $ka[$slot] = "$rs $rq\n";
    } else {
        $ka[$slot] = "$a[1] $a[2]\n";
    }
}
# With -m, zcat is closed mid-stream, so its exit status is ignored
if (!close($fa) && !$early) { die "zcat failed reading input\n"; }
if ($fb && !close($fb) && !$early) { die "zcat failed reading R2 input\n"; }
sub write_out { my ($path, $recs) = @_; open(my $o, ">", $path) or die "cannot write $path: $!\n"; print $o @$recs; close($o) or die "error writing $path\n"; }
if ($mode eq "fastq") { write_out("${out}_R1.fastq", \@ka); write_out("${out}_R2.fastq", \@kb); }
else { write_out("${out}.qf", \@ka); }
print "$i " . scalar(@ka) . "\n";
'

fmt_n() { printf "%'d" "$1"; }

# -----------------------------------------------------------------------------
# Scan each sample
# -----------------------------------------------------------------------------
scan_note=""
if [ "$max_scan" -gt 0 ]; then scan_note=" from the first $(fmt_n "$max_scan")"; fi
echo "Sampling up to $(fmt_n "$n_sample") reads per sample${scan_note}." >&2
skipped=()
k=0
for i in "${sample_names_array[@]}"; do
    k=$((k + 1))
    t0=$(date +%s)
    printf "  [%d/%d] %s ... " "$k" "${#sample_names_array[@]}" "$i" >&2

    find_sample_fastqs "$i"
    r1_files=(); r2_files=()
    for f in ${sample_fastqs[@]+"${sample_fastqs[@]}"}; do
        if [[ "$f" == *_R2_001.fastq* ]]; then r2_files+=("$f"); else r1_files+=("$f"); fi
    done

    err=""
    if [ "$merge" == "TRUE" ]; then
        qf_ends="LR"; mode="fastq"
        if [ "${#r1_files[@]}" -ne 1 ] || [ "${#r2_files[@]}" -ne 1 ]; then
            err="merge=TRUE needs exactly one R1 and one R2 file (found ${#r1_files[@]} R1, ${#r2_files[@]} R2)"
        fi
    elif [ "$singleend" == "FALSE" ]; then
        qf_ends="LR"; mode="concat"
        if [ "${#r1_files[@]}" -eq 0 ] || [ "${#r2_files[@]}" -eq 0 ]; then
            err="paired-end run but R1 or R2 is missing"
        fi
    else
        if [ "${#r1_files[@]}" -gt 0 ] && [ "${#r2_files[@]}" -gt 0 ]; then
            err="single-end run but both R1 and R2 files are present"
        elif [ "${#r2_files[@]}" -gt 0 ]; then
            qf_ends="R"; mode="rc"
        else
            qf_ends="L"; mode="plain"
        fi
    fi
    if [ -n "$err" ]; then
        echo "skipped: ${err}" >&2
        skipped+=("$i"); continue
    fi

    out="${tmp_dir}/sub"
    if ! counts=$("${perl_filepath}" -e "$sampler_perl" "$n_sample" "$max_scan" "$mode" "$out" \
                  "${#r1_files[@]}" ${r1_files[@]+"${r1_files[@]}"} ${r2_files[@]+"${r2_files[@]}"}); then
        echo "skipped: subsampling failed" >&2
        skipped+=("$i"); continue
    fi
    read -r n_scanned n_kept <<< "$counts"

    pear_un=0; pear_disc=0
    if [ "$merge" == "TRUE" ]; then
        if ! "${pear_filepath}" -f "${out}_R1.fastq" -r "${out}_R2.fastq" -o "${tmp_dir}/pear" \
               -y "${memory}" -j 2 -v "${pear_overlap}" -g "${pear_stattest}" -p "${pear_pvalue}" \
               > "${tmp_dir}/pear.log" 2>&1; then
            echo "skipped: PEAR failed:" >&2
            cat "${tmp_dir}/pear.log" >&2
            skipped+=("$i"); continue
        fi
        pear_un=$(( $(wc -l < "${tmp_dir}/pear.unassembled.forward.fastq") / 4 ))
        pear_disc=$(( $(wc -l < "${tmp_dir}/pear.discarded.fastq") / 4 ))
        paste -d '\t' - - - - < "${tmp_dir}/pear.assembled.fastq" | cut -f2,4 --output-delimiter=' ' > "${out}.qf"
    fi

    if ! "$qfscan" scan "$i" "$qf_ends" "$end_slop_value" "$q_floor" "$q_cutoff" "$cutoff_pct" \
           "$n_scanned" "$n_kept" "$pear_un" "$pear_disc" "$grid_csv" "$samples_csv" < "${out}.qf"; then
        echo "skipped: scan failed" >&2
        skipped+=("$i"); continue
    fi
    rm -f -- "${tmp_dir}"/sub* "${tmp_dir}"/pear.*

    echo "$(fmt_n "$n_kept") of $(fmt_n "$n_scanned") reads sampled ($(( $(date +%s) - t0 ))s)" >&2
done

if [ ! -s "$samples_csv" ]; then
    echo "Error: no samples could be scanned." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
header() {
    echo "QUALITY FILTER PREVIEW (estimates from a random subsample; not exact)"
    if [ "$merge" == "TRUE" ]; then
        echo "Mode: merge=TRUE (PEAR overlap=${pear_overlap:-}, test=${pear_stattest:-}, p=${pear_pvalue:-}), end_slop=${end_slop}"
    elif [ "$singleend" == "FALSE" ]; then
        echo "Mode: paired-end, R1 + revcomp(R2) (merge=FALSE), end_slop=${end_slop}"
    else
        echo "Mode: single-end, end_slop=${end_slop}"
    fi
    echo "Current settings: q_floor=${q_floor}, q_cutoff=${q_cutoff}, cutoff_pct=${cutoff_pct}"
    if [ "${#skipped[@]}" -gt 0 ]; then
        echo "Skipped samples (would also fail in the real run): ${skipped[*]}"
    fi
}

render_args=("$grid_csv" "$samples_csv" "$q_floor" "$q_cutoff" "$cutoff_pct" "$end_slop" "$merge")

{ header; "$qfscan" render "${render_args[@]}" 0 0 1; } > "$report_txt"

use_color=0
if [ "$color_opt" -eq 1 ] && [ -t 1 ] && [ -z "${NO_COLOR:-}" ] && [ "${TERM:-dumb}" != "dumb" ]; then
    use_color=1
fi
use_utf8=0
if [[ "$(locale charmap 2> /dev/null)" == "UTF-8" ]]; then use_utf8=1; fi

echo ""
header
"$qfscan" render "${render_args[@]}" "$use_color" "$use_utf8" "$show_all"
echo ""
echo "Saved: logs/filter_preview.txt (every sample), logs/filter_preview_grid.csv, logs/filter_preview_samples.csv"
