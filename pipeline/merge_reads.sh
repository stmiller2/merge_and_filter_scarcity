#!/bin/bash
set -uo pipefail
# extglob must be enabled before the functions that use +(...)/?(...) are parsed
shopt -s extglob
# NOTE: -e is intentionally NOT set here. Failures are checked explicitly at
# each step so that one bad sample can be logged and skipped (via `continue`)
# without aborting the whole batch run.

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
exec &> run_progress.log
mkdir -p "${data_filepath}/csvs/"
rm -f "${data_filepath}/Fastq/paste_fastq_files_here"

if [ "$reorganize" == "TRUE" ]; then
    mkdir -p "${data_filepath}/csvs/raw"/{good_reads,poor_reads,info,combined}
fi

# -----------------------------------------------------------------------------
# Derive sample names from Fastq/ filenames
# -----------------------------------------------------------------------------
derive_sample_names() {
    local fastq_dir="$1"
    local f base
    local -a names=()
    shopt -s nullglob
    # Match R1 or R2 so single-end runs on R2-only data are detected too.
    for f in "${fastq_dir}"/*_S[0-9]*_R[12]_001.fastq.gz \
             "${fastq_dir}"/*_S[0-9]*_R[12]_001.fastq \
             "${fastq_dir}"/*/*_S[0-9]*_R[12]_001.fastq.gz \
             "${fastq_dir}"/*/*_S[0-9]*_R[12]_001.fastq; do
        base=$(basename "$f")
        base="${base%_S[0-9]*_R[12]_001.fastq.gz}"
        base="${base%_S[0-9]*_R[12]_001.fastq}"
        names+=("$base")
    done
    shopt -u nullglob
    [[ "${#names[@]}" -gt 0 ]] || return 0
    printf '%s\n' "${names[@]}" | awk 'NF && !seen[$0]++'
}

mapfile -t sample_names_array < <(derive_sample_names "${data_filepath}/Fastq")
if [ "${#sample_names_array[@]}" -eq 0 ]; then
    echo "Error: no fastq files found in ${data_filepath}/Fastq matching {sample}_S{#}_R[12]_001.fastq[.gz] (with or without a _L001_-style lane segment)."
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

# Function to format reads. Accepts one or more input files (e.g. sample split across lanes)
format_reads() {
    echo -e "\n$(date '+%I:%M%p') -- FORMATTING READS"
    if ! cat -- "$@" | paste -d '\t' - - - - | awk -F '\t' '{print $2, $4}' > combined.fastq; then
        echo -e "$(date '+%I:%M%p') -- ERROR: failed to format reads"
        return 1
    fi
    echo -e "$(date '+%I:%M%p') -- READS FORMATTED"
}

# Collect sample $1's fastqs in directory $2 into the array sample_fastqs.
# Strict pattern: {sample}_S{digits}[_L{digits}]_R{1|2}_001.fastq[.gz], so
# sample 'A' does not pick up files belonging to sample 'A_2'.
find_sample_fastqs() {
    local s="$1" dir="$2"
    shopt -s nullglob
    sample_fastqs=( "${dir}/${s}"_S+([0-9])?(_L+([0-9]))_R[12]_001.fastq?(.gz) )
    shopt -u nullglob
}

# Reverse-complement one or more fastq files into $1 (quality strings reversed too).
revcomp_fastq() {
    local out="$1"; shift
    if ! cat -- "$@" | paste -d '\t' - - - - > rc.4col.tsv; then
        rm -f rc.4col.tsv
        return 1
    fi
    paste -d '\t' \
        <(cut -f1 rc.4col.tsv) \
        <(cut -f2 rc.4col.tsv | rev | tr 'ACGTNacgtn' 'TGCANTGCAN') \
        <(cut -f3 rc.4col.tsv) \
        <(cut -f4 rc.4col.tsv | rev) \
        | tr '\t' '\n' > "$out"
    local status=$?
    rm -f rc.4col.tsv
    return $status
}

failed_samples=()

# Iterate through samples
for i in "${sample_names_array[@]}"; do
	echo -e "\033[1m$(printf %80s |tr " " "=")\033[0m\n"
	echo -e "\033[1m$(printf %$(((80-(18+${#i}))/2))s |tr " " " ")PROCESSING SAMPLE $i\033[0m\n"
	echo -e "\033[1m$(printf %80s |tr " " "=")\033[0m"

    # Prepare sample directory
    mkdir -p "${data_filepath}/Fastq/${i}/"
    cd "${data_filepath}/Fastq/${i}/" || {
        echo "ERROR: cannot access Fastq/${i}/, skipping sample."
        failed_samples+=("$i"); continue
    }

    # Move this sample's input files in, unless a previous run already did so
    # (keeps the pipeline safely re-runnable).
    find_sample_fastqs "$i" ..
    if [ "${#sample_fastqs[@]}" -gt 0 ]; then
        if ! mv -- "${sample_fastqs[@]}" .; then
            echo "ERROR: failed to move input files for sample ${i}, skipping."
            failed_samples+=("$i"); continue
        fi
    fi
    find_sample_fastqs "$i" .
    if [ "${#sample_fastqs[@]}" -eq 0 ]; then
        echo "ERROR: no input files found for sample ${i} (checked Fastq/ and Fastq/${i}/), skipping."
        failed_samples+=("$i"); continue
    fi

    # Merge or concatenate reads
    if [ "$merge" == "TRUE" ]; then
        log start "MERGING PAIRED-END READS"
        pear_output=$("${pear_filepath}" -f *_R1_001.fastq.gz -r *_R2_001.fastq.gz -o combined \
            -y "${memory}" -j "${cpus}" -v "${pear_overlap}" -g "${pear_stattest}" -p "${pear_pvalue}" 2>&1)
        pear_status=$?
        echo "$pear_output" | grep -E -m 3 "Assembled reads|Discarded reads|Not assembled reads" | awk 'NF' | sed 's/^/           /'
        if [ $pear_status -ne 0 ]; then
            echo "ERROR: PEAR failed for sample ${i} (exit code ${pear_status}):"
            echo "$pear_output"
            failed_samples+=("$i"); continue
        fi
        log stop "READS MERGED"
        format_reads combined.assembled.fastq || { failed_samples+=("$i"); continue; }
    elif [ "$singleend" == "FALSE" ]; then
        log start "CREATING REVERSE COMPLEMENT OF READ 2"
        if ! revcomp_fastq R2_rc.fastq *_R2_001.fastq; then
            echo "ERROR: reverse-complement step failed for sample ${i}, skipping."
            failed_samples+=("$i"); continue
        fi
        log stop "REVERSE COMPLEMENTED"
        log start "CONCATENATING PAIRED-END READS"
        if ! paste -d "" *_R1_001.fastq R2_rc.fastq > combined.assembled.fastq; then
            echo "ERROR: concatenation step failed for sample ${i}, skipping."
            failed_samples+=("$i"); continue
        fi
        log stop "READS MERGED"
        format_reads combined.assembled.fastq || { failed_samples+=("$i"); continue; }
        rm -f R2_rc.fastq
    else
        # Single-end. Decompress any .gz left in Fastq/<sample>/ by an earlier
        # merge=TRUE run (process_ngs.sh only gunzips top-level Fastq/).
        gz_files=()
        for f in "${sample_fastqs[@]}"; do
            [[ "$f" == *.gz ]] && gz_files+=("$f")
        done
        if [ "${#gz_files[@]}" -gt 0 ] && ! gunzip -f -- "${gz_files[@]}"; then
            echo "ERROR: failed to decompress fastq files for sample ${i}, skipping."
            failed_samples+=("$i"); continue
        fi
        find_sample_fastqs "$i" .

        # Split by read (multiple lanes of the same read are concatenated).
        r1_files=(); r2_files=()
        for f in "${sample_fastqs[@]}"; do
            if [[ "$f" == *_R2_001.fastq ]]; then r2_files+=("$f"); else r1_files+=("$f"); fi
        done

        if [ "${#r1_files[@]}" -gt 0 ] && [ "${#r2_files[@]}" -gt 0 ]; then
            echo "ERROR: sample ${i} has both R1 and R2 fastq files in single-end mode, skipping."
            failed_samples+=("$i"); continue
        elif [ "${#r2_files[@]}" -gt 0 ]; then
            log start "SINGLE-END READ 2: CREATING REVERSE COMPLEMENT"
            if ! revcomp_fastq R2_rc.fastq "${r2_files[@]}"; then
                echo "ERROR: reverse-complement step failed for sample ${i}, skipping."
                failed_samples+=("$i"); continue
            fi
            log stop "REVERSE COMPLEMENTED"
            format_reads R2_rc.fastq || { failed_samples+=("$i"); continue; }
            rm -f R2_rc.fastq
        else
            format_reads "${r1_files[@]}" || { failed_samples+=("$i"); continue; }
        fi
    fi

    # Quality filtering
    if ! cp "${data_filepath}/pipeline/QF3.out" .; then
        echo "ERROR: could not copy QF3.out for sample ${i}, skipping."
        failed_samples+=("$i"); continue
    fi
    log start "FILTERING FOR HIGH QUALITY READS"
    if ! ./QF3.out; then
        echo "ERROR: quality filtering failed for sample ${i} (missing/empty combined.fastq?), skipping."
        rm -f QF3.out
        failed_samples+=("$i"); continue
    fi
    log stop "QUALITY FILTERING COMPLETE"
    rm -f QF3.out

    # Stats calculation
    if ! cp "${data_filepath}/pipeline/get_stats.sh" .; then
        echo "ERROR: could not copy get_stats.sh for sample ${i}, skipping."
        failed_samples+=("$i"); continue
    fi
    log start "GETTING ASSEMBLY AND FILTERING STATISTICS"
    chmod +x get_stats.sh
    if ! ./get_stats.sh; then
        echo "ERROR: stats calculation failed for sample ${i}, skipping."
        rm -f get_stats.sh
        failed_samples+=("$i"); continue
    fi
    log stop "READ FATES WRITTEN TO INFO.CSV"
    rm -f get_stats.sh

    # Reorganize files
    if [ "$reorganize" == "TRUE" ]; then
        log start "REORGANIZING FILES"
        if mv good_reads.csv "${data_filepath}/csvs/raw/good_reads/good_reads_${i}.csv" &&
           mv poor_reads.csv "${data_filepath}/csvs/raw/poor_reads/poor_reads_${i}.csv" &&
           mv info.csv "${data_filepath}/csvs/raw/info/info_${i}.csv" &&
           mv combined.fastq "${data_filepath}/csvs/raw/combined/combined_${i}.fastq"; then
            log stop "REORGANIZATION COMPLETE"
        else
            echo "ERROR: reorganizing output files failed for sample ${i}."
            failed_samples+=("$i"); continue
        fi
    else
        log start "NO REORGANIZATION REQUESTED. FIND GOOD_READS WITHIN FASTQ SUBDIRECTORIES."
    fi

    echo -e "\nSAMPLE $i COMPLETE\n\n"
done

# Finish
end_time=$(date --utc +%s)
runtime=$((end_time - start_time))
echo -e "\n$(printf %80s |tr " " "=")"
if [ "${#failed_samples[@]}" -gt 0 ]; then
    echo "Failed samples: ${failed_samples[*]}"
    log stop "Done WITH ERRORS. Total runtime: $(date -u -d "@${runtime}" +%H:%M:%S)"
    exit 1
else
    log stop "Done. All samples completed successfully. Total runtime: $(date -u -d "@${runtime}" +%H:%M:%S)"
fi