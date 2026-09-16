#!/bin/bash
set -uo pipefail
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
# Derive sample names from Fastq/ filenames (same logic used for preflight
# validation in process_ngs.sh). Looks one directory level down too, since a
# prior run may have already moved a sample's files into Fastq/<sample>/.
# -----------------------------------------------------------------------------
derive_sample_names() {
    local fastq_dir="$1"
    local f base
    local -a names=()
    shopt -s nullglob
    for f in "${fastq_dir}"/*_S*_L001_R1_001.fastq.gz \
             "${fastq_dir}"/*_S*_L001_R1_001.fastq \
             "${fastq_dir}"/*/*_S*_L001_R1_001.fastq.gz \
             "${fastq_dir}"/*/*_S*_L001_R1_001.fastq; do
        base=$(basename "$f")
        base="${base%_S*_L001_R1_001.fastq.gz}"
        base="${base%_S*_L001_R1_001.fastq}"
        names+=("$base")
    done
    shopt -u nullglob
    printf '%s\n' "${names[@]}" | awk 'NF && !seen[$0]++'
}

mapfile -t sample_names_array < <(derive_sample_names "${data_filepath}/Fastq")
if [ "${#sample_names_array[@]}" -eq 0 ]; then
    echo "Error: no fastq files found in ${data_filepath}/Fastq matching {sample}_S{#}_L001_R1_001.fastq[.gz]."
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

# Function to format reads. Accepts one or more input files -- a single-end
# sample can legitimately be split across several (e.g. multiple lanes), and
# passing a glob straight through as a single "$1" would silently drop every
# file after the first match.
format_reads() {
    echo -e "\n$(date '+%I:%M%p') -- FORMATTING READS"
    if ! cat -- "$@" | paste -d '\t' - - - - | awk -F '\t' '{print $2, $4}' > combined.fastq; then
        echo -e "$(date '+%I:%M%p') -- ERROR: failed to format reads"
        return 1
    fi
    echo -e "$(date '+%I:%M%p') -- READS FORMATTED"
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
    if compgen -G "../${i}_"* > /dev/null 2>&1; then
        if ! mv ../"${i}"_* .; then
            echo "ERROR: failed to move input files for sample ${i}, skipping."
            failed_samples+=("$i"); continue
        fi
    elif ! compgen -G "${i}_"* > /dev/null 2>&1; then
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
        if ! awk '{
        	if(NR%4==2){
            	seq=$0
            	gsub("A","t",seq); gsub("T","a",seq); gsub("G","c",seq); gsub("C","g",seq)
            	n = length(seq); rc=""
            	for(i=n;i>=1;i--){rc=rc toupper(substr(seq,i,1))}
            	print rc
        	} else {print}
   		}' *_R2_001.fastq > R2_rc.fastq; then
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
        format_reads "${i}"_* || { failed_samples+=("$i"); continue; }
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
