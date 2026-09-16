#!/bin/bash
set -euo pipefail

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <params_file.env>"
  exit 1
fi

params_file=$1

# Check if the params file exists
if [ ! -f "$params_file" ]; then
  echo "Error: Params file '$params_file' not found."
  exit 1
fi

# Get parameters from params file
source "$params_file"

# -----------------------------------------------------------------------------
# Derive sample names from Fastq/ filenames instead of a manually maintained
# params.env list. Expected naming convention (see README):
#   {samplename}_S{#}_L001_R1_001.fastq.gz  (and _R2_ for paired-end)
# Also looks one directory level down, since a prior run may have already
# moved a sample's files into Fastq/<sample>/ (see merge_reads.sh).
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

sample_names_array=()
if [[ -d "${data_filepath:-}/Fastq" ]]; then
    mapfile -t sample_names_array < <(derive_sample_names "${data_filepath}/Fastq")
fi

if [ "${#sample_names_array[@]}" -gt 0 ]; then
    detected_samples_str="${sample_names_array[*]}"
else
    detected_samples_str="(none detected)"
fi

# Print parameters to the terminal
echo "================ Parameter Settings ================"
echo "Data Filepath       : ${data_filepath:-}"
echo "Pear Filepath       : ${pear_filepath:-}"
echo "Pear Overlap        : ${pear_overlap:-}"
echo "Pear Stat Test      : ${pear_stattest:-}"
echo "Pear P-value        : ${pear_pvalue:-}"
echo "Merge               : ${merge:-}"
echo "Single-end          : ${singleend:-}"
echo "Compiler Filepath   : ${compiler_filepath:-}"
echo "CPUs                : ${cpus:-}"
echo "Memory              : ${memory:-}"
echo "Disk                : ${disk:-}"
echo "Q Floor             : ${q_floor:-}"
echo "Q Cutoff            : ${q_cutoff:-}"
echo "Cutoff Percent      : ${cutoff_pct:-}"
echo "Sample Names        : ${detected_samples_str} (auto-detected from Fastq/)"
echo "Reorganize          : ${reorganize:-}"
echo "Notify Email        : ${notify_email:-(none - notifications disabled)}"
echo "===================================================="

# -----------------------------------------------------------------------------
# PREFLIGHT VALIDATION
# Collect every problem before failing, so a bad config is reported all at
# once instead of one condor-queue-wasting attempt at a time.
# -----------------------------------------------------------------------------
preflight_errors=()

check_true_false() {
    local name="$1" value="$2"
    if [[ "$value" != "TRUE" && "$value" != "FALSE" ]]; then
        preflight_errors+=("Parameter '${name}' must be TRUE or FALSE (got '${value}').")
    fi
}

check_positive_int() {
    local name="$1" value="$2"
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        preflight_errors+=("Parameter '${name}' must be a positive integer (got '${value}').")
    fi
}

if [[ -z "${data_filepath:-}" ]]; then
    preflight_errors+=("data_filepath is not set.")
elif [[ ! -d "$data_filepath" ]]; then
    preflight_errors+=("data_filepath ('${data_filepath}') does not exist or is not a directory.")
elif [[ ! -d "${data_filepath}/Fastq" ]]; then
    preflight_errors+=("Fastq directory not found at '${data_filepath}/Fastq'.")
fi

check_true_false merge "${merge:-}"
check_true_false singleend "${singleend:-}"
check_true_false reorganize "${reorganize:-}"

if [[ "${merge:-}" == "TRUE" && "${singleend:-}" == "TRUE" ]]; then
    preflight_errors+=("merge=TRUE implies paired-end reads, but singleend is also TRUE. Set one of them to FALSE.")
fi

if [[ "${merge:-}" == "TRUE" ]]; then
    if [[ -z "${pear_filepath:-}" ]]; then
        preflight_errors+=("pear_filepath is not set (required when merge=TRUE).")
    elif [[ ! -x "$pear_filepath" ]]; then
        preflight_errors+=("pear_filepath ('${pear_filepath}') does not exist or is not executable.")
    fi
    check_positive_int pear_overlap "${pear_overlap:-}"
fi

if [[ -z "${compiler_filepath:-}" ]]; then
    preflight_errors+=("compiler_filepath is not set.")
elif [[ ! -x "$compiler_filepath" ]]; then
    preflight_errors+=("compiler_filepath ('${compiler_filepath}') does not exist or is not executable.")
fi

check_positive_int cpus "${cpus:-}"
check_positive_int q_floor "${q_floor:-}"
check_positive_int q_cutoff "${q_cutoff:-}"

if ! [[ "${cutoff_pct:-}" =~ ^[0-1](\.[0-9]+)?$ ]]; then
    preflight_errors+=("cutoff_pct should be a decimal between 0 and 1 (got '${cutoff_pct:-}').")
fi

if ! [[ "${memory:-}" =~ ^[0-9]+[A-Za-z]+$ ]]; then
    preflight_errors+=("memory should look like '2G' (got '${memory:-}').")
fi
if ! [[ "${disk:-}" =~ ^[0-9]+[A-Za-z]+$ ]]; then
    preflight_errors+=("disk should look like '60G' (got '${disk:-}').")
fi

if ! command -v condor_submit > /dev/null 2>&1; then
    preflight_errors+=("condor_submit was not found on PATH. Is HTCondor installed/loaded on this machine?")
fi

# notify_email is optional. If set, HTCondor emails that address on job
# failure; if left blank, notification is disabled entirely.
if [[ -n "${notify_email:-}" && ! "${notify_email}" =~ ^[^[:space:]@]+@[^[:space:]@]+\.[^[:space:]@]+$ ]]; then
    preflight_errors+=("notify_email ('${notify_email}') doesn't look like a valid email address. Leave it blank to disable notifications.")
fi

if [[ "${#sample_names_array[@]}" -eq 0 ]]; then
    preflight_errors+=("No sample fastq files found in '${data_filepath:-}/Fastq'. Expected names like {sample}_S{#}_L001_R1_001.fastq.gz.")
else
    for s in "${sample_names_array[@]}"; do
        # Check top-level Fastq/ (compressed or already-decompressed from a
        # previous run) and Fastq/<sample>/ (already moved by merge_reads.sh).
        if ! compgen -G "${data_filepath}/Fastq/${s}_S*_L001_R1_001.fastq*" > /dev/null 2>&1 \
           && ! compgen -G "${data_filepath}/Fastq/${s}/${s}_S*_L001_R1_001.fastq*" > /dev/null 2>&1; then
            preflight_errors+=("Sample '${s}': R1 fastq file not found.")
        fi
        if [[ "${merge:-}" == "TRUE" || "${singleend:-}" == "FALSE" ]]; then
            if ! compgen -G "${data_filepath}/Fastq/${s}_S*_L001_R2_001.fastq*" > /dev/null 2>&1 \
               && ! compgen -G "${data_filepath}/Fastq/${s}/${s}_S*_L001_R2_001.fastq*" > /dev/null 2>&1; then
                preflight_errors+=("Sample '${s}': paired-end run but R2 fastq file not found.")
            fi
        fi
    done
fi

if [[ "${#preflight_errors[@]}" -gt 0 ]]; then
    echo "================ Preflight Validation Failed ================" >&2
    for err in "${preflight_errors[@]}"; do
        echo " - $err" >&2
    done
    echo "===============================================================" >&2
    exit 1
fi

echo "Preflight validation passed."

# -----------------------------------------------------------------------------
# Decompress fastqs when the downstream steps need plain-text input.
# Skipped harmlessly if there's nothing to decompress (e.g. on a rerun).
# -----------------------------------------------------------------------------
if [[ "$merge" == "FALSE" || "$singleend" == "TRUE" ]]; then
    cd "${data_filepath}/Fastq"
    if compgen -G "*.fastq.gz" > /dev/null 2>&1; then
        # -f is required here, not just a nice-to-have: files pulled from the
        # GLBRC Data Catalog via `datasync files pull` land as symlinks, and
        # gunzip can refuse/misbehave on a symlinked source without -f.
        gunzip -f -- *.fastq.gz
    fi
fi

cd "${data_filepath}/pipeline"

# merge_reads.sh is the file HTCondor directly executes on the remote node
# (see submit_template.sub's `executable =` line). On a shared filesystem
# like Scarcity's, HTCondor runs it in place from disk rather than a
# transferred sandbox copy, so it must carry the execute bit *on disk* --
# git does not preserve/restore this reliably across clones, and a missing
# +x here surfaces as a cryptic "errno=13: Permission denied" from HTCondor
# with no other clue. Set it unconditionally rather than just checking.
if ! chmod +x merge_reads.sh; then
    echo "Error: could not set the execute bit on merge_reads.sh." >&2
    exit 1
fi

# Prepare C++ filtering script from template.
# NOTE: the template is intentionally kept in place (not deleted) so this
# script can be re-run without re-cloning the repo.
sed -e "s/{{Q_FLOOR}}/$q_floor/" \
    -e "s/{{Q_CUTOFF}}/$q_cutoff/" \
    -e "s/{{CUTOFF_PCT}}/$cutoff_pct/" \
    QF3_template.cpp > QF3.cpp

# Compile the quality filtering script (static compiling allows this to run through a cluster job)
if ! "${compiler_filepath}" -static-libstdc++ -o QF3.out QF3.cpp; then
    echo "Error: compilation of QF3.cpp failed. Check that compiler_filepath ('${compiler_filepath}') is a working g++ executable." >&2
    exit 1
fi

# notify_email is optional: if unset, disable HTCondor email notifications
# entirely but still give notify_user a harmless, non-empty value.
if [[ -n "${notify_email:-}" ]]; then
    notification_value="ERROR"
    notify_user_value="${notify_email}"
else
    notification_value="NEVER"
    notify_user_value="$(whoami)"
fi

# Prepare condor submit file from template.
# NOTE: the template is intentionally kept in place (not deleted) so this
# script can be re-run without re-cloning the repo.
sed -e "s|{{ARGS}}|$params_file|g" \
    -e "s|{{INPUT_FILES}}|${data_filepath}/${params_file}|g" \
    -e "s|{{INITIALDIR}}|${data_filepath}/pipeline|g" \
    -e "s|{{CPUS}}|${cpus}|g" \
    -e "s|{{MEMORY}}|${memory}|g" \
    -e "s|{{DISK}}|${disk}|g" \
    -e "s|{{NOTIFICATION}}|${notification_value}|g" \
    -e "s|{{NOTIFY_USER}}|${notify_user_value}|g" \
    submit_template.sub > process_ngs.sub

# Submit the job with condor_submit
if ! condor_submit process_ngs.sub; then
    echo "Error: condor_submit failed. Check process_ngs.sub and your HTCondor environment." >&2
    exit 1
fi
