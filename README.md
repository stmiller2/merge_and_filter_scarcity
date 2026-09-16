# NGS Processing Pipeline - Version 10
**Updated:** 09/16/2026

**Contact:** stmiller2@wisc.edu

---

## Overview
This pipeline facilitates the merging of paired-end reads and quality filtering for Next-Generation Sequencing (NGS) data. It also supports processing single-end reads or paired-end reads without overlap. Follow the steps below to process your experiment efficiently.

---

## Usage Instructions

### STEP 1
Clone the `merge_and_filter_starterpack` GitHub repo and rename it to reflect your experiment.
```
git clone https://github.com/stmiller2/merge_and_filter_starterpack
mv merge_and_filter_starterpack/ my_experiment/
```

### STEP 2
Move all your fastq files from the sequencer to the `Fastq/` directory. Filenames should follow this format:
- `{samplename}_S{#}_L001_R1_001.fastq.gz` – Forward read  
- `{samplename}_S{#}_L001_R2_001.fastq.gz` – Reverse read (paired-end only)

Sample names are auto-detected from these filenames (there is nothing to type into `params.env`), so it's important that they match this pattern exactly.

### STEP 3
Set parameters in the `params.env` file with your preferred command line file editor:
```
vi params.env
```

| PARAMETER         | VALUE      | DESCRIPTION |
|-------------------|------------|-------------|
| data_filepath     | filepath   | Path to the file you are running this pipeline in |
| pear_filepath     | filepath   | Path to the PEAR executable |
| pear_overlap      | int        | Minimum overlap for merging (see PEAR docs) – suggest 10 |
| pear_stattest     | int        | PEAR merging statistical test (see PEAR docs) – suggest 1 |
| pear_pvalue       | float      | P-value for PEAR statistical test. Set to 1.0 to disable the test |
| merge             | TRUE/FALSE | Typically TRUE – merges reads with PEAR. If FALSE, PEAR parameters are ignored and R1 is concatenated to the reverse complement of R2 |
| singleend         | TRUE/FALSE | FALSE for merging paired-end reads; TRUE for formatting & filtering single-end reads |
| compiler_filepath | filepath   | Path to g++ executable (C++ compiler) |
| cpus              | int        | Number of cores requested – suggest 8 |
| memory            | int + unit | Amount of RAM requested – suggest 2G |
| disk              | int + unit | Amount of disk space requested – suggest 60G |
| q_floor           | int        | Reads with any bases below this value are discarded |
| q_cutoff          | int        | Reads with too many bases below this value are discarded (see cutoff_pct) |
| cutoff_pct        | float      | Proportion of bases that must be at or above q_cutoff for a read to be kept |
| reorganize        | TRUE/FALSE | If TRUE, all final `good_reads.csv` files (and others) will be reorganized into a single directory |

### STEP 4
Move the entire directory to your scratch folder:
```
cp -r {directory name} /scratch/{wiscID}
```

### STEP 5
Start the processing run

If you cloned the GitHub repository, you will first need to make process_ngs.sh executable with `chmod +x`.
```
chmod +x process_ngs.sh
./process_ngs.sh params.env
```

Params will print to the console -- double check everything is set correctly. 

### STEP 6
Monitor the run progress
```
tail -f -n 50 run_progress.log
```
(exit using ctrl+C)

### STEP 7
When finished, move the entire directory back to the fileserver.

### STEP 8
Publish CNS

---

## Changelog

### Version 10 — 09/16/2026
- Sample names are now auto-detected from filenames in `Fastq/`; the `sample_names` parameter has been removed from `params.env`.
- Added preflight validation to `process_ngs.sh`: it now checks that paths and executables exist, that parameter values are well-formed, and that every detected sample has the expected fastq files, reporting all problems at once before compiling anything or touching the condor queue.
- Added error handling throughout `process_ngs.sh` and `merge_reads.sh` (PEAR, compilation, `condor_submit`, file moves, quality filtering, stats) so a failed step is reported clearly instead of silently producing empty or incorrect output. A sample that fails is logged and skipped rather than silently corrupting output or (previously) not being distinguishable from a successful run.
- Fixed a bug where a single-end sample split across multiple files (e.g. multiple sequencing lanes) would only have its first file processed; all matching files are now concatenated before formatting.
- The pipeline is now safely re-runnable: `process_ngs.sh` no longer deletes `QF3_template.cpp`/`submit_template.sub` after use, and `merge_reads.sh` no longer fails if a sample's fastq files were already moved into its `Fastq/<sample>/` subdirectory by a previous run.
- `get_stats.sh` now overwrites `info.csv` instead of appending to it, so reruns don't produce duplicated stats rows.
- `QF3.out` now reports an error instead of silently writing empty output files if `combined.fastq` is missing or contains no reads.
- Quoted `compiler_filepath` in the compile step so paths containing spaces work correctly.
- `.gitignore` now excludes runtime-generated files (`run_progress.log`, `process_ngs.sub`, `QF3.cpp`, condor logs) and the `Fastq/`/`csvs/` data directories.
