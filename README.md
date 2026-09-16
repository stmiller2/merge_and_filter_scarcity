# NGS Processing Pipeline - Version 10
**Updated:** 09/16/2026

**Contact:** stmiller2@wisc.edu

## Overview
This pipeline facilitates the merging of paired-end reads and quality filtering for Next-Generation Sequencing (NGS) data. It also supports processing single-end reads or paired-end reads without overlap. It runs as an HTCondor job on **Scarcity**, the WEI/GLBRC compute cluster.

## Usage Instructions

### STEP 1 — Connect to Scarcity
Connect to the WEI/GLBRC VPN (GlobalProtect, portal `weiglbrc.vpn.wisc.edu`, login with GLBRC credentials). Then SSH into the submit node:
```
ssh YOUR_USERNAME@scarcity-submit.glbrc.org
```
You'll land in your home directory, `/home/GLBRCORG/YOUR_USERNAME`.

### STEP 2 — Clone the pipeline
Within your home directory:
```
git clone https://github.com/stmiller2/merge_and_filter_starterpack
mv merge_and_filter_starterpack/ my_experiment/
cd my_experiment
```

### STEP 3 — Place fastqs within the Fastq/ directory
Fastq filenames should look like this:
- `{samplename}_S{#}_L001_R1_001.fastq.gz` – Forward read  
- `{samplename}_S{#}_L001_R2_001.fastq.gz` – Reverse read (paired-end only)

Transfer your files to Scarcity with one of the following:
1. If your data is from a UWBC core sequencing run, it's already synced into the Data Catalog automatically (nightly). Pull the Fastqs from data catalog into the Fastq/ directory with:
   ```
   cd Fastq
   datasync files pull --dataset-id YOUR_DATASET_ID --target .
   ```
2. If your data is from an in-house sequencing run, navigate to the location of the files **on your local computer/fileserver** and transfer to scarcity with scp:
   ```
   scp sample1_S1_L001_R1_001.fastq.gz USERNAME@scarcity-submit.glbrc.org:/home/GLBRCORG/username/my_experiment/Fastq/
   ```
   As with cp, to secure copy a directory, use scp -r.

### STEP 4 — Set parameters
Edit `params.env` with your preferred editor, e.g.:
```
vi params.env
```

| PARAMETER         | VALUE      | DESCRIPTION |
|-------------------|------------|-------------|
| data_filepath     | filepath   | Absolute path to this pipeline's directory, e.g. `/home/GLBRCORG/YOUR_USERNAME/my_experiment` |
| pear_filepath     | filepath   | Path to the PEAR executable - `/mnt/cephfs/bifxapps/bin/pear` |
| pear_overlap      | int        | Minimum overlap for merging (see PEAR docs) – suggest 10 |
| pear_stattest     | int        | PEAR merging statistical test (see PEAR docs) – suggest 1 |
| pear_pvalue       | float      | P-value for PEAR statistical test. Set to 1.0 to disable the test |
| merge             | TRUE/FALSE | TRUE for overlapping paired-end reads – merges reads with PEAR. If FALSE, PEAR parameters are ignored and R1 is concatenated to the reverse complement of R2 |
| singleend         | TRUE/FALSE | FALSE for merging paired-end reads; TRUE for formatting & filtering single-end reads |
| compiler_filepath | filepath   | Path to g++ executable (C++ compiler) - `/usr/bin/g++` |
| cpus              | int        | Number of cores requested – suggest 8 |
| memory            | int + K/M/G | Amount of RAM requested – suggest 4G. **Needs K/M/G suffix, not KB/MB/GB** |
| disk              | int + unit | Amount of disk space requested – suggest 60GB |
| q_floor           | int        | Reads with any bases below this value are discarded |
| q_cutoff          | int        | Reads with too many bases below this value are discarded (see cutoff_pct) |
| cutoff_pct        | float      | Proportion of bases that must be at or above q_cutoff for a read to be kept |
| reorganize        | TRUE/FALSE | If TRUE, all final `good_reads.csv` files (and others) will be reorganized into a single directory (recommended) |
| notify_email      | email      | Optional. If set, HTCondor emails this address if the job errors out. Leave blank to disable. |

### STEP 5 — Start processing 
```
chmod +x process_ngs.sh
./process_ngs.sh params.env
```
Params will print to the console — double check everything is set correctly before it submits to the queue.

### STEP 6 — Monitor the run
View a live-updating tail of the run progress with:
```
tail -f -n 50 run_progress.log
```
(exit using ctrl+C). You can also check `condor_q` for overall job status (`RUN`/`IDLE`/`HELD`/etc.).

### STEP 7 — Transfer data off of scarcity
Home directory is **not backed up**, so once a run finishes, move anything you want to keep somewhere permanent, such as:
- Copy final results into `/mnt/bigdata/processed_data/YOUR_USERNAME/` (request a folder from `helpdesk@energy.wisc.edu` if you don't have one yet).
- Push results to the GLBRC Data Catalog: `datasync files push --dataset-id YOUR_DATASET_ID --path ./csvs`
- scp to your local device: **from your local device,** `scp -r USERNAME@scarcity-submit.glbrc.org:/home/GLBRCORG/username/my_experiment/csvs/ .` (copies `csvs` folder from Scarcity to your current directory)

### STEP 8
Publish CNS

---

## Acknowledgements
Initial scripts and framework are based on work from Anthony Meger, modified by me, and in many cases reviewed/updated by an LLM. I always carefully review LLM-generated code and take responsibility for the content here.

## Changelog

### Version 10 — 09/16/2026
- Sped up the `merge=FALSE` (reverse-complement + concatenate) path in `merge_reads.sh`. 
- Fixed bug in `merge=FALSE` where the R2 quality string wasn't reversed to match reverse-complemented reads. This didn't affect which reads passed or failed filtering (that only depends on the *set* of quality values present, not their order), but it did make the reported failure position in `poor_reads.csv`'s "E1" column wrong for reads that failed via the q_floor check.
- Ported the pipeline from BCC to **Scarcity**. No changes needed to core logic but the following were updated:
  - For `merge=FALSE`, `process_ngs.sh`'s decompression step now uses `gunzip -f` instead of plain `gunzip` (otherwise this breaks with symlinks from data catalog).
  - `pipeline/submit_template.sub` has minor updates to fit Scarcity's options
  - README updates for scarcity workflow
- Sample names are now auto-detected from filenames in `Fastq/` and don't need to be set in `params.env`
- Added preflight validation to `process_ngs.sh` to check that paths and executables exist, parameters are valid, and that every detected sample has the expected fastq files. Reports all problems at once before submitting to HTCondor.
- Added more informative error messages throughout all scripts.
- Made it easier to rerun the pipeline after error: `process_ngs.sh` doesn't delete templates `QF3_template.cpp`/`submit_template.sub` anymore and `get_stats.sh` overwrites `info.csv` instead of concatenating.
- Added a changelog :)
