# NGS Processing Pipeline - Version 11.0
**Updated:** 10/19/2026

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
Within your home directory, clone the pipeline:
```
git clone https://github.com/stmiller2/merge_and_filter_scarcity
```
and rename if desired:
```
mv merge_and_filter_scarcity/ my_experiment/
cd my_experiment
```

### STEP 3 — Place fastqs within the Fastq/ directory
Fastq filenames should look like this:
- `{samplename}_S{#}_L001_R1_001.fastq.gz` – Forward read  
- `{samplename}_S{#}_L001_R2_001.fastq.gz` – Reverse read

For single-end runs (`singleend=TRUE`), provide only one read per sample (R1 *or* R2, not both). R2 is reverse-complemented before filtering so output is in R1 orientation.

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
| perl_filepath     | filepath   | Path to perl executable - `/usr/bin/perl` (Only needed for reverse-complementing when merge=FALSE) |
| cpus              | int        | Number of cores requested – suggest 8 |
| memory            | int + K/M/G | Amount of RAM requested – suggest 4G. **Needs K/M/G suffix, not KB/MB/GB** |
| disk              | int + unit | Amount of disk space requested – suggest 60GB |
| q_floor           | int        | Reads with any bases below this value are discarded |
| q_cutoff          | int        | Reads with too many bases below this value are discarded (see cutoff_pct) |
| cutoff_pct        | float      | Proportion of bases that must be at or above q_cutoff for a read to be kept |
| end_slop          | int or FALSE | Number of bases at sequencing-read 5' ends exempted from quality filtering. Suggest FALSE unless necessary (inspect fastqs for early-cycle N's) |
| trim_slop         | TRUE/FALSE | Remove end_slop bases from output sequences - suggest TRUE if using end_slop. If FALSE, they are kept in the output but still not used for quality filtering |
| space_saver       | TRUE/FALSE | Deletes intermediate files (`work/<sample>/`) once run finishes and gzips `poor_reads` csvs |
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
tail -f logs/run_progress.log
```
(exit using ctrl+C). You can also check `condor_q` for overall job status (`RUN`/`IDLE`/`HELD`/etc.).

### Output
```
my_experiment/
├── Fastq/                       raw inputs
├── csvs/
│   ├── good_reads/              good_reads_<sample>.csv
│   ├── poor_reads/              poor_reads_<sample>.csv  (gzipped if space_saver)
│   └── info/                    info_<sample>.csv  (read counts at each step)
├── logs/                        run_progress.log, condor.{log,out,err}, params_used.env
└── work/<sample>/               intermediate files; removed if space_saver (unless sample failed)
```

`poor_reads` columns: sequence, failure code, value. `E1` = a base below `q_floor` (value = 0-based index of the first such base); `E2` = too few bases at/above `q_cutoff` (value = actual proportion); `E3` = read too short after removing `end_slop` (value = read length).

### STEP 7 — Transfer data off of scarcity
Home directory is **not backed up**, so once a run finishes, move anything you want to keep somewhere permanent, such as:
- Copy final results into `/mnt/bigdata/processed_data/YOUR_USERNAME/` (request a folder from `helpdesk@energy.wisc.edu` if you don't have one yet).
- Push results to the GLBRC Data Catalog: `datasync files push --dataset-id YOUR_DATASET_ID --path ./csvs`
- scp to your local device: **from your local device,** `scp -r USERNAME@scarcity-submit.glbrc.org:/home/GLBRCORG/username/my_experiment/ .` (copies entire run folder from Scarcity to your current directory)

### STEP 8
Publish CNS

---

## Acknowledgements
Initial scripts and framework are based on work from Anthony Meger, modified by me, and in many cases reviewed/updated by an LLM. I always carefully review LLM-generated code and take responsibility for the content here.

## Changelog

### Version 11 — 10/19/2026
- Reworked file organization esp. for large files:
  - `Fastq/` directory is not modified anymore, input files aren't moved into subdirectories
  - Intermediate files go in `work/<sample>/`; results go in `csvs/{good_reads,poor_reads,info}/`; logs go in `logs/`
  - `space_saver` param deletes each sample's `work/` dir on success and gzips `poor_reads`
- Reads are now streamed straight into the quality filtering script so `combined.fastq` is no longer written
- Removed `get_stats.sh`, read counting is now in `merge_reads.sh`
- New `perl_filepath` param (needed for fast R2 reverse-complementing when `merge=FALSE`)
- New quality-filtering params `end_slop` and `trim_slop`: exempt (and optionally trim) N bases at sequencing-read 5' ends from filtering
  - Useful for position-dependent N's in sequencing data due to sequencer focusing during first few cycles
- For `merge=FALSE` runs, fastqs are no longer decompressed, input is streamed with `zcat -f` instead.
- Updates related to single-end mode:
  - Fixed bug where R2 files weren't detected
  - R2 gets reverse-complemented in single-end mode
  - Preflight checks give an error if both R1 and R2 are provided when `singleend=TRUE`
- Improved filename matching so samples with shared prefix (e.g. `A` and `A_2`) don't collide

### Version 10 — 09/16/2026
- Ported the pipeline from BCC to **Scarcity**. Core logic is the same with the following changes:
  - For `merge=FALSE`, `process_ngs.sh` now decompresses with `gunzip -f` (otherwise this breaks with symlinks from data catalog)
  - `pipeline/submit_template.sub` has minor updates to fit Scarcity's options
  - Minor changes to filepath handling to fit Scarcity's directory structure
  - README updates for scarcity workflow
- Other quality of life improvements:
   - Faster reverse-complementing & concatenation when `merge=FALSE`
   - Sample names are now auto-detected from filenames in `Fastq/` and don't need to be set in `params.env`
   - Preflight check gives all errors at once before submitting to HTCondor
   - Added more informative error messages throughout all scripts
   - Made it easier to rerun the pipeline after error: `process_ngs.sh` doesn't delete templates `QF3_template.cpp`/`submit_template.sub` anymore and `get_stats.sh` overwrites `info.csv` instead of concatenating
   - Added a changelog :)