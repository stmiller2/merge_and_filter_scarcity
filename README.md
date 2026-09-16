# NGS Processing Pipeline - Version 13
**Updated:** 09/16/2026

**Contact:** stmiller2@wisc.edu

---

## Overview
This pipeline facilitates the merging of paired-end reads and quality filtering for Next-Generation Sequencing (NGS) data. It also supports processing single-end reads or paired-end reads without overlap. It runs as an HTCondor job on **Scarcity**, the WEI/GLBRC compute cluster. Follow the steps below to process your experiment efficiently.

---

## Usage Instructions

### STEP 1 — Get access to Scarcity
If you don't already have a Scarcity account, email the Computational Biology Lead, Kevin Myers (`kmyers2@wisc.edu`), with your estimated storage needs and a short description of your workflow. You'll be sent some onboarding material and then granted access.

### STEP 2 — Connect to Scarcity
If you're off the WEI/GLBRC network, connect to the VPN first (GlobalProtect, portal `weiglbrc.vpn.wisc.edu`, your GLBRC credentials). Then SSH into the submit node — **never run jobs directly on it**, it only has enough resources to host the login shell and hand jobs off to the pool:
```
ssh YOUR_USERNAME@scarcity-submit.glbrc.org
```
You'll land in your home directory, `/home/GLBRCORG/YOUR_USERNAME`. This is where you'll do all your work below. Note that this directory is **not backed up** — it's working space, not long-term storage (see Step 9).

### STEP 3 — Clone the pipeline
Do this directly on `scarcity-submit`, inside your home directory, and rename it to reflect your experiment:
```
git clone https://github.com/stmiller2/merge_and_filter_starterpack
mv merge_and_filter_starterpack/ my_experiment/
cd my_experiment
```

### STEP 4 — Get your sequencing data onto Scarcity
Filenames must follow this format (sample names are auto-detected from them — there is nothing to type into `params.env`):
- `{samplename}_S{#}_L001_R1_001.fastq.gz` – Forward read  
- `{samplename}_S{#}_L001_R2_001.fastq.gz` – Reverse read (paired-end only)

For any non-trivial amount of sequencing data, the preferred route is the **GLBRC Data Catalog** rather than copying raw files straight into your home directory:
1. If your data is from a GLBRC-funded JGI or UW-Biotech sequencing run, it's already synced into the Data Catalog automatically (nightly) — you don't need to upload anything.
2. Otherwise, push it there yourself with the `datasync` CLI (one-time setup: create an App Token at `data.glbrc.org` under your profile, then run `datasync` and paste it when prompted):
   ```
   datasync datasets create --name "My Experiment" --description "Short description" --path /path/to/your/fastqs
   ```
3. From inside this pipeline's `Fastq/` directory, pull the files down:
   ```
   cd Fastq
   datasync files pull --dataset-id YOUR_DATASET_ID --target .
   ```

**Important:** `datasync files pull` creates *symlinks* back into the Data Catalog's storage rather than physical copies — this is intentional (it keeps large raw fastqs out of your unbacked-up home directory), and this pipeline already handles it correctly (`process_ngs.sh` decompresses with `gunzip -f`, which is required for symlinked sources). You don't need to do anything special beyond pulling the files into `Fastq/`.

If you only have a handful of small files, plain SFTP/SCP into `Fastq/` also works fine (WinSCP or FileZilla for a GUI, or `scp -r your_dir YOUR_USERNAME@scarcity-submit.glbrc.org:~/my_experiment/Fastq/` from the command line).

### STEP 5 — Set up PEAR and a compiler
PEAR and a C++ compiler are **not** part of Scarcity's default cluster-wide software (check `/opt/bifxapps` if you want to confirm). The recommended approach on Scarcity is a personal Conda environment:
```
conda create -n ngs_pipeline -c bioconda -c conda-forge pear gxx_linux-64
conda activate ngs_pipeline
which pear
which g++
```
Use the paths these print in `params.env` below (they typically live at `~/.conda/envs/ngs_pipeline/bin/pear` and `~/.conda/envs/ngs_pipeline/bin/g++` or similar).

### STEP 6 — Set parameters
Edit `params.env` with your preferred editor:
```
vi params.env
```

| PARAMETER         | VALUE      | DESCRIPTION |
|-------------------|------------|-------------|
| data_filepath     | filepath   | Absolute path to this pipeline's directory, e.g. `/home/GLBRCORG/YOUR_USERNAME/my_experiment` |
| pear_filepath     | filepath   | Path to the PEAR executable — from your Conda env (Step 5) |
| pear_overlap      | int        | Minimum overlap for merging (see PEAR docs) – suggest 10 |
| pear_stattest     | int        | PEAR merging statistical test (see PEAR docs) – suggest 1 |
| pear_pvalue       | float      | P-value for PEAR statistical test. Set to 1.0 to disable the test |
| merge             | TRUE/FALSE | Typically TRUE – merges reads with PEAR. If FALSE, PEAR parameters are ignored and R1 is concatenated to the reverse complement of R2 |
| singleend         | TRUE/FALSE | FALSE for merging paired-end reads; TRUE for formatting & filtering single-end reads |
| compiler_filepath | filepath   | Path to g++ executable (C++ compiler) — from your Conda env (Step 5) |
| cpus              | int        | Number of cores requested – suggest 8 |
| memory            | int + unit | Amount of RAM requested – suggest 4GB |
| disk              | int + unit | Amount of disk space requested – suggest 60GB |
| q_floor           | int        | Reads with any bases below this value are discarded |
| q_cutoff          | int        | Reads with too many bases below this value are discarded (see cutoff_pct) |
| cutoff_pct        | float      | Proportion of bases that must be at or above q_cutoff for a read to be kept |
| reorganize        | TRUE/FALSE | If TRUE, all final `good_reads.csv` files (and others) will be reorganized into a single directory |
| notify_email      | email      | Optional. If set, HTCondor emails this address if the job errors out. Leave blank to disable. |

### STEP 7 — Start the processing run
```
chmod +x process_ngs.sh
./process_ngs.sh params.env
```
Params (and preflight validation results) will print to the console — double check everything is set correctly before it submits to the queue.

### STEP 8 — Monitor the run
```
tail -f -n 50 run_progress.log
```
(exit using ctrl+C). You can also check `condor_q` for overall job status (`RUN`/`IDLE`/`HELD`/etc.).

### STEP 9 — Preserve your results
Your home directory is **not backed up**. Once a run finishes, move anything you want to keep out of it:
- For long-term, backed-up storage: copy final results into `/mnt/bigdata/processed_data/YOUR_USERNAME/` (request a folder from `helpdesk@energy.wisc.edu` if you don't have one yet).
- To share results back through the Data Catalog: `datasync files push --dataset-id YOUR_DATASET_ID --path ./csvs`

### STEP 10
Publish CNS

---

## Changelog

### Version 13 — 09/16/2026
- Fixed the next Scarcity submission failure: after the executable-permission fix, jobs failed with `merge_reads.sh: line 15: params.env: No such file or directory` followed by `data_filepath: unbound variable`. The condor submit file's `arguments` was set to the bare filename `params.env`, which is resolved relative to `initialdir` (`.../pipeline/`) — one directory above where `params.env` actually lives. This was silently fine on CHTC, where HTCondor transfers the job into a private sandbox and both files land together, but Scarcity runs jobs in place on its shared filesystem, so the bare filename pointed at the wrong directory. `process_ngs.sh` now passes `params.env`'s full absolute path as the job argument instead, which is correct regardless of whether the execute model transfers files or runs in place.

### Version 12 — 09/16/2026
- Fixed a job-submission failure on Scarcity: HTCondor reported `Failed to execute '.../pipeline/merge_reads.sh' ... errno=13: Permission denied`. `merge_reads.sh` is the file HTCondor directly runs on the remote node (`executable =` in `submit_template.sub`), and it was tracked without the execute bit set — nothing in the pipeline or the old README ever `chmod +x`'d it (only `process_ngs.sh` got that treatment, but that's not the file Condor executes remotely). On a shared filesystem like Scarcity's, HTCondor runs the executable in place from disk rather than a transferred sandbox copy, so the on-disk permission bit matters and git alone doesn't guarantee it survives a clone/copy. `process_ngs.sh` now unconditionally `chmod +x`'s `merge_reads.sh` right before compiling/submitting, so this can't recur regardless of how the repo got onto the cluster.

### Version 11 — 09/16/2026
- Ported the pipeline (and this README) from CHTC to **Scarcity**, the WEI/GLBRC HTCondor cluster. No changes were needed to the core job logic — both are plain HTCondor pools with a shared home-directory filesystem — but the following were updated:
  - `params.env` defaults now point at a Scarcity home directory (`/home/GLBRCORG/...`) instead of a CHTC scratch path, and at a personal Conda environment for `pear_filepath`/`compiler_filepath`, since PEAR and a C++ compiler aren't part of Scarcity's default cluster-wide software.
  - `process_ngs.sh`'s decompression step now runs `gunzip -f` instead of plain `gunzip`. This is required, not cosmetic: files pulled from the GLBRC Data Catalog via `datasync files pull` land as symlinks, and Scarcity's own documentation notes that expanding those in place can fail without `-f`.
  - `pipeline/submit_template.sub` now sets `output = condor.out` (previously only `log`/`error` were set) and supports optional email notification on failure via a new `notify_email` parameter in `params.env` (blank by default — leaves notifications disabled).
  - README rewritten around the Scarcity workflow: requesting access, connecting (VPN + `scarcity-submit.glbrc.org`), pulling sequencing data from the Data Catalog via `datasync` (symlinks, not copies), setting up PEAR/g++ via Conda, and moving final results to backed-up storage (`/mnt/bigdata/processed_data/`) since the home directory isn't backed up.
  - `pipeline/submit_template.sub` now also gets its `output`/`notification`/`notify_user` filled in (see above).
- Fixed a latent bug found while testing the above: preflight's per-sample fastq check only recognized a compressed `.fastq.gz` file sitting directly in `Fastq/`, not an already-decompressed `.fastq` file there (it only recognized a decompressed file one directory down, in `Fastq/<sample>/`). A rerun of `process_ngs.sh` after `gunzip` had already run once would wrongly report "file not found" for samples whose files were sitting right there. Both locations now accept either extension.

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
