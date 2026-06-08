# nf-gwascleaner - A Nextflow GWAS summary statistics cleaner pipeline (v1.0.0-prototype)

An automated, low-footprint Nextflow pipeline designed to ingest, quality-control, standardise, liftover, and coordinate-sort heterogeneous GWAS summary statistics into an analysis-ready format.

Last updated: 2026-06-08
Current version: **v1 prototype (Minimum Viable Product)**. 

---

## Overview

Genome-wide association analyses (GWAS) summary statistics are one of the workhorses of modern medical genetics. These flat text files contain results for marginal by-variant association tests, with multiple downstream research applications, such as polygenic risk scores, LD-score regression, fine-mapping, and functional genomics integration. GWAS summary statistics come in an array of column formats, coordinate systems, and header nomenclatures, and manually pre-processing them is time-consuming and error-prone.

This pipeline provides a unified, automated, and streaming architecture that dynamically assesses raw GWAS files in multiple genome builds and in different formats and produces standardised and anlysis-ready **hg38** datasets, setting them up for modern genomic databases. This is the Nextflow DSL2 version (and future expansion track) of the bash/R GWAS tools processing pipeline ([GWAS_tools](https://github.com/GRealesM/GWAS_tools)).

---

## Pipeline Roadmap & Scope

`nf-gwascleaner` is purpose-built for basic structural GWAS summary statistic harmonisation.

### What It Does:
* Renames heterogeneous headers to unified terminology 
    - `CHR` = Chromosome, 
    - `BP`  = Position (base pairs)
    - `SNPID` = Variant ID (eg. rsid), 
    - `REF` = Reference allele, 
    - `ALT` = Effect allele (to which beta/OR refer to),
    - `BETA` = Effect size (log(OR)), 
    - `SE` = Standard error, 
    - `P` = P-value
* Check minimum information exist to derive the above columns if any missing (see below).
* Checks for missing data in key columns (`BETA`/`OR`, and `SE`) and stops if above a threshold (Default: 50%). 
* Splits unified `CPRA` (Chrom:Pos:Ref:Alt) strings into separate core columns if structural coordinate columns are missing.
* Reconstructs missing `SNPID` headers from clean `CHR:BP` strings when variant IDs are missing.
* Derives missing `BETA` columns dynamically from Odds Ratios (OR).
* Computes missing Standard Errors (`SE`) from `BETA` and `P-value` using an inverse normal distribution survival function.
* Computes missing `P-values` backwards from `BETA` and `SE`.
* Automatically upgrades `hg18` and `hg19` builds to `hg38` coordinates.
* Generates persistent `.tsv` quality reports listing files that failed QC thresholds or tracking variant drop rates during liftOver.

> ⚠️ **A Crucial Note on REF/ALT Assignment**
> While every effort has been made to capture all historical variations of effect (to which the beta/OR refers) and reference alleles encoding, edge cases exist. For instance, nomenclatures like `A1/A2` are notorious: `A1` represents the effect allele (ALT) in some consortia but denotes the reference allele (REF) in others. 
> 
> To ensure absolute analytical safety, **users are strongly encouraged to double-check their allele columns.** If in doubt, manually rename your raw columns to `REF` (for reference) and `ALT` (for effect allele) before running the pipeline. Check `bin/inspect_gwas.py` for how `nf-gwascleaner` interprets each column name.


### What It Cannot Do (yet):
* **Allele Alignment / Strand Flipping:** It does not align alleles against a reference genome sequence (e.g., checking for triallelic sites, matching to the positive strand, or resolving ambiguous A/T or C/G SNPs).
* **Imputation / Fine-Mapping:** It does not infer missing genotypes or calculate linkage disequilibrium.
* **Sample Size / Frequency Imputation:** It passes through available allele frequencies (`ALT_FREQ`) and sample sizes (`N`) but does not infer them if missing.


---

## Prerequisites & Dependencies

Software infrastructure is natively managed by Nextflow. You do not need to manually install dependencies like scipy or liftOver on your host environment.

- Nextflow (DSL2 compatible)
- Conda or Mamba (for automated process containerisation)

Nextflow will automatically isolate process execution profiles, dynamically download the required version of ucsc-liftover from bioconda, and isolate Python library environments seamlessly behind the scenes.

## Installation & Usage

Clone the repository:

``` bash
git clone git@github.com:GRealesM/nf-gwascleaner.git
cd nf-gwascleaner
``` 
Configure your inputs in nextflow.config or supply them directly as command-line arguments:

``` bash
nextflow run main.nf \
  --gwas_dir /path/to/raw_gwas_files/ \
  --outdir /path/to/results/
```

To resume a failed run or process additional files without re-computing past steps, you can harness Nextflow's caching engine:

```bash
nextflow run main.nf -resume
```

**System capability notice**
This pipeline relies on native Unix process execution architectures and optimized Linux binaries (such as GNU `sort`). 

* **macOS / Linux:** Fully supported out of the box.
* **Windows 10/11:** **Not supported natively** via CMD or PowerShell. Windows users must run this pipeline inside **Windows Subsystem for Linux (WSL)** (Ubuntu 22.04 LTS or later recommended) with Conda/Mamba configured within the WSL environment.

### 🐍 Running Scripts Natively (Optional)

Software environments are fully automated by Nextflow during execution. However, if you wish to run or test the underlying Python core utilities inside `bin/` directly from your command line, ensure your local environment utilizes a matching footprint:

```bash
# Explicitly matching the prototype development environment
conda create -n gwas-clean-local python=3.14 scipy=1.17.1 -c conda-forge
conda activate gwas-clean-local
```

## Quality & Operational Reports

The pipeline isolates corrupted inputs using automated routing boundaries. If a file displays over 50% missing data in critical parameters or completely lacks the minimum attributes required for standardisation, the file is routed into a graceful soft-fail state.

Centralized audit reports are dynamically written to the output workspace:

`results/reports/pipeline_qc_failures.tsv`: Lists specific input files flagged with structural errors alongside descriptive reason strings.

`results/reports/liftover_drop_metrics.tsv`: Log summary reporting the exact count of coordinates that failed mapping parameters across individual datasets.

Finished outputs are compiled uniformly in `results/Final_GWAS/`.

## Roadmap (Future Version Enhancements)

[ ] Add extra liftover options for specific target builds (eg. hg19, T2T)

[ ] Add tools to address further formatting issues that GWAS summary statistics files may come with.

[ ] Add option to set custom NA threshold limit.

[ ] Add option for linear to log(OR) scaling in case-control datasets (will require N0 and N1 numbers to be supplied).

[ ] Add option adjust sdY for quantitative datasets, such that variance = 1 (useful for integration of multiple quantitative GWAS summary statistics).

[ ] Integrate native Docker/Singularity container configurations alongside Conda profiles for enterprise cloud scaling.

[ ] Expanded QC reports to check files beyond formatting.
