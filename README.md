# nf-gwascleaner
A tool to process GWAS summary statistics into a single format using Nextflow

# GWAS Clean Pipeline (v1.0.0-prototype)

An automated, low-footprint Nextflow pipeline designed to ingest, quality-control, standardise, liftover, and coordinate-sort heterogeneous GWAS summary statistics into an analysis-ready format.

---

## 🧬 Overview

Genetic association data from repositories like the GWAS Catalog come in an array of column formats, coordinate systems, and header nomenclatures. Pre-processing these files manually for downstream analyses (like Polygenic Risk Score generation, LD-Score Regression, or fine-mapping) is time-consuming and error-prone.

This pipeline provides a unified, automated, and streaming architecture that dynamically assesses raw GWAS files in multiple genome builds and produces standardised **hg38** datasets.

---

## 🚀 Key Features


* **Single-Pass Streaming ETL:** Modifying columns, computing metrics, and rewriting massive (10M+ row) data files in separate steps strains standard hardware I/O. The `STANDARDISE` process handles string manipulation, renaming, formatting, and mathematical conversions in an $O(N)$ streaming loop with an absolute minimal memory footprint.
* **Smart Parameter Inspection:** The pipeline runs an initial lightweight `INSPECTOR` process that checks file delimiters, maps heterogeneous headers using a built-in alias dictionary, tracks missing data ratios across a 500,000-line lookup, and infers the active genome build via coordinate overlaps against a built-in manifest.
* **Robust LiftOver Integration:** Standard implementations of UCSC `liftOver` can drop or scramble variants if duplicated RSIDs or structural variations are encountered in the 4th column of a standard BED file. This pipeline completely side-steps this issue by dynamically using internal row-indices as ephemeral BED identifiers, translating the coordinates, and seamlessly re-merging them back into the data stream.
* **Memory-Safe Coordinate Sorting:** Sorting millions of genomic coordinates in memory (e.g., via standard Python or R dataframes) can easily crash standard laptops or limited HPC nodes. This pipeline offloads sorting to native Linux `sort` optimized using an external merge-sort routine (`LC_ALL=C`), placing chromosomes (`1-22`, `X`, `Y`, `MT`) and numerical base-pairs in flawless physical sequence with zero memory overhead.

---

## 🛠️ Pipeline Roadmap & Scope

This repository represents a functional **v1 prototype (Minimum Viable Product)**. It is purpose-built for basic structural harmonisation but deliberately leaves complex biological filtering to downstream workflows.

### What It Does:
* Renames heterogeneous headers to unified terminology (`CHR`, `BP`, `SNPID`, `REF`, `ALT`, `BETA`, `SE`, `P`).
* Splits unified `CPRA` (Chrom:Pos:Ref:Alt) strings into separate core columns if structural coordinate columns are missing.
* Reconstructs missing `SNPID` headers from clean `CHR:BP` strings when variant IDs are missing.
* Derives missing `BETA` columns dynamically from Odds Ratios (`log(OR)`).
* Computes missing Standard Errors (`SE`) from `BETA` and `P-value` using an inverse normal distribution survival function.
* Computes missing `P-values` backwards from `BETA` and `SE`.
* Automatically upgrades `hg18` and `hg19` builds to `hg38` coordinates.
* Generates persistent `.tsv` quality reports listing files that failed QC thresholds or tracking variant drop rates during liftOver.

### What It Cannot Do (Out of Scope for v1):
* **Allele Alignment / Strand Flipping:** It does not align alleles against a reference genome sequence (e.g., checking for triallelic sites, matching to the positive strand, or resolving ambiguous A/T or C/G SNPs).
* **Imputation / Fine-Mapping:** It does not infer missing genotypes or calculate linkage disequilibrium.
* **Sample Size / Frequency Imputation:** It passes through available allele frequencies (`ALT_FREQ`) and sample sizes (`N`) but does not infer them if missing.

---

## 📊 Directory Structure

```text
gwas-clean-pipeline/
├── main.nf                 # Main Nextflow pipeline orchestrator
├── nextflow.config         # Runtime configuration & Conda environment controls
├── assets/
│   ├── Manifest_build_translator.tsv  # Reference list for genome build sniffing
│   ├── hg18ToHg38.over.chain.gz       # UCSC LiftOver chain file
│   └── hg19ToHg38.over.chain.gz       # UCSC LiftOver chain file
├── bin/                    # Scripts automatically injected into the environment $PATH
│   ├── inspect_gwas.py
│   ├── standardise_gwas.py
│   └── liftover_gwas.py
└── .gitignore
```

## ⚙️ Prerequisites & Dependencies
Software infrastructure is natively managed by Nextflow. You do not need to manually install dependencies like scipy or liftOver on your host environment.

- Nextflow (DSL2 compatible)
- Conda or Mamba (for automated process containerisation)

Nextflow will automatically isolate process execution profiles, dynamically download the required version of ucsc-liftover from bioconda, and isolate Python library environments seamlessly behind the scenes.

## 💻 Usage

Clone the repository:

``` bash
git clone git@github.com:your_username/gwas-clean-pipeline.git
cd gwas-clean-pipeline
``` 
Configure your inputs in nextflow.config or supply them directly as command-line arguments:

``` bash
nextflow run main.nf \
  --gwas_dir /path/to/raw_gwas_files/ \
  --outdir /path/to/results/
```

To resume a failed run or process additional files without re-computing past steps, harness Nextflow's caching engine:

```bash
nextflow run main.nf -resume
```

### 💻 System Compatibility Notice (Windows Users)

This pipeline relies on native Unix process execution architectures and optimized Linux binaries (such as GNU `sort`). 

* **macOS / Linux:** Fully supported out of the box.
* **Windows 10/11:** **Not supported natively** via CMD or PowerShell. Windows users must run this pipeline inside **Windows Subsystem for Linux (WSL)** (Ubuntu 22.04 LTS or later recommended) with Conda/Mamba configured within the WSL environment.

## 📈 Quality & Operational Reports

The pipeline isolates corrupted inputs using automated routing boundaries. If a file displays over 50% missing data in critical parameters or completely lacks the minimum attributes required for standardisation, the file is routed into a graceful soft-fail state.

Centralized audit reports are dynamically written to the output workspace:

`results/reports/pipeline_qc_failures.tsv`: Lists specific input files flagged with structural errors alongside descriptive reason strings.

`results/reports/liftover_drop_metrics.tsv`: Log summary reporting the exact count of coordinates that failed mapping parameters across individual datasets.

Finished outputs are compiled uniformly in `results/Final_GWAS/`.

## 🗺️ Roadmap (Future Version Enhancements)
[ ] Add extra options for specific liftover transformations

[ ] Add tools to address further issues that GWAS summary statistics files may come with.

[ ] Integrate native Docker/Singularity container configurations alongside Conda profiles for enterprise cloud scaling.
