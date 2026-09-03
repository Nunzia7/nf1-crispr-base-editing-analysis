# MAGeCK CRISPR analysis workflow

Public-friendly scripts for a CRISPR screen analysis workflow using MAGeCK, BEstimate, MuSA variant annotations, Python, and R.

This repository contains code only. Raw sequencing data, large intermediate files, cluster logs, and generated plots/results should not be committed.

## Workflow overview

1. `01_concatenate_fastq_lanes.sh`  
   Concatenate lane-level R1 FASTQ files into one FASTQ per sample.

2. `02_run_mageck_count.sh`  
   Run `mageck count` to quantify sgRNA abundance.

3. `03_run_mageck_test.sh`  
   Run `mageck test` for treatment-vs-control comparisons.

4. `04_generate_qc_report.sh`  
   Generate a Markdown QC/performance report from FASTQ, library, MAGeCK count, and MAGeCK test outputs.

5. `05_plot_mageck_sgrna_volcano.sh` and `plot_mageck_sgrna_volcano.py`  
   Generate sgRNA-level volcano plots and hit tables from MAGeCK `*.sgrna_summary.txt` files.

6. `06_run_bestimate.sbatch`  
   Example SLURM script to run BEstimate on a CRISPR guide library.

7. `07_merge_mageck_bestimate.py`  
   Command-line Python script for merging MAGeCK and BEstimate outputs.

8. `08_plot_bestimate_results.R`  
   R script for downstream BEstimate/MAGeCK plots and candidate-guide tables.

9. `09_join_musa_predictions_with_screen_stats.sbatch` and `join_musa_predictions_with_screen_stats.py`  
   Join MuSA variant-prediction/annotation output with MAGeCK/BEEstimate/VEP LFC/FDR annotations.

10. `10_filter_musa_variants.R`  
    Filter joined variant-prediction tables into enriched/depleted candidate categories.

11. `11_plot_musa_variant_scores.R`  
    Plot filtered NF1 variants/guides along protein coordinates using prediction scores.

## Requirements

Command-line tools:

- Bash
- SLURM, if running `*.sbatch` scripts on a cluster
- MAGeCK
- BEstimate, for base-editing guide annotation/downstream interpretation
- Python 3
- R

Python packages:

```bash
pip install -r requirements.txt
```

R packages used by `08_plot_bestimate_results.R`:

```r
install.packages(c("data.table", "dplyr", "stringr", "tidyr", "ggplot2", "scales"))
```

## BEstimate base-editing annotation input

This workflow also uses **BEstimate** outputs for base-editor guide annotation and interpretation:

> Dinçer, C., Fussing, B., Garnett, M.J. et al. **BEstimate: a computational tool for the design and interpretation of CRISPR base editing experiments.** *Genome Biology* 27, 191 (2026). https://doi.org/10.1186/s13059-026-04077-z

BEestimate is a computational pipeline for CRISPR base-editing experiments. It identifies base-editor gRNA target sites and provides on-target/off-target predictions plus functional, structural, and clinical annotations for installed variants. In this repository, BEstimate-derived guide annotations are merged with MAGeCK screen statistics so candidate guides can be prioritized by both screening signal and predicted editing consequences.

## MuSA variant annotation input

The `after_musa/` input files are expected to come from **MuSA**:

> Scognamiglio, D., Bonetti, E., Moroni, A. et al. **MuSA: a Nextflow pipeline for deep, reproducible annotation and clinical ranking of genomic variants.** *BMC Bioinformatics* 27, 199 (2026). https://doi.org/10.1186/s12859-026-06513-0

MuSA, Multi-Source variant Annotation, is an nf-core-compliant Nextflow workflow for germline variant annotation and clinical ranking. In this CRISPR workflow, MuSA-derived variant prediction/annotation tables are joined to MAGeCK/BEEstimate guide-level results so variants can be prioritized together with screen statistics such as LFC and FDR.

## Expected project layout

The scripts are configurable, but by default they expect a structure similar to:

```text
project/
├── input/
│   ├── library_complete.csv
│   ├── negative_controls.txt
│   └── NGN_library.csv
├── data/
│   ├── raw_fastq/
│   └── processed_fastq/
├── after_musa/
│   └── predictions_musa.csv        # MuSA-derived variant annotation/prediction table
└── analysis/
```

## Running the scripts

All paths can be overridden with environment variables. Example:

```bash
PROJECT_DIR=/path/to/project \
ANALYSIS_DIR=/path/to/project/analysis \
INPUT_DIR=/path/to/project/data/raw_fastq \
OUTPUT_DIR=/path/to/project/data/processed_fastq \
bash 01_concatenate_fastq_lanes.sh
```

MAGeCK count example:

```bash
PROJECT_DIR=/path/to/project \
ANALYSIS_DIR=/path/to/project/analysis \
MAGECK_ENV=mageck \
bash 02_run_mageck_count.sh
```

MAGeCK test example:

```bash
PROJECT_DIR=/path/to/project \
ANALYSIS_DIR=/path/to/project/analysis \
bash 03_run_mageck_test.sh
```

Merge MAGeCK and BEstimate example:

```bash
PROJECT_DIR=/path/to/project \
ANALYSIS_DIR=/path/to/project/analysis \
python 07_merge_mageck_bestimate.py --comparison starv_vs_UN
```

Join MuSA variant annotations/predictions example:

```bash
PROJECT_DIR=/path/to/project \
ANALYSIS_ROOT=/path/to/project/analysis \
python join_musa_predictions_with_screen_stats.py \
  --comparison starv_vs_UN \
  --musa /path/to/project/after_musa/predictions_musa.csv \
  --out /path/to/project/analysis/after_musa/musa_predictions_with_starv_vs_UN_screen_stats_expanded.csv
```

Volcano plot example:

```bash
TEST_DIR=/path/to/project/analysis/mageck_test_trim_22 \
OUT_DIR=/path/to/project/analysis/mageck_test_trim_22/volcano_plots_NF1_only \
TARGET_GENE=NF1 \
bash 05_plot_mageck_sgrna_volcano.sh
```

## Notes for public release

- Do not commit FASTQ/BAM/VCF files or other large sequencing outputs.
- Do not commit SLURM `.out`/`.err` logs.
- Keep local HPC paths outside the code by passing them through environment variables.
- Start with a private repository if you are unsure whether sample metadata can be public, then switch to public after review.
