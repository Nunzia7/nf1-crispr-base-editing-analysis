#!/bin/bash
#SBATCH --job-name=mageck_count
#SBATCH --output=mageck_count_%j.out
#SBATCH --error=mageck_count_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

################################################################################
# CRISPR/Cas9 HT NGS pipeline - MAGeCK count
#
# This script counts sgRNA abundance from concatenated FASTQ files.
# Paths are configurable so the repository can be public/reusable.
#
# Usage example:
#   PROJECT_DIR=/path/to/project ANALYSIS_DIR=/path/to/analysis bash 02_run_mageck_count.sh
#
# Optional environment variables:
#   MAGECK_ENV      conda environment name or prefix for MAGeCK; default: mageck
#   PROJECT_DIR     project root containing input/; default: current directory
#   ANALYSIS_DIR    folder containing FASTQs and receiving output; default: ./analysis
#   SGRNA_LIST      library annotation CSV
#   NEGATIVE_CONTROL_SGRNA_LIST  text file with negative-control sgRNA IDs
################################################################################

set -euo pipefail

MAGECK_ENV="${MAGECK_ENV:-mageck}"
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT_DIR}/analysis}"
OUT_DIR="${OUT_DIR:-${ANALYSIS_DIR}}"

mkdir -p "${OUT_DIR}"

# Activate MAGeCK environment if conda is available. If MAGeCK is already on PATH,
# this section can be skipped by setting MAGECK_ENV="".
if [[ -n "${MAGECK_ENV}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${MAGECK_ENV}"
  else
    echo "WARNING: conda not found; assuming mageck is already available on PATH" >&2
  fi
fi

FASTQ_FILES=(
  "${ANALYSIS_DIR}/Sample_S68582_HCC_starv.fastq.gz"
  "${ANALYSIS_DIR}/Sample_S68581_HCC_UN.fastq.gz"
)

SGRNA_LIST="${SGRNA_LIST:-${PROJECT_DIR}/input/library_complete.csv}"
NEGATIVE_CONTROL_SGRNA_LIST="${NEGATIVE_CONTROL_SGRNA_LIST:-${PROJECT_DIR}/input/negative_controls.txt}"
SAMPLE_LABELS="${SAMPLE_LABELS:-S68582_HCC_starv,S68581_HCC_UN}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-mageck_count_trim_22_aug_2026}"

mageck count \
  -l "${SGRNA_LIST}" \
  --fastq "${FASTQ_FILES[@]}" \
  --control-sgrna "${NEGATIVE_CONTROL_SGRNA_LIST}" \
  -n "${OUT_DIR}/${OUTPUT_PREFIX}" \
  --trim-5 "${TRIM_5:-22}" \
  --sgrna-len "${SGRNA_LEN:-20}" \
  --sample-label "${SAMPLE_LABELS}"
