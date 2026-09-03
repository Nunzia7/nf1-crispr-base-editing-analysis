#!/bin/bash
#SBATCH --job-name=mageck_test_trim22
#SBATCH --output=mageck_test_trim22_%j.out
#SBATCH --error=mageck_test_trim22_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

################################################################################
# CRISPR/Cas9 HT NGS pipeline - MAGeCK test
#
# Runs a treatment-vs-control comparison using the count table produced by
# 02_run_mageck_count.sh. By default this example compares starved vs untreated cells.
################################################################################

set -euo pipefail

MAGECK_ENV="${MAGECK_ENV:-mageck}"
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT_DIR}/analysis}"
OUT_DIR="${OUT_DIR:-${ANALYSIS_DIR}/mageck_test_trim_22}"
mkdir -p "${OUT_DIR}"

if [[ -n "${MAGECK_ENV}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${MAGECK_ENV}"
  else
    echo "WARNING: conda not found; assuming mageck is already available on PATH" >&2
  fi
fi

COUNT_FILE="${COUNT_FILE:-${ANALYSIS_DIR}/mageck_count_trim_22_aug_2026.count.txt}"
FIXED_COUNT_FILE="${FIXED_COUNT_FILE:-${OUT_DIR}/mageck_count_trim_22_aug_2026.fixed_gene.count.txt}"
NEGATIVE_CONTROLS="${NEGATIVE_CONTROLS:-${PROJECT_DIR}/input/negative_controls.txt}"
TREATMENT_LABEL="${TREATMENT_LABEL:-S68582_HCC_starv}"
CONTROL_LABEL="${CONTROL_LABEL:-S68581_HCC_UN}"
COMPARISON_NAME="${COMPARISON_NAME:-starv_vs_UN}"

# MAGeCK expects each sgRNA to have a gene label. Assign missing gene labels to
# NonTargeting in a local copy of the count file.
awk -F '\t' 'BEGIN{OFS="\t"} NR==1 {print; next} $2=="" {$2="NonTargeting"} {print}' \
  "${COUNT_FILE}" > "${FIXED_COUNT_FILE}"

mageck test \
  -k "${FIXED_COUNT_FILE}" \
  -t "${TREATMENT_LABEL}" \
  -c "${CONTROL_LABEL}" \
  --norm-method control \
  --control-sgrna "${NEGATIVE_CONTROLS}" \
  -n "${OUT_DIR}/${COMPARISON_NAME}"
