#!/bin/bash
#SBATCH --job-name=mageck_sgrna_volcano
#SBATCH --output=mageck_sgrna_volcano_%j.out
#SBATCH --error=mageck_sgrna_volcano_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Wrapper around plot_mageck_sgrna_volcano.py.
# Usage:
#   TEST_DIR=/path/to/mageck_test_trim_22 OUT_DIR=/path/to/plots bash 05_plot_mageck_sgrna_volcano.sh

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-${SCRIPT_DIR}/plot_mageck_sgrna_volcano.py}"

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
  echo "ERROR: cannot find plot_mageck_sgrna_volcano.py at: ${PYTHON_SCRIPT}" >&2
  exit 1
fi

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT_DIR}/analysis}"
TEST_DIR="${TEST_DIR:-${ANALYSIS_DIR}/mageck_test_trim_22}"
OUT_DIR="${OUT_DIR:-${TEST_DIR}/volcano_plots_NF1_only}"

LFC_THRESHOLD="${LFC_THRESHOLD:-1}"
FDR_THRESHOLD="${FDR_THRESHOLD:-0.05}"
TOP_N="${TOP_N:-30}"
TARGET_GENE="${TARGET_GENE:-NF1}"

echo "[$(date)] Starting sgRNA volcano plot generation"
echo "TEST_DIR=${TEST_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "TARGET_GENE=${TARGET_GENE}"

python "${PYTHON_SCRIPT}" \
  --test-dir "${TEST_DIR}" \
  --out-dir "${OUT_DIR}" \
  --target-gene "${TARGET_GENE}" \
  --lfc-threshold "${LFC_THRESHOLD}" \
  --fdr-threshold "${FDR_THRESHOLD}" \
  --top-n "${TOP_N}"

echo "[$(date)] sgRNA volcano plot generation completed"
