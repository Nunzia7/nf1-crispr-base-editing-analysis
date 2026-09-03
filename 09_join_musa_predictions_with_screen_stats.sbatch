#!/bin/bash
#SBATCH --job-name=join_musa_screen_stats
#SBATCH --output=logs/join_musa_screen_stats_%x_%j.out
#SBATCH --error=logs/join_musa_screen_stats_%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

###############################################################################
# Join variant predictions with MAGeCK/BEstimate/VEP LFC/FDR on an HPC/SLURM cluster
#
# Default run:
#   sbatch 09_join_musa_predictions_with_screen_stats.sbatch
#
# Override comparison if needed:
#   sbatch --export=ALL,COMPARISON=starv_vs_UN 09_join_musa_predictions_with_screen_stats.sbatch
#
# Override python/env if needed:
#   sbatch --export=ALL,PYTHON=/path/to/python 09_join_musa_predictions_with_screen_stats.sbatch
###############################################################################

echo "[$(date)] Job starting on ${SLURMD_NODENAME:-unknown_node}"

# Public/reusable defaults. Override these paths with --export or environment variables.
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
COMPARISON=${COMPARISON:-starv_vs_UN}

ANALYSIS_ROOT=${ANALYSIS_ROOT:-${PROJECT_DIR}/analysis}
MUSA_FILE=${MUSA_FILE:-${PROJECT_DIR}/after_musa/predictions_musa.csv}
OUT_FILE=${OUT_FILE:-${ANALYSIS_ROOT}/after_musa/musa_predictions_with_${COMPARISON}_screen_stats_expanded.csv}
SCRIPT=${SCRIPT:-${PROJECT_DIR}/join_musa_predictions_with_screen_stats.py}

mkdir -p "logs"
mkdir -p "$(dirname "${OUT_FILE}")"

cd "${PROJECT_DIR}"

PYTHON=${PYTHON:-python}

echo "PROJECT_DIR: ${PROJECT_DIR}"
echo "COMPARISON: ${COMPARISON}"
echo "ANALYSIS_ROOT: ${ANALYSIS_ROOT}"
echo "MUSA_FILE: ${MUSA_FILE}"
echo "OUT_FILE: ${OUT_FILE}"
echo "SCRIPT: ${SCRIPT}"
echo "PYTHON: ${PYTHON}"

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${SCRIPT}"

echo "[$(date)] Running join"
"${PYTHON}" "${SCRIPT}" \
  --comparison "${COMPARISON}" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --musa "${MUSA_FILE}" \
  --out "${OUT_FILE}"

echo "[$(date)] Output file: ${OUT_FILE}"
ls -lh "${OUT_FILE}"

echo "[$(date)] Job finished"
