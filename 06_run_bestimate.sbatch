#!/bin/bash
#SBATCH --job-name=bestimate
#SBATCH --output=bestimate_%j.out
#SBATCH --error=bestimate_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

###############################################################################
# BEstimate SLURM script
#
# Usage examples:
#   sbatch run_bestimate.sbatch
#
#   # Override parameters at submission time:
#   sbatch --export=ALL,EDIT=C,EDIT_TO=T,OFILE=NF1_NGN_CBE_act4-8 run_bestimate.sbatch
#
# Notes:
# - This script assumes you already have a conda env named "bestimate".
# - By default it runs WITHOUT -vep (faster; no external API calls).
#   To enable VEP (needs internet access from compute nodes), submit with:
#     sbatch --export=ALL,USE_VEP=1 run_bestimate.sbatch
###############################################################################

# -----------------------------
# Job context
# -----------------------------

echo "[$(date)] Job starting on ${SLURMD_NODENAME:-unknown_node}"

# Capture the submission directory so we can refer to input/output with absolute paths
SUBMIT_DIR=${SLURM_SUBMIT_DIR:-"$(pwd)"}
echo "Submit directory: ${SUBMIT_DIR}"

# Project/run paths. Override these variables at submission time for your cluster/project.
PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
ANALYSIS_DIR=${ANALYSIS_DIR:-${PROJECT_DIR}/analysis}
BESTIMATE_REPO=${BESTIMATE_REPO:-${PROJECT_DIR}/external/BEstimate}

# -----------------------------
# Parameters (override via --export)
# -----------------------------
GENE=${GENE:-NF1}
ASSEMBLY=${ASSEMBLY:-GRCh38}

LIBRARY_FILE=${LIBRARY_FILE:-${PROJECT_DIR}/input/NGN_library.csv}

PAMSEQ=${PAMSEQ:-NGN}
PAMWIN=${PAMWIN:-21-23}
ACTWIN=${ACTWIN:-4-8}

# protospacer+PAM length. Your CRISPR_PAM_Sequence is 23 nt (20+3)
PROTOLEN=${PROTOLEN:-23}

# Pick one of:
#   ABE: EDIT=A EDIT_TO=G
#   CBE: EDIT=C EDIT_TO=T
EDIT=${EDIT:-A}
EDIT_TO=${EDIT_TO:-G}

OUTDIR=${OUTDIR:-${ANALYSIS_DIR}/bestimate_out_plus_vep_full}
OFILE=${OFILE:-${GENE}_${PAMSEQ}_${EDIT}2${EDIT_TO}_act${ACTWIN}}

USE_VEP=${USE_VEP:-1}

# Workaround switch (SBATCH-only; no edits to BEstimate source):
#
# Some setups hit:
#   ValueError: Columns must be same length as key
# inside BEstimate's `extract_grna_sites()`.
#
# BEstimate will SKIP `extract_grna_sites()` if the file
#   ${OUTDIR}/${OFILE}_crispr_df.csv
# already exists.
#
# If USE_LIBRARY_AS_CRISPR_DF=1, we pre-create that file from your
# `-library_file` input so BEstimate reads it instead of recomputing.
USE_LIBRARY_AS_CRISPR_DF=${USE_LIBRARY_AS_CRISPR_DF:-1}

# Make paths absolute (important because we will `cd` into the BEstimate repo folder)
LIBRARY_FILE_ABS="${LIBRARY_FILE}"
if [[ "${LIBRARY_FILE_ABS}" != /* ]]; then
  LIBRARY_FILE_ABS="${SUBMIT_DIR}/${LIBRARY_FILE_ABS}"
fi

OUTDIR_ABS="${OUTDIR}"
if [[ "${OUTDIR_ABS}" != /* ]]; then
  OUTDIR_ABS="${SUBMIT_DIR}/${OUTDIR_ABS}"
fi

mkdir -p "${OUTDIR_ABS}"

# Override if your env lives somewhere else:
ENV_PREFIX=${ENV_PREFIX:-"${HOME}/.conda/envs/bestimate"}

if [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "ERROR: ENV_PREFIX does not exist: ${ENV_PREFIX}"
  echo "Submit with: sbatch --export=ALL,ENV_PREFIX=/full/path/to/env run_bestimate.sbatch"
  exit 1
fi

export PATH="${ENV_PREFIX}/bin:${PATH}"
PYTHON="${ENV_PREFIX}/bin/python"

if [[ ! -x "${PYTHON}" ]]; then
  echo "ERROR: python not found/executable at: ${PYTHON}"
  exit 1
fi

echo "Using ENV_PREFIX: ${ENV_PREFIX}"
echo "Using python: ${PYTHON}"
"${PYTHON}" -V

echo "Pandas version: $(${PYTHON} -c 'import pandas as pd; print(pd.__version__)')"

# -----------------------------
# Run BEstimate
# -----------------------------

# Optional workaround: pre-create ${OFILE}_crispr_df.csv from your library file
if [[ "${USE_LIBRARY_AS_CRISPR_DF}" == "1" ]]; then
  CRISPR_DF_OUT="${OUTDIR_ABS}/${OFILE}_crispr_df.csv"
  echo "USE_LIBRARY_AS_CRISPR_DF=1 -> writing ${CRISPR_DF_OUT} from ${LIBRARY_FILE_ABS}"

  "${PYTHON}" - <<PY
import pandas as pd

lib = pd.read_csv(r"${LIBRARY_FILE_ABS}")

required = [
    "Hugo_Symbol",
    "CRISPR_PAM_Sequence",
    "gRNA_Target_Sequence",
    "Location",
    "Direction",
    "Gene_ID",
    "Transcript_ID",
    "Exon_ID",
    "guide_in_CDS",
]

missing = [c for c in required if c not in lib.columns]
if missing:
    raise SystemExit(f"Missing required columns in library file: {missing}")

crispr_df = lib[required].copy()

# Columns BEstimate normally appends in extract_grna_sites()
crispr_df["gRNA_flanking_sequences"] = None
crispr_df["Poly_T"] = crispr_df["gRNA_Target_Sequence"].astype(str).str.contains(r"T{4,}", regex=True)
seq = crispr_df["gRNA_Target_Sequence"].astype(str)
crispr_df["GC%"] = ((seq.str.count("G") + seq.str.count("C")) * 100.0 / seq.str.len()).fillna(0)

# Make sure guide_in_CDS is boolean-ish
crispr_df["guide_in_CDS"] = crispr_df["guide_in_CDS"].astype(bool)

crispr_df.to_csv(r"${CRISPR_DF_OUT}", index=False)
print("Wrote", r"${CRISPR_DF_OUT}", "rows=", len(crispr_df))
PY
fi

# IMPORTANT:
# BEstimate-main/BEstimate/BEstimate.py expects to find data files via
# paths like: os.getcwd() + "/../data/H_sapiens_interfaces.txt"
# Therefore we run it FROM the BEstimate-main/BEstimate directory.
cd "${BESTIMATE_REPO}"
echo "Running from: $(pwd)"

VEP_FLAG=()
if [[ "${USE_VEP}" == "1" ]]; then
  VEP_FLAG=(-vep)
fi

# Prefer installed CLI in the env if available, otherwise fall back to repo script
if [[ -x "${ENV_PREFIX}/bin/BEstimate" ]]; then
  CMD=(
    "${ENV_PREFIX}/bin/BEstimate"
    -gene "${GENE}" -assembly "${ASSEMBLY}"
    -pamseq "${PAMSEQ}" -pamwin "${PAMWIN}"
    -actwin "${ACTWIN}" -protolen "${PROTOLEN}"
    -library_file "${LIBRARY_FILE_ABS}"
    -edit "${EDIT}" -edit_to "${EDIT_TO}"
    -o "${OUTDIR_ABS}" -ofile "${OFILE}"
    "${VEP_FLAG[@]}"
  )
else
  CMD=(
    "${PYTHON}" BEstimate.py
    -gene "${GENE}" -assembly "${ASSEMBLY}"
    -pamseq "${PAMSEQ}" -pamwin "${PAMWIN}"
    -actwin "${ACTWIN}" -protolen "${PROTOLEN}"
    -library_file "${LIBRARY_FILE_ABS}"
    -edit "${EDIT}" -edit_to "${EDIT_TO}"
    -o "${OUTDIR_ABS}" -ofile "${OFILE}"
    "${VEP_FLAG[@]}"
  )
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

echo "[$(date)] Job finished"

