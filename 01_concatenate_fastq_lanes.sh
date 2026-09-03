#!/bin/bash
# Concatenate lane-level R1 FASTQ files into one FASTQ per sample.
#
# Expected input layout:
#   ${INPUT_DIR}/Sample_*/<sample>_L00*_R1_001.fastq.gz
#
# Usage:
#   INPUT_DIR=/path/to/raw_fastq OUTPUT_DIR=/path/to/analysis bash 01_concatenate_fastq_lanes.sh

set -euo pipefail

INPUT_DIR="${INPUT_DIR:-data/raw_fastq}"
OUTPUT_DIR="${OUTPUT_DIR:-data/processed_fastq}"

mkdir -p "${OUTPUT_DIR}"

for sample_dir in "${INPUT_DIR}"/Sample_*; do
    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")
    echo "Processing ${sample_name}"

    # Sort lane files to keep a reproducible concatenation order.
    mapfile -t r1_files < <(find "${sample_dir}" -maxdepth 1 -name '*_L00*_R1_001.fastq.gz' | sort)

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "WARNING: no R1 FASTQ files found in ${sample_dir}" >&2
        continue
    fi

    cat "${r1_files[@]}" > "${OUTPUT_DIR}/${sample_name}.fastq.gz"
done
