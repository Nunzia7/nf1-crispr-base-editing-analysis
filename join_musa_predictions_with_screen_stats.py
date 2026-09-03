#!/usr/bin/env python3
"""Join variant prediction output with MAGeCK/BEstimate/VEP LFC and FDR.

This performs a one-to-many join:
- MuSA prediction table key: vcf_id
- merged file key: HGVS

If a HGVS location appears multiple times in the merged file, the matching
row from the MuSA prediction table is duplicated once per merged match.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

REQUESTED_MERGED_COLUMNS = [
    "sgrna",
    "Gene",
    "control_count",
    "treatment_count",
    "control_mean",
    "treat_mean",
    "LFC",
    "control_var",
    "adj_var",
    "score",
    "p.low",
    "p.high",
    "p.twosided",
    "FDR",
    "high_in_treatment",
    "CRISPR_PAM_Sequence",
    "curated_Domain",
    "Protein_Change",
]

OPTIONAL_MERGED_COLUMN_ALIASES = {
    "HGVSp": ["HGVSp"],
    "Amino_acids": ["Amino_acids"],
}

OUTPUT_MERGED_COLUMNS = [
    *REQUESTED_MERGED_COLUMNS,
    *OPTIONAL_MERGED_COLUMN_ALIASES.keys(),
]

# Public default: use a relative analysis/ directory. Override with --analysis-root.
DEFAULT_ANALYSIS_ROOT = Path("analysis")
DEFAULT_COMPARISON = "starv_vs_UN"

def output_colname(source_col: str, output_prefix: str) -> str:
    return f"{output_prefix}{source_col}"


def default_merged_path(analysis_root: Path, comparison: str) -> Path:
    return (
        analysis_root
        / "mageck_test_x_bestimate_nf1"
        / comparison
        / "merged_mageck_test_bestimate_plus_vep.csv"
    )


def default_output_path(comparison: str) -> Path:
    return Path(f"musa_predictions_with_{comparison}_screen_stats_expanded.csv")


def build_lookup(merged_path: Path, output_prefix: str):
    lookup = defaultdict(list)
    merged_rows = 0
    merged_non_na = 0

    with merged_path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"HGVS", *REQUESTED_MERGED_COLUMNS}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in merged file: {sorted(missing)}")
        available_columns = set(reader.fieldnames or [])

        optional_column_source = {
            output_col: next(
                (alias for alias in aliases if alias in available_columns),
                None,
            )
            for output_col, aliases in OPTIONAL_MERGED_COLUMN_ALIASES.items()
        }

        for row in reader:
            merged_rows += 1
            key = (row.get("HGVS") or "").strip()
            if not key or key == "NA":
                continue
            merged_non_na += 1
            match = {
                output_colname(col, output_prefix): (row.get(col) or "").strip()
                for col in REQUESTED_MERGED_COLUMNS
            }
            for output_col, source_col in optional_column_source.items():
                match[output_colname(output_col, output_prefix)] = (
                    (row.get(source_col) or "").strip() if source_col else ""
                )
            lookup[key].append(match)

    counts = Counter({k: len(v) for k, v in lookup.items()})
    stats = {
        "merged_rows": merged_rows,
        "merged_non_na_hgvs_rows": merged_non_na,
        "unique_non_na_hgvs_keys": len(lookup),
        "hgvs_keys_with_duplicates": sum(1 for v in counts.values() if v > 1),
        "max_rows_per_hgvs": max(counts.values()) if counts else 0,
    }
    return lookup, stats


def join_expanded(musa_path: Path, output_path: Path, lookup, output_prefix: str):
    musa_rows = 0
    output_rows = 0
    matched_musa_rows = 0
    unmatched_musa_rows = 0
    extra_rows_due_to_duplicates = 0

    with musa_path.open(newline="") as fin, output_path.open("w", newline="") as fout:
        reader = csv.DictReader(fin)
        if "vcf_id" not in (reader.fieldnames or []):
            raise ValueError("Missing column vcf_id in MuSA predictions file")

        added_fieldnames = [output_colname(col, output_prefix) for col in OUTPUT_MERGED_COLUMNS]
        fieldnames = list(reader.fieldnames) + added_fieldnames
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            musa_rows += 1
            key = (row.get("vcf_id") or "").strip()
            matches = lookup.get(key)

            if matches:
                matched_musa_rows += 1
                extra_rows_due_to_duplicates += max(0, len(matches) - 1)
                for match in matches:
                    out_row = dict(row)
                    out_row.update(match)
                    writer.writerow(out_row)
                    output_rows += 1
            else:
                unmatched_musa_rows += 1
                out_row = dict(row)
                for col in OUTPUT_MERGED_COLUMNS:
                    out_row[output_colname(col, output_prefix)] = ""
                writer.writerow(out_row)
                output_rows += 1

    return {
        "original_musa_rows": musa_rows,
        "matched_musa_rows": matched_musa_rows,
        "unmatched_musa_rows": unmatched_musa_rows,
        "extra_rows_due_to_duplicated_hgvs": extra_rows_due_to_duplicates,
        "expanded_output_rows": output_rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--musa", default="predictions_musa.csv", type=Path, help="MuSA-derived prediction/annotation CSV")
    parser.add_argument(
        "--comparison",
        default=DEFAULT_COMPARISON,
        help=(
            "Comparison name used to derive the default merged file, output file, "
            "and output-column prefix. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--analysis-root",
        default=DEFAULT_ANALYSIS_ROOT,
        type=Path,
        help=(
            "Analysis directory containing mageck_test_x_bestimate_nf1/. "
            "Used only when --merged is not provided. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--merged",
        default=None,
        type=Path,
        help=(
            "Merged MAGeCK/BEstimate/VEP CSV. If omitted, derived from "
            "--analysis-root and --comparison."
        ),
    )
    parser.add_argument(
        "--out",
        default=None,
        type=Path,
        help="Output CSV. If omitted, derived from --comparison.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Prefix for added output columns. If omitted, uses "
            "'<comparison>_'."
        ),
    )
    args = parser.parse_args()

    merged_path = args.merged or default_merged_path(args.analysis_root, args.comparison)
    output_path = args.out or default_output_path(args.comparison)
    output_prefix = args.output_prefix or f"{args.comparison}_"

    lookup, merged_stats = build_lookup(merged_path, output_prefix)
    join_stats = join_expanded(args.musa, output_path, lookup, output_prefix)

    print(f"Merged input: {merged_path.resolve()}")
    print(f"Output prefix: {output_prefix}")
    print(f"Created: {output_path.resolve()}")
    for label, stats in (("Merged stats", merged_stats), ("Join stats", join_stats)):
        print(f"\n{label}:")
        for key, value in stats.items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    main()



# Example:
# python join_musa_predictions_with_screen_stats.py --comparison starv_vs_UN
