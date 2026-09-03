#!/usr/bin/env python3
"""Merge BEstimate guide annotations with MAGeCK sgRNA-level results.

The script combines:
- BEstimate output tables: ``*_crispr_df.csv`` and ``*_edit_df.csv``
- an annotated guide library, usually ``input/NGN_library.csv``
- MAGeCK test ``*.sgrna_summary.txt`` files

Outputs are written per comparison and include:
- ``merged_mageck_test_bestimate.tsv``: full merged table
- ``top_enriched.csv`` / ``top_depleted.csv``: ranked guide tables
- ``clean_single_edit_cds_FDRltX.csv``: significant guides with cleaner editing profile

All paths are configurable for public/reusable use. Example:

    python 07_merge_mageck_bestimate.py \
      --project-dir /path/to/project \
      --analysis-dir /path/to/project/analysis \
      --comparison starv_vs_UN
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional

import pandas as pd


def read_bestimate_tables(
    crispr_df: Path,
    edit_df: Path,
    annotated: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read BEstimate and guide annotation tables."""
    crispr = pd.read_csv(crispr_df)
    edit = pd.read_csv(edit_df)
    ann = pd.read_csv(annotated)
    return crispr, edit, ann


def read_mageck_sgrna_summary(path: Path) -> pd.DataFrame:
    """Read a tab-delimited MAGeCK sgRNA summary file."""
    return pd.read_csv(path, sep="\t")


def build_bestimate_per_guide(
    crispr: pd.DataFrame,
    edit: pd.DataFrame,
    ann: pd.DataFrame,
) -> pd.DataFrame:
    """Create one BEstimate annotation row per CRISPR_PAM_Sequence."""
    required_cols = [
        "CRISPR_PAM_Sequence",
        "gRNA_Target_Sequence",
        "Location",
        "Direction",
        "guide_in_CDS",
        "Poly_T",
        "GC%",
        "Transcript_ID",
        "Exon_ID",
    ]
    missing = [col for col in required_cols if col not in crispr.columns]
    if missing:
        raise ValueError(f"crispr_df is missing required columns: {missing}")

    base_cols = required_cols + (["Hugo_Symbol"] if "Hugo_Symbol" in crispr.columns else [])
    base = crispr[base_cols].drop_duplicates("CRISPR_PAM_Sequence").copy()

    if not {"CRISPR_PAM_Sequence", "Edit_Location"}.issubset(edit.columns):
        raise ValueError("edit_df must contain CRISPR_PAM_Sequence and Edit_Location")

    agg = {"Edit_Location": "count"}
    if "Edit_in_Exon" in edit.columns:
        agg["Edit_in_Exon"] = "sum"
    if "Edit_in_CDS" in edit.columns:
        agg["Edit_in_CDS"] = "sum"

    per = (
        edit.groupby("CRISPR_PAM_Sequence", dropna=False)
        .agg(agg)
        .rename(
            columns={
                "Edit_Location": "n_edit_positions",
                "Edit_in_Exon": "n_edit_in_exon",
                "Edit_in_CDS": "n_edit_in_cds",
            }
        )
        .reset_index()
    )

    out = base.merge(per, on="CRISPR_PAM_Sequence", how="left")
    out["n_edit_positions"] = out["n_edit_positions"].fillna(0).astype(int)

    if "n_edit_in_exon" not in out.columns:
        out["n_edit_in_exon"] = 0
    if "n_edit_in_cds" not in out.columns:
        out["n_edit_in_cds"] = 0

    out["has_edit_in_exon"] = out["n_edit_in_exon"].fillna(0).astype(int) > 0
    out["has_edit_in_cds"] = out["n_edit_in_cds"].fillna(0).astype(int) > 0

    # Add guide ID/gene information from the annotated library when available.
    if "CRISPR_PAM_Sequence" in ann.columns:
        ann_key = "CRISPR_PAM_Sequence"
    elif "gRNA_Target_Sequence" in ann.columns:
        ann_key = "gRNA_Target_Sequence"
    else:
        ann_key = None

    if ann_key:
        ann_cols = [ann_key] + [col for col in ["ID", "Gene", "Hugo_Symbol"] if col in ann.columns]
        ann_small = ann[ann_cols].drop_duplicates(ann_key)
        out = out.merge(ann_small, left_on=ann_key, right_on=ann_key, how="left", suffixes=("", "_library"))

    return out


def merge_mageck_with_bestimate(
    mageck: pd.DataFrame,
    best_per_guide: pd.DataFrame,
    join_left: str = "sgrna",
    join_right: str = "ID",
) -> pd.DataFrame:
    """Merge MAGeCK sgRNA rows with one-row-per-guide BEstimate annotations."""
    if join_left not in mageck.columns:
        raise ValueError(f"MAGeCK table is missing join column: {join_left}")
    if join_right not in best_per_guide.columns:
        raise ValueError(f"BEestimate table is missing join column: {join_right}")

    return mageck.merge(
        best_per_guide,
        left_on=join_left,
        right_on=join_right,
        how="left",
        indicator=True,
    )


def add_cleanliness_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Add simple filtering flags for candidate guide prioritization."""
    out = df.copy()
    if "n_edit_positions" in out.columns:
        out["single_edit"] = out["n_edit_positions"].fillna(0).astype(int) == 1
    if "has_edit_in_cds" in out.columns:
        out["edit_in_cds"] = out["has_edit_in_cds"].fillna(False).astype(bool)
    if "Poly_T" in out.columns:
        out["no_poly_t"] = ~out["Poly_T"].fillna(False).astype(bool)
    if "GC%" in out.columns:
        gc = pd.to_numeric(out["GC%"], errors="coerce")
        out["GC_ok"] = gc.between(30, 75)
    flags = [col for col in ["single_edit", "edit_in_cds", "no_poly_t", "GC_ok"] if col in out.columns]
    if flags:
        out["clean_single_edit_cds"] = out[flags].all(axis=1)
    return out


def write_outputs(
    merged: pd.DataFrame,
    outdir: Path,
    focus_gene: Optional[str],
    fdr_threshold: float,
    top_n: int,
) -> dict[str, Path]:
    """Write merged and ranked output tables."""
    outdir.mkdir(parents=True, exist_ok=True)

    merged_path = outdir / "merged_mageck_test_bestimate.tsv"
    merged.to_csv(merged_path, sep="\t", index=False)

    df = merged.copy()
    if focus_gene and "Gene" in df.columns:
        df = df[df["Gene"].astype(str) == focus_gene].copy()

    merged_only = df[df["_merge"] == "both"].copy() if "_merge" in df.columns else df

    enriched = merged_only.sort_values(["FDR", "LFC"], ascending=[True, False]).head(top_n)
    depleted = merged_only.sort_values(["FDR", "LFC"], ascending=[True, True]).head(top_n)

    enriched_path = outdir / "top_enriched.csv"
    depleted_path = outdir / "top_depleted.csv"
    enriched.to_csv(enriched_path, index=False)
    depleted.to_csv(depleted_path, index=False)

    clean_sig = merged_only.copy()
    if "FDR" in clean_sig.columns:
        clean_sig = clean_sig[pd.to_numeric(clean_sig["FDR"], errors="coerce") < fdr_threshold]
    if "clean_single_edit_cds" in clean_sig.columns:
        clean_sig = clean_sig[clean_sig["clean_single_edit_cds"] == True]
    if "GC_ok" in clean_sig.columns:
        clean_sig = clean_sig[clean_sig["GC_ok"].fillna(False)]
    clean_sig = clean_sig.sort_values(["FDR", "LFC"], ascending=[True, False]).head(top_n)

    clean_sig_path = outdir / f"clean_single_edit_cds_FDRlt{fdr_threshold:g}.csv"
    clean_sig.to_csv(clean_sig_path, index=False)

    return {
        "merged": merged_path,
        "top_enriched": enriched_path,
        "top_depleted": depleted_path,
        "clean_sig": clean_sig_path,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    default_project_dir = Path(os.environ.get("PROJECT_DIR", "."))
    default_analysis_dir = Path(os.environ.get("ANALYSIS_DIR", default_project_dir / "analysis"))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, default=default_project_dir)
    parser.add_argument("--analysis-dir", type=Path, default=default_analysis_dir)
    parser.add_argument("--bestimate-dir", type=Path, default=None)
    parser.add_argument("--comparison", action="append", default=None, help="Comparison name. Can be provided multiple times.")
    parser.add_argument("--crispr-df", type=Path, default=None)
    parser.add_argument("--edit-df", type=Path, default=None)
    parser.add_argument("--annotated-library", type=Path, default=None)
    parser.add_argument("--mageck-dir", type=Path, default=None)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--focus-gene", default="NF1", help="Set to empty string to keep all genes.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--join-left", default="sgrna")
    parser.add_argument("--join-right", default="ID")
    return parser.parse_args()


def main() -> None:
    """Run the merge for one or more comparisons."""
    args = parse_args()

    project_dir = args.project_dir
    analysis_dir = args.analysis_dir
    bestimate_dir = args.bestimate_dir or analysis_dir / "bestimate_out_plus_vep_full"
    comparisons = args.comparison or ["starv_vs_UN"]

    crispr_df = args.crispr_df or bestimate_dir / "NF1_NGN_A2G_act4-8_crispr_df.csv"
    edit_df = args.edit_df or bestimate_dir / "NF1_NGN_A2G_act4-8_edit_df.csv"
    annotated = args.annotated_library or project_dir / "input" / "NGN_library.csv"
    mageck_dir = args.mageck_dir or analysis_dir / "mageck_test_trim_22"
    out_root = args.out_root or analysis_dir / "mageck_test_x_bestimate_nf1"
    focus_gene = args.focus_gene or None

    crispr, edit, ann = read_bestimate_tables(crispr_df, edit_df, annotated)
    best_per_guide = build_bestimate_per_guide(crispr=crispr, edit=edit, ann=ann)

    for comparison in comparisons:
        mageck_sgrna_summary = mageck_dir / f"{comparison}.sgrna_summary.txt"
        outdir = out_root / comparison

        mageck = read_mageck_sgrna_summary(mageck_sgrna_summary)
        merged = merge_mageck_with_bestimate(
            mageck=mageck,
            best_per_guide=best_per_guide,
            join_left=args.join_left,
            join_right=args.join_right,
        )
        merged = add_cleanliness_flags(merged)

        print(f"\nComparison: {comparison}")
        print("Merge counts:")
        if "_merge" in merged.columns:
            print(merged["_merge"].value_counts(dropna=False).to_string())

        paths = write_outputs(
            merged=merged,
            outdir=outdir,
            focus_gene=focus_gene,
            fdr_threshold=args.fdr_threshold,
            top_n=args.top_n,
        )
        for label, path in paths.items():
            print(f"{label}: {path}")


if __name__ == "__main__":
    main()
