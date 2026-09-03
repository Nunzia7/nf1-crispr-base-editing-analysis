#!/usr/bin/env python3
"""Generate NF1 sgRNA-level volcano plots and significant guide tables.

Reads MAGeCK ``*.sgrna_summary.txt`` files from an input directory and writes
plots/tables to an output directory. This keeps the analysis code in a normal
Python file instead of embedding Python inside a bash script.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def classify_hit(row: pd.Series, lfc_threshold: float, fdr_threshold: float) -> str:
    """Classify one sgRNA row based on LFC/FDR thresholds."""
    if pd.isna(row["FDR"]) or pd.isna(row["LFC"]):
        return "NA"
    if row["FDR"] < fdr_threshold and row["LFC"] <= -lfc_threshold:
        return "significant_depletion"
    if row["FDR"] < fdr_threshold and row["LFC"] >= lfc_threshold:
        return "significant_enrichment"
    return "other"


def process_file(
    input_file: Path,
    out_root: Path,
    target_gene: str,
    lfc_threshold: float,
    fdr_threshold: float,
    top_n: int,
) -> dict:
    """Process one MAGeCK sgRNA summary file."""
    comparison = input_file.name.replace(".sgrna_summary.txt", "")
    outdir = out_root / comparison
    outdir.mkdir(parents=True, exist_ok=True)

    df_all = pd.read_csv(input_file, sep="\t")
    required = ["sgrna", "Gene", "LFC", "FDR"]
    missing = [col for col in required if col not in df_all.columns]
    if missing:
        raise ValueError(f"{input_file}: missing columns {missing}; available={list(df_all.columns)}")

    n_all = len(df_all)
    df = df_all[df_all["Gene"].astype(str) == target_gene].copy()
    if df.empty:
        raise ValueError(f"{input_file}: no sgRNAs found with Gene == {target_gene!r}")

    df["LFC"] = pd.to_numeric(df["LFC"], errors="coerce")
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df = df.dropna(subset=["LFC", "FDR"]).copy()
    df["minus_log10_FDR"] = -np.log10(df["FDR"].clip(lower=1e-300))
    df["hit_class"] = df.apply(
        classify_hit,
        axis=1,
        lfc_threshold=lfc_threshold,
        fdr_threshold=fdr_threshold,
    )
    df["absLFC"] = df["LFC"].abs()

    df.to_csv(outdir / f"{comparison}.annotated_sgrna_summary.csv", index=False)

    plt.figure(figsize=(8, 6))
    plt.scatter(df["LFC"], df["FDR"], s=14, alpha=0.6)
    plt.axvline(-lfc_threshold, linestyle="--")
    plt.axvline(lfc_threshold, linestyle="--")
    plt.axhline(fdr_threshold, linestyle="--")
    plt.xlabel("LFC")
    plt.ylabel("FDR")
    plt.title(f"{comparison}: {target_gene} sgRNA LFC vs FDR")
    plt.tight_layout()
    plt.savefig(outdir / "lfc_vs_fdr_raw.png", dpi=300)
    plt.close()

    colors = {
        "significant_depletion": "tab:red",
        "significant_enrichment": "tab:green",
        "other": "tab:gray",
    }
    plt.figure(figsize=(8, 6))
    for hit_class, sub in df.groupby("hit_class"):
        plt.scatter(
            sub["LFC"],
            sub["minus_log10_FDR"],
            s=16,
            alpha=0.7,
            color=colors.get(hit_class, "tab:gray"),
            label=f"{hit_class} (n={len(sub)})",
        )
    plt.axvline(-lfc_threshold, linestyle="--")
    plt.axvline(lfc_threshold, linestyle="--")
    plt.axhline(-np.log10(fdr_threshold), linestyle="--")
    plt.xlabel("LFC")
    plt.ylabel("-log10(FDR)")
    plt.title(f"{comparison}: {target_gene} sgRNA volcano-like plot")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / "lfc_vs_minuslog10fdr.png", dpi=300)
    plt.savefig(outdir / "lfc_vs_minuslog10fdr.pdf", dpi=300)
    plt.close()

    df.sort_values("LFC").head(top_n).to_csv(outdir / f"top_{top_n}_depleted_sgrnas.csv", index=False)
    df.sort_values("LFC", ascending=False).head(top_n).to_csv(outdir / f"top_{top_n}_enriched_sgrnas.csv", index=False)

    sig_dep = df[df["hit_class"] == "significant_depletion"].sort_values("LFC")
    sig_enr = df[df["hit_class"] == "significant_enrichment"].sort_values("LFC", ascending=False)
    sig_dep.to_csv(outdir / "significant_depleted_sgrnas.csv", index=False)
    sig_enr.to_csv(outdir / "significant_enriched_sgrnas.csv", index=False)

    gene_summary = (
        df.groupby("Gene")
        .agg(
            n_guides=("sgrna", "count"),
            mean_LFC=("LFC", "mean"),
            median_LFC=("LFC", "median"),
            min_FDR=("FDR", "min"),
            n_sig_dep=("hit_class", lambda x: (x == "significant_depletion").sum()),
            n_sig_enr=("hit_class", lambda x: (x == "significant_enrichment").sum()),
        )
        .reset_index()
    )
    gene_summary.to_csv(outdir / "gene_level_summary_from_sgrnas.csv", index=False)

    n_dep = len(sig_dep)
    n_enr = len(sig_enr)
    n_other = int((df["hit_class"] == "other").sum())
    with (outdir / "report.txt").open("w") as handle:
        handle.write(f"MAGeCK sgRNA quick report: {comparison}\n")
        handle.write("=" * 60 + "\n\n")
        handle.write(f"Input file: {input_file}\n")
        handle.write(f"Target gene: {target_gene}\n")
        handle.write(f"Total sgRNAs in input file: {n_all}\n")
        handle.write(f"{target_gene} sgRNAs plotted/exported: {len(df)}\n")
        handle.write(f"LFC threshold: {lfc_threshold}\n")
        handle.write(f"FDR threshold: {fdr_threshold}\n")
        handle.write(f"Significant depletion: {n_dep}\n")
        handle.write(f"Significant enrichment: {n_enr}\n")
        handle.write(f"Other/non-significant: {n_other}\n")

    print(f"Done {comparison}: outputs written to {outdir}")
    print(df["hit_class"].value_counts(dropna=False).to_string())

    return {
        "comparison": comparison,
        "target_gene": target_gene,
        "n_input_sgrnas": n_all,
        "n_sgrnas": len(df),
        "n_significant_depletion": n_dep,
        "n_significant_enrichment": n_enr,
        "n_other": n_other,
        "median_LFC": df["LFC"].median(),
        "min_FDR": df["FDR"].min(),
        "outdir": str(outdir),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate sgRNA volcano plots from MAGeCK summaries.")
    parser.add_argument("--test-dir", required=True, type=Path, help="Directory with *.sgrna_summary.txt files")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory where outputs will be written")
    parser.add_argument("--target-gene", default="NF1")
    parser.add_argument("--lfc-threshold", default=1.0, type=float)
    parser.add_argument("--fdr-threshold", default=0.05, type=float)
    parser.add_argument("--top-n", default=30, type=int)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sgrna_files = sorted(args.test_dir.glob("*.sgrna_summary.txt"))
    if not sgrna_files:
        raise FileNotFoundError(f"No *.sgrna_summary.txt files found in {args.test_dir}")

    summary_rows = [
        process_file(
            input_file,
            args.out_dir,
            args.target_gene,
            args.lfc_threshold,
            args.fdr_threshold,
            args.top_n,
        )
        for input_file in sgrna_files
    ]

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.out_dir / "sgrna_volcano_summary.tsv", sep="\t", index=False)
    print(f"Combined summary: {args.out_dir / 'sgrna_volcano_summary.tsv'}")


if __name__ == "__main__":
    main()