#!/usr/bin/env python3
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _load_snps(path: str) -> set:
    if not os.path.exists(path):
        return set()
    df = pd.read_csv(path, sep="\t")
    return set(df["SNP"].astype(str)) if "SNP" in df.columns else set()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Manhattan style view with method highlights")
    parser.add_argument("--matched-tsv", required=True)
    parser.add_argument("--finemap-tsv", required=True)
    parser.add_argument("--susier-tsv", required=True)
    parser.add_argument("--cojo-tsv", required=True)
    parser.add_argument("--plot-png", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.plot_png)
    _mkdir_for(args.done_file)

    df = pd.read_csv(args.matched_tsv, sep="\t")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce").clip(lower=1e-300)
    df = df.dropna(subset=["BP", "P"]).copy()
    df["minus_log10_p"] = -np.log10(df["P"])
    df["SNP"] = df["SNP"].astype(str)

    finemap = _load_snps(args.finemap_tsv)
    susier = _load_snps(args.susier_tsv)
    cojo = _load_snps(args.cojo_tsv)

    df["category"] = "background"
    df.loc[df["SNP"].isin(finemap), "category"] = "finemap"
    df.loc[df["SNP"].isin(susier), "category"] = "susier"
    df.loc[df["SNP"].isin(cojo), "category"] = "cojo"

    colors = {
        "background": "#9CA3AF",
        "finemap": "#D97706",
        "susier": "#2563EB",
        "cojo": "#059669",
    }

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for category, sub in df.groupby("category", sort=False):
        ax.scatter(sub["BP"], sub["minus_log10_p"], s=12, alpha=0.8, c=colors.get(category, "#9CA3AF"), label=category)

    ax.axhline(-np.log10(5e-8), color="#DC2626", linestyle="--", linewidth=1)
    ax.set_xlabel("Position (BP)")
    ax.set_ylabel("-log10(P)")
    ax.set_title("Locus Manhattan Plot with FINEMAP/SuSiE/COJO highlights")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.plot_png, dpi=150)
    plt.close(fig)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
