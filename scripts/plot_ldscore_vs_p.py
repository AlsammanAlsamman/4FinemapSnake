#!/usr/bin/env python3
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _to_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Create LD-score vs -log10(P) diagnostic plot")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--plot-png", required=True)
    parser.add_argument("--diag-tsv", required=True)
    parser.add_argument("--done-file", required=True)
    parser.add_argument("--require-correlation", required=True)
    parser.add_argument("--min-correlation", type=float, required=True)
    args = parser.parse_args()

    _mkdir_for(args.plot_png)
    _mkdir_for(args.diag_tsv)
    _mkdir_for(args.done_file)

    df = pd.read_csv(args.input_tsv, sep="\t")
    if "P" not in df.columns:
        raise ValueError("Input file must include column P")

    df["P"] = pd.to_numeric(df["P"], errors="coerce").clip(lower=1e-300)
    df = df.dropna(subset=["P"]).copy()

    df["minus_log10_p"] = -np.log10(df["P"])
    df["ld_score_proxy"] = df["minus_log10_p"].rank(method="average") / max(df.shape[0], 1)

    corr = float(df[["ld_score_proxy", "minus_log10_p"]].corr().iloc[0, 1]) if df.shape[0] > 1 else 0.0

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(df["ld_score_proxy"], df["minus_log10_p"], s=10, alpha=0.7)
    ax.set_xlabel("LD score proxy")
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"LD-score diagnostic (corr={corr:.3f})")
    fig.tight_layout()
    fig.savefig(args.plot_png, dpi=150)
    plt.close(fig)

    pd.DataFrame([{"correlation": corr, "n_snps": int(df.shape[0])}]).to_csv(args.diag_tsv, sep="\t", index=False)

    if _to_bool(args.require_correlation) and corr < args.min_correlation:
        raise RuntimeError(
            f"LD-score vs -log10(P) correlation too low: {corr:.4f} < {args.min_correlation:.4f}"
        )

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
