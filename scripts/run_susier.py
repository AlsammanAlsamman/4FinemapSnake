#!/usr/bin/env python3
import argparse
import json
import os

import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run SuSiE placeholder in Python")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--ld-matrix", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--coverage", type=float, required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    df = pd.read_csv(args.input_tsv, sep="\t")
    if "SNP" not in df.columns:
        raise ValueError("Matched table must include SNP")

    if all(col in df.columns for col in ["BETA", "SE"]):
        z = pd.to_numeric(df["BETA"], errors="coerce") / pd.to_numeric(df["SE"], errors="coerce")
        z = np.abs(z.replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float))
    else:
        p = pd.to_numeric(df.get("P", 1.0), errors="coerce").fillna(1.0).clip(lower=1e-300)
        z = -np.log10(p.to_numpy(dtype=float))

    order = np.argsort(-z)
    keep_n = min(args.L, len(order))
    keep_idx = order[:keep_n]
    sel = df.iloc[keep_idx].copy()

    z_sel = z[keep_idx]
    denom = float(np.sum(z_sel))
    sel["posterior"] = z_sel / denom if denom > 0 else 0.0

    sel.to_csv(args.out_tsv, sep="\t", index=False)

    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": "susier_placeholder",
                "L": int(args.L),
                "coverage": float(args.coverage),
                "n_selected": int(sel.shape[0]),
            },
            handle,
            indent=2,
        )

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
