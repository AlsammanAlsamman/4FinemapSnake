#!/usr/bin/env python3
import argparse
import os

import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _tag(df: pd.DataFrame, method: str) -> pd.DataFrame:
    out = df.copy()
    out["method"] = method
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Export per-locus summary table and Excel")
    parser.add_argument("--target", required=True)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--finemap-tsv", required=True)
    parser.add_argument("--susier-tsv", required=True)
    parser.add_argument("--cojo-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-xlsx", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.out_xlsx)
    _mkdir_for(args.done_file)

    finemap = _tag(pd.read_csv(args.finemap_tsv, sep="\t"), "finemap")
    susier = _tag(pd.read_csv(args.susier_tsv, sep="\t"), "susier")
    cojo = _tag(pd.read_csv(args.cojo_tsv, sep="\t"), "cojo")

    summary = pd.concat([finemap, susier, cojo], ignore_index=True, sort=False)
    summary.insert(0, "target", args.target)
    summary.insert(1, "locus", args.locus)

    summary.to_csv(args.out_tsv, sep="\t", index=False)
    with pd.ExcelWriter(args.out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
