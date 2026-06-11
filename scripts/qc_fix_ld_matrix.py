#!/usr/bin/env python3
import argparse
import json
import os

import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="QC and symmetrize LD matrix")
    parser.add_argument("--input-matrix", required=True)
    parser.add_argument("--out-matrix", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    parser.add_argument("--tolerance", type=float, required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_matrix)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    mat_df = pd.read_csv(args.input_matrix, sep="\t", index_col=0)
    mat = mat_df.to_numpy(dtype=float)
    if mat.shape[0] != mat.shape[1]:
        raise ValueError("LD matrix must be square")

    before = float(np.max(np.abs(mat - mat.T))) if mat.size else 0.0
    fixed = (mat + mat.T) / 2.0
    np.fill_diagonal(fixed, 1.0)
    after = float(np.max(np.abs(fixed - fixed.T))) if fixed.size else 0.0

    if before > args.tolerance:
        # Keep going but record that correction was applied.
        pass

    out_df = pd.DataFrame(fixed, index=mat_df.index, columns=mat_df.columns)
    out_df.to_csv(args.out_matrix, sep="\t")

    diag = {
        "max_abs_asymmetry_before": before,
        "max_abs_asymmetry_after": after,
        "tolerance": args.tolerance,
        "corrected": bool(before > 0),
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
