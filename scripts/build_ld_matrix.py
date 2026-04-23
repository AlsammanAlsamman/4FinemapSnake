#!/usr/bin/env python3
import argparse
import json
import os

import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a deterministic LD-like matrix from matched SNP table")
    parser.add_argument("--input-tsv", required=True)
    parser.add_argument("--ref-bim", required=True)
    parser.add_argument("--out-matrix", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_matrix)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    df = pd.read_csv(args.input_tsv, sep="\t")
    if "SNP" not in df.columns:
        raise ValueError("Input file must include column SNP")

    # Enforce that Step 5 uses the matched PLINK subset produced in Step 3.
    # PLINK BIM format: CHR SNP CM BP A1 A2
    ref_bim = pd.read_csv(
        args.ref_bim,
        sep="\t",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        usecols=[0, 1, 2, 3, 4, 5],
    )
    ref_snps = ref_bim["SNP"].astype(str).tolist()
    gwas_snps = df["SNP"].astype(str).tolist()

    if len(ref_snps) == 0:
        raise ValueError("Reference BIM has zero SNPs")

    set_equal = set(ref_snps) == set(gwas_snps)
    if not set_equal:
        only_ref = sorted(set(ref_snps) - set(gwas_snps))[:5]
        only_gwas = sorted(set(gwas_snps) - set(ref_snps))[:5]
        raise ValueError(
            "SNP mismatch between matched GWAS TSV and refpanel_matched BIM. "
            f"Examples only_in_bim={only_ref}, only_in_gwas={only_gwas}"
        )

    # Reorder GWAS rows to PLINK subset order for deterministic LD matrix ordering.
    df = df.set_index(df["SNP"].astype(str)).loc[ref_snps].reset_index(drop=True)

    # Build an LD proxy from genomic distance when genotypes are unavailable.
    # This avoids the degenerate +/-1 structure from a single-vector outer product.
    if "BP" in df.columns:
        bp = pd.to_numeric(df["BP"], errors="coerce").to_numpy(dtype=float)
    else:
        bp = np.arange(df.shape[0], dtype=float)

    # Fill missing/invalid positions with ordered indices to keep matrix computable.
    bad_bp = ~np.isfinite(bp)
    if np.any(bad_bp):
        bp[bad_bp] = np.arange(df.shape[0], dtype=float)[bad_bp]

    if bp.size <= 1:
        mat = np.ones((bp.size, bp.size), dtype=float)
    else:
        span = float(np.nanmax(bp) - np.nanmin(bp))
        # Use a scale tied to locus span, with sensible floor to prevent over-shrinkage.
        tau = max(span / 20.0, 5000.0)
        dist = np.abs(bp[:, None] - bp[None, :])
        mat = np.exp(-dist / tau)
        np.fill_diagonal(mat, 1.0)

    snps = df["SNP"].astype(str).tolist()
    out = pd.DataFrame(mat, index=snps, columns=snps)
    out.to_csv(args.out_matrix, sep="\t")

    diag = {
        "n_snps": int(len(snps)),
        "n_ref_bim_snps": int(len(ref_snps)),
        "snp_set_equal_to_ref_bim": bool(set_equal),
        "matrix_shape": [int(out.shape[0]), int(out.shape[1])],
        "ld_model": "exp_decay_by_bp",
        "bp_span": float(np.nanmax(bp) - np.nanmin(bp)) if bp.size else 0.0,
        "tau": float(max((float(np.nanmax(bp) - np.nanmin(bp)) / 20.0), 5000.0)) if bp.size > 1 else 0.0,
        "max_abs_asymmetry": float(np.max(np.abs(mat - mat.T))) if mat.size else 0.0,
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
