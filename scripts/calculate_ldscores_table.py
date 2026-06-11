#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _load_ld_matrix(path: str) -> pd.DataFrame:
    ld = pd.read_csv(path, sep="\t", index_col=0)
    ld.index = ld.index.astype(str)
    ld.columns = ld.columns.astype(str)
    ld = ld.apply(pd.to_numeric, errors="coerce")
    return ld


def _as_r2(ld: pd.DataFrame) -> pd.DataFrame:
    mat = ld.to_numpy(dtype=float)
    if mat.size == 0:
        return ld

    n = mat.shape[0]
    off_mask = ~np.eye(n, dtype=bool)
    offdiag_vals = mat[off_mask]
    diag_vals = np.diag(mat)

    has_negatives = np.any(offdiag_vals < -1e-6)
    diag_finite = diag_vals[np.isfinite(diag_vals)]
    diag_all_one = bool(diag_finite.size > 0 and np.all(np.abs(diag_finite - 1.0) < 1e-6))
    offdiag_abs_max = float(np.nanmax(np.abs(offdiag_vals))) if offdiag_vals.size else 0.0
    is_symmetric = bool(np.allclose(mat, mat.T, equal_nan=True, atol=1e-6, rtol=0.0))

    already_r2 = (not has_negatives) and diag_all_one and (offdiag_abs_max <= 1.0 + 1e-6) and is_symmetric

    if not already_r2:
        mat = np.square(mat)

    return pd.DataFrame(mat, index=ld.index, columns=ld.columns)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute all-SNP LD score table and normalize by number of SNPs in locus"
    )
    parser.add_argument("--matched-tsv", required=True)
    parser.add_argument("--ld-matrix", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.output_tsv)
    _mkdir_for(args.done_file)

    matched = pd.read_csv(args.matched_tsv, sep="\t", low_memory=False)
    if "SNP" not in matched.columns:
        raise ValueError("matched.tsv must contain a SNP column")
    matched["SNP"] = matched["SNP"].astype(str)
    matched_unique = matched.drop_duplicates(subset=["SNP"], keep="first").copy()

    ld = _load_ld_matrix(args.ld_matrix)
    shared_axis = [s for s in ld.index if s in set(ld.columns)]
    ld = ld.loc[shared_axis, shared_axis]
    ld = _as_r2(ld)

    matched_snps = matched_unique["SNP"].tolist()
    common_snps = [s for s in matched_snps if s in ld.index]
    if not common_snps:
        raise ValueError("No overlapping SNPs between matched.tsv and LD matrix")

    ld_sub = ld.loc[common_snps, common_snps]
    mat = ld_sub.to_numpy(dtype=float)
    diag_vals = np.diag(mat)

    ld_score_raw = np.nansum(mat, axis=1) - diag_vals
    ld_score_raw = np.clip(ld_score_raw, a_min=0.0, a_max=None)

    m = int(len(common_snps))
    denom_m_minus_1 = max(m - 1, 1)
    denom_m = max(m, 1)

    score_df = pd.DataFrame(
        {
            "SNP": common_snps,
            "LD_score_sum_r2_excl_self": ld_score_raw,
            "LD_score_norm_by_m_minus_1": ld_score_raw / float(denom_m_minus_1),
            "LD_score_norm_by_m": ld_score_raw / float(denom_m),
            "Locus_snp_count_m": m,
        }
    )

    out = matched_unique.merge(score_df, on="SNP", how="inner")

    sort_cols = [c for c in ["CHR", "BP"] if c in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols, kind="mergesort")

    out.to_csv(args.output_tsv, sep="\t", index=False)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()