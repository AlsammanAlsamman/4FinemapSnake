#!/usr/bin/env python3
"""
export_snp_master_table.py
Produce four output files per locus:
  snp_master_table.tsv  — one row per SNP with GWAS stats, fine-mapping PIPs,
                          LD score, allele frequencies and per-dataset eQTL columns.
  ld_matrix_r2.tsv      — pairwise r² for credible-set SNPs only.
  ld_matrix_r.tsv       — signed Pearson r for credible-set SNPs only.
  ld_matrix_Dprime.tsv  — D' for credible-set SNPs only.

All numeric values are written at full floating-point precision.
Missing values are written as NA; no cells are left empty.
"""

import argparse
import os
import struct

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _ci(col: str, df: pd.DataFrame) -> str | None:
    """Case-insensitive column lookup; return actual column name or None."""
    want = col.lower()
    for c in df.columns:
        if str(c).lower() == want:
            return str(c)
    return None


def _ci_any(candidates: list[str], df: pd.DataFrame) -> str | None:
    """Return first matching column from a priority list."""
    for name in candidates:
        found = _ci(name, df)
        if found is not None:
            return found
    return None


# ---------------------------------------------------------------------------
# PLINK BED reader  (SNP-major order, v1 magic)
# ---------------------------------------------------------------------------
# Encoding (2 bits per sample, LSB-first within each byte):
#   00 → hom-A1  (dosage 0)
#   01 → missing (dosage NaN)
#   10 → het     (dosage 1)
#   11 → hom-A2  (dosage 2)
_BED_DECODE = np.array([0.0, np.nan, 1.0, 2.0])


def _read_bed(bed_path: str, bim_df: pd.DataFrame, fam_df: pd.DataFrame,
              snp_subset: set | None = None) -> tuple[np.ndarray, list[str]]:
    """
    Read a PLINK BED file (SNP-major) and return:
      (dosage_matrix, snp_ids)
    dosage_matrix shape: (n_snps_returned, n_samples)
    Dosage = number of A2 alleles (0 / 1 / 2 / NaN for missing).
    snp_ids: list of SNP IDs in the same row order as dosage_matrix.
    """
    n_samples = len(fam_df)
    n_per_row = (n_samples + 3) // 4  # bytes per SNP row

    all_snps = bim_df["SNP"].astype(str).tolist()
    if snp_subset is None:
        indices = list(range(len(all_snps)))
        out_snps = all_snps
    else:
        indices = [i for i, s in enumerate(all_snps) if s in snp_subset]
        out_snps = [all_snps[i] for i in indices]

    if not indices:
        return np.empty((0, n_samples), dtype=float), []

    dosage = np.full((len(indices), n_samples), np.nan, dtype=float)

    with open(bed_path, "rb") as fh:
        magic = fh.read(3)
        if magic != bytes([0x6C, 0x1B, 0x01]):
            raise ValueError(f"Invalid PLINK BED magic bytes: {magic.hex()}")

        for out_i, snp_i in enumerate(indices):
            fh.seek(3 + snp_i * n_per_row)
            raw = np.frombuffer(fh.read(n_per_row), dtype=np.uint8)
            unpacked = np.empty(n_per_row * 4, dtype=float)
            for shift in range(4):
                bits = (raw >> (2 * shift)) & 3
                unpacked[shift::4] = _BED_DECODE[bits]
            dosage[out_i] = unpacked[:n_samples]

    return dosage, out_snps


# ---------------------------------------------------------------------------
# LD computation
# ---------------------------------------------------------------------------

def _compute_ld(dosage: np.ndarray, snp_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    From a (n_snps × n_samples) dosage matrix, return three square DataFrames:
      r_df       — signed Pearson r
      r2_df      — r²
      dprime_df  — D′  (estimated from unphased dosage data)
    All diagonals = 1.  Missing genotypes are dropped pairwise.
    """
    n = len(snp_ids)

    r_mat = np.ones((n, n), dtype=float)
    r2_mat = np.ones((n, n), dtype=float)
    dp_mat = np.ones((n, n), dtype=float)

    for i in range(n):
        gi = dosage[i]
        pi = np.nanmean(gi) / 2.0          # freq of A2 allele

        for j in range(i + 1, n):
            gj = dosage[j]
            pj = np.nanmean(gj) / 2.0

            # Pairwise complete observations
            mask = np.isfinite(gi) & np.isfinite(gj)
            n_obs = int(mask.sum())

            if n_obs < 2 or pi <= 0.0 or pi >= 1.0 or pj <= 0.0 or pj >= 1.0:
                r_val = 0.0
                r2_val = 0.0
                dp_val = 0.0
            else:
                gi_c = gi[mask] - gi[mask].mean()
                gj_c = gj[mask] - gj[mask].mean()
                num = float(np.dot(gi_c, gj_c))
                denom = float(np.sqrt(np.dot(gi_c, gi_c) * np.dot(gj_c, gj_c)))
                r_val = num / denom if denom > 0.0 else 0.0
                r_val = float(np.clip(r_val, -1.0, 1.0))
                r2_val = r_val * r_val

                # D′ estimate via r and allele frequencies
                # D = r * sqrt(p_i*(1-p_i) * p_j*(1-p_j))
                pi_obs = float(np.nanmean(dosage[i][mask]) / 2.0)
                pj_obs = float(np.nanmean(dosage[j][mask]) / 2.0)
                D = r_val * float(np.sqrt(pi_obs * (1.0 - pi_obs) * pj_obs * (1.0 - pj_obs)))
                if D > 0.0:
                    D_max = min(pi_obs * (1.0 - pj_obs), pj_obs * (1.0 - pi_obs))
                elif D < 0.0:
                    D_max = min(pi_obs * pj_obs, (1.0 - pi_obs) * (1.0 - pj_obs))
                else:
                    D_max = 0.0
                dp_val = float(np.clip(D / D_max, -1.0, 1.0)) if D_max > 0.0 else 0.0

            r_mat[i, j] = r_mat[j, i] = r_val
            r2_mat[i, j] = r2_mat[j, i] = r2_val
            dp_mat[i, j] = dp_mat[j, i] = dp_val

    r_df = pd.DataFrame(r_mat, index=snp_ids, columns=snp_ids)
    r2_df = pd.DataFrame(r2_mat, index=snp_ids, columns=snp_ids)
    dp_df = pd.DataFrame(dp_mat, index=snp_ids, columns=snp_ids)

    r_df.index.name = "SNP_rsID"
    r2_df.index.name = "SNP_rsID"
    dp_df.index.name = "SNP_rsID"

    return r_df, r2_df, dp_df


# ---------------------------------------------------------------------------
# write_tsv: always NA for missing, full precision, no rounding
# ---------------------------------------------------------------------------

def _write_tsv(df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame to a tab-separated file with NA for missing values."""
    _mkdir_for(path)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format=None)


def _write_square_tsv(df: pd.DataFrame, path: str) -> None:
    """Write a square LD matrix (index = SNP_rsID) to TSV with header row."""
    _mkdir_for(path)
    df.to_csv(path, sep="\t", index=True, na_rep="NA", float_format=None)


# ---------------------------------------------------------------------------
# eQTL helpers
# ---------------------------------------------------------------------------

_BETA_CANDIDATES  = ["beta_eQTL", "beta", "b", "effect_size", "effect", "slope"]
_PVAL_CANDIDATES  = ["pval_eQTL", "pval", "p_value", "pvalue", "p.value", "p", "nom_pval", "pval_nominal"]
_RSID_CANDIDATES  = ["rsid", "snp", "snpid", "markername", "varid", "variant_id", "rs_id"]


def _pivot_eqtls(combined: pd.DataFrame, snp_set: set[str]) -> dict[str, dict[str, dict]]:
    """
    Returns nested dict:
      { dataset_name: { rsid: { "beta": v, "pval": v } } }
    """
    if combined.empty:
        return {}

    result: dict[str, dict[str, dict]] = {}

    # Ensure source_dataset column exists
    if "source_dataset" not in combined.columns:
        return {}

    for ds, grp in combined.groupby("source_dataset", sort=False):
        ds_str = str(ds)
        # Prefer matched_rsid because it is aligned to summary SNP IDs by step 14.
        rsid_col = "matched_rsid" if "matched_rsid" in grp.columns else None
        if rsid_col is None:
            rsid_col = _ci_any(_RSID_CANDIDATES, grp)
        if rsid_col is None:
            continue

        beta_col = _ci_any(_BETA_CANDIDATES, grp)
        pval_col = _ci_any(_PVAL_CANDIDATES, grp)

        grp2 = grp.copy()
        grp2["_snp"] = grp2[rsid_col].astype(str)
        grp2 = grp2[grp2["_snp"].isin(snp_set)]
        if grp2.empty:
            result[ds_str] = {}
            continue

        grp2["_beta"] = pd.to_numeric(grp2[beta_col], errors="coerce") if beta_col else np.nan
        grp2["_pval"] = pd.to_numeric(grp2[pval_col], errors="coerce") if pval_col else np.nan
        grp2["_abs_beta"] = grp2["_beta"].abs()

        # Choose one row per SNP per dataset: smallest p-value, then largest |beta|.
        grp2 = grp2.sort_values(["_snp", "_pval", "_abs_beta"], ascending=[True, True, False], na_position="last")
        best = grp2.drop_duplicates(subset=["_snp"], keep="first")

        per_snp: dict[str, dict] = {}
        for _, row in best.iterrows():
            snp = str(row["_snp"])
            per_snp[snp] = {
                "beta": row["_beta"] if pd.notna(row["_beta"]) else pd.NA,
                "pval": row["_pval"] if pd.notna(row["_pval"]) else pd.NA,
            }
        result[ds_str] = per_snp

    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Export SNP master table and LD matrices")
    parser.add_argument("--matched-tsv",   required=True)
    parser.add_argument("--afreq",         required=True)
    parser.add_argument("--bed",           required=True)
    parser.add_argument("--bim",           required=True)
    parser.add_argument("--fam",           required=True)
    parser.add_argument("--finemap-tsv",   required=True)
    parser.add_argument("--susier-tsv",    required=True)
    parser.add_argument("--eqtl-combined", required=True)
    parser.add_argument("--eqtl-manifest", required=True)
    parser.add_argument("--out-master",    required=True)
    parser.add_argument("--out-ld-r2",     required=True)
    parser.add_argument("--out-ld-r",      required=True)
    parser.add_argument("--out-ld-dprime", required=True)
    parser.add_argument("--done-file",     required=True)
    args = parser.parse_args()

    for p in [args.out_master, args.out_ld_r2, args.out_ld_r, args.out_ld_dprime, args.done_file]:
        _mkdir_for(p)

    # ------------------------------------------------------------------
    # 1. Load GWAS matched table
    # ------------------------------------------------------------------
    matched = pd.read_csv(args.matched_tsv, sep="\t")
    matched["SNP"] = matched["SNP"].astype(str)
    matched["A1"] = matched["A1"].astype(str).str.upper()
    matched["A2"] = matched["A2"].astype(str).str.upper()

    p_col   = _ci_any(["p", "pvalue", "p_value"], matched)
    beta_col = _ci_any(["beta", "b"], matched)
    se_col   = _ci_any(["se"], matched)
    chr_col  = _ci_any(["chr", "chrom", "chromosome"], matched)
    bp_col   = _ci_any(["bp", "pos", "position"], matched)

    if p_col is None or beta_col is None or se_col is None:
        raise ValueError("matched.tsv must contain P, BETA, SE columns")

    matched["_P"]    = pd.to_numeric(matched[p_col],    errors="coerce")
    matched["_BETA"] = pd.to_numeric(matched[beta_col], errors="coerce")
    matched["_SE"]   = pd.to_numeric(matched[se_col],   errors="coerce")

    # GWAS derived stats
    matched["_OR"]        = np.exp(matched["_BETA"])
    matched["_log10p"]    = -np.log10(matched["_P"].clip(lower=np.finfo(float).tiny))
    matched["_chi2"]      = (matched["_BETA"] / matched["_SE"]) ** 2
    matched["_ci95_low"]  = np.exp(matched["_BETA"] - 1.96 * matched["_SE"])
    matched["_ci95_high"] = np.exp(matched["_BETA"] + 1.96 * matched["_SE"])

    # ------------------------------------------------------------------
    # 2. Load allele frequencies
    # ------------------------------------------------------------------
    afreq = pd.read_csv(args.afreq, sep="\t")
    id_col  = _ci_any(["id", "#id"], afreq)
    ref_col = _ci_any(["ref"], afreq)
    alt_col = _ci_any(["alt"], afreq)
    af_col  = _ci_any(["alt_freqs", "alt_freq", "af", "freq"], afreq)

    if id_col:
        afreq = afreq.rename(columns={id_col: "ID"})
    if ref_col:
        afreq = afreq.rename(columns={ref_col: "REF"})
    if alt_col:
        afreq = afreq.rename(columns={alt_col: "ALT"})
    if af_col:
        afreq = afreq.rename(columns={af_col: "ALT_FREQS"})

    afreq["ID"] = afreq["ID"].astype(str)
    afreq["REF"] = afreq["REF"].astype(str).str.upper()
    afreq["ALT"] = afreq["ALT"].astype(str).str.upper()

    matched = matched.merge(
        afreq[["ID", "REF", "ALT", "ALT_FREQS"]],
        left_on="SNP", right_on="ID", how="left",
    )

    def _af_for_allele(row, allele_col):
        a = row[allele_col].upper()
        alt = str(row.get("ALT", "")).upper()
        ref = str(row.get("REF", "")).upper()
        af  = pd.to_numeric(row.get("ALT_FREQS"), errors="coerce")
        if pd.isna(af):
            return pd.NA
        if a == alt:
            return float(af)
        if a == ref:
            return float(1.0 - af)
        return pd.NA

    matched["AF_EA"] = matched.apply(lambda r: _af_for_allele(r, "A1"), axis=1)
    matched["AF_OA"] = matched.apply(lambda r: _af_for_allele(r, "A2"), axis=1)

    # ------------------------------------------------------------------
    # 3. Load FINEMAP credible set
    # ------------------------------------------------------------------
    finemap = pd.read_csv(args.finemap_tsv, sep="\t")
    finemap["SNP"] = finemap["SNP"].astype(str)
    fp_col = _ci_any(["pip", "prob"], finemap)
    fc_col = _ci_any(["cs", "credibleset", "cs_id"], finemap)

    if fp_col:
        finemap["_fm_pip"] = pd.to_numeric(finemap[fp_col], errors="coerce")
    else:
        finemap["_fm_pip"] = pd.NA

    finemap["_fm_cs"] = 1   # all rows in credible-set file are in the CS
    if fc_col:
        finemap["_fm_cs"] = finemap[fc_col].where(finemap[fc_col].notna(), 1)

    fm_keep = finemap[["SNP", "_fm_pip", "_fm_cs"]].drop_duplicates("SNP", keep="first")

    # ------------------------------------------------------------------
    # 4. Load SuSiE credible set
    # ------------------------------------------------------------------
    susier = pd.read_csv(args.susier_tsv, sep="\t")
    susier["SNP"] = susier["SNP"].astype(str)
    sp_col = _ci_any(["pip"], susier)
    sc_col = _ci_any(["cs", "credibleset", "cs_id", "cs_index"], susier)

    if sp_col:
        susier["_su_pip"] = pd.to_numeric(susier[sp_col], errors="coerce")
    else:
        susier["_su_pip"] = pd.NA

    if sc_col:
        susier["_su_cs"] = susier[sc_col]
    else:
        susier["_su_cs"] = pd.NA

    su_keep = susier[["SNP", "_su_pip", "_su_cs"]].drop_duplicates("SNP", keep="first")

    # ------------------------------------------------------------------
    # 5. Merge fine-mapping into base table
    # ------------------------------------------------------------------
    base = matched.merge(fm_keep, on="SNP", how="left")
    base = base.merge(su_keep,  on="SNP", how="left")

    # Composite PIP = max(FINEMAP_PIP, SuSiE_PIP), treating NA as 0
    def _composite(row):
        fm  = pd.to_numeric(row["_fm_pip"], errors="coerce")
        su  = pd.to_numeric(row["_su_pip"], errors="coerce")
        vals = [v for v in [fm, su] if pd.notna(v)]
        return float(max(vals)) if vals else pd.NA

    base["_composite_pip"] = base.apply(_composite, axis=1)

    # CS membership: prefer SuSiE CS label, fall back to FINEMAP indicator
    def _cs_membership(row):
        su_cs = row.get("_su_cs")
        fm_cs = row.get("_fm_cs")
        if pd.notna(su_cs) and str(su_cs).strip() not in ("", "nan"):
            return str(su_cs)
        if pd.notna(fm_cs) and str(fm_cs).strip() not in ("", "nan"):
            return str(fm_cs)
        return pd.NA

    base["_cs_membership"] = base.apply(_cs_membership, axis=1)

    # CS size: count SNPs per CS label
    cs_sizes: dict[str, int] = {}
    for label, grp in base.dropna(subset=["_cs_membership"]).groupby("_cs_membership"):
        cs_sizes[str(label)] = len(grp)

    def _cs_size(row):
        m = row["_cs_membership"]
        if pd.isna(m):
            return pd.NA
        return cs_sizes.get(str(m), pd.NA)

    base["_cs_size"] = base.apply(_cs_size, axis=1)

    # ------------------------------------------------------------------
    # 6. Compute LD from PLINK BED/BIM/FAM
    # ------------------------------------------------------------------
    bim_df = pd.read_csv(
        args.bim, sep="\t", header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    )
    fam_df = pd.read_csv(
        args.fam, sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHEN"],
    )

    all_locus_snps = set(bim_df["SNP"].astype(str))
    dosage_all, snp_ids_all = _read_bed(args.bed, bim_df, fam_df, snp_subset=all_locus_snps)

    # LD score: Σ r² over all locus SNPs (standard LDSC definition, including self → +1)
    ld_score_map: dict[str, float] = {}
    if dosage_all.shape[0] > 0:
        # Compute r² row-sums efficiently
        # Standardise each row (mean-center, unit-variance), ignoring NaN pairwise
        n_all = dosage_all.shape[0]
        # We compute r for each pair only once; this is O(n²) but n is typically ≤ a few thousand
        for i, snp_i in enumerate(snp_ids_all):
            gi = dosage_all[i]
            r2_sum = 1.0  # self r² = 1
            for j in range(n_all):
                if j == i:
                    continue
                gj = dosage_all[j]
                mask = np.isfinite(gi) & np.isfinite(gj)
                n_obs = int(mask.sum())
                if n_obs < 2:
                    continue
                gi_c = gi[mask] - gi[mask].mean()
                gj_c = gj[mask] - gj[mask].mean()
                denom = float(np.sqrt(np.dot(gi_c, gi_c) * np.dot(gj_c, gj_c)))
                if denom > 0.0:
                    r = float(np.dot(gi_c, gj_c)) / denom
                    r2_sum += r * r
            ld_score_map[snp_i] = r2_sum

    base["_ld_score"] = base["SNP"].map(ld_score_map)

    # ------------------------------------------------------------------
    # 7. Credible-set SNP LD matrices (r, r², D′)
    # ------------------------------------------------------------------
    cs_snps = base.loc[base["_cs_membership"].notna(), "SNP"].astype(str).unique().tolist()

    if cs_snps:
        cs_set = set(cs_snps)
        dosage_cs, snp_ids_cs = _read_bed(args.bed, bim_df, fam_df, snp_subset=cs_set)
        r_df, r2_df, dp_df = _compute_ld(dosage_cs, snp_ids_cs)
    else:
        # No credible set — produce empty 0×0 matrices with a stub header
        empty_idx: list[str] = []
        r_df  = pd.DataFrame(index=empty_idx, columns=empty_idx)
        r2_df = pd.DataFrame(index=empty_idx, columns=empty_idx)
        dp_df = pd.DataFrame(index=empty_idx, columns=empty_idx)
        r_df.index.name = r2_df.index.name = dp_df.index.name = "SNP_rsID"

    # ------------------------------------------------------------------
    # 8. eQTL columns
    # ------------------------------------------------------------------
    try:
        eqtl_combined = pd.read_csv(args.eqtl_combined, sep="\t", low_memory=False)
    except Exception:
        eqtl_combined = pd.DataFrame()

    try:
        eqtl_manifest = pd.read_csv(args.eqtl_manifest, sep="\t")
        datasets = [str(d) for d in eqtl_manifest["dataset"].tolist() if pd.notna(d)]
    except Exception:
        datasets = []

    snp_set_all = set(base["SNP"].astype(str))
    eqtl_by_ds = _pivot_eqtls(eqtl_combined, snp_set_all)

    # Build eQTL columns for each dataset in manifest order
    eqtl_col_defs: list[tuple[str, str]] = []   # (col_name, attr)
    for ds in datasets:
        for attr in ["beta_eQTL", "pval_eQTL", "log10p_eQTL", "OR_eQTL"]:
            eqtl_col_defs.append((f"{attr}@{ds}", ds))

    def _safe_log10p(pval):
        v = pd.to_numeric(pval, errors="coerce")
        if pd.isna(v) or v <= 0.0:
            return pd.NA
        return float(-np.log10(v))

    eqtl_cols_data: dict[str, pd.Series] = {}
    for col_name, ds in eqtl_col_defs:
        attr = col_name.split("@")[0]  # e.g. "beta_eQTL"
        ds_data = eqtl_by_ds.get(ds, {})

        def _get_eqtl(snp, _ds_data=ds_data, _attr=attr):
            rec = _ds_data.get(str(snp), {})
            if _attr == "beta_eQTL":
                return rec.get("beta", pd.NA)
            if _attr == "pval_eQTL":
                return rec.get("pval", pd.NA)
            if _attr == "log10p_eQTL":
                return _safe_log10p(rec.get("pval", pd.NA))
            if _attr == "OR_eQTL":
                beta = pd.to_numeric(rec.get("beta", pd.NA), errors="coerce")
                return float(np.exp(beta)) if pd.notna(beta) else pd.NA
            return pd.NA

        eqtl_cols_data[col_name] = base["SNP"].map(_get_eqtl)

    if eqtl_cols_data:
        base = pd.concat([base, pd.DataFrame(eqtl_cols_data, index=base.index)], axis=1)

    # ------------------------------------------------------------------
    # 9. Assemble final master table columns
    # ------------------------------------------------------------------
    chr_val = chr_col if chr_col else "CHR"
    bp_val  = bp_col  if bp_col  else "BP"

    master_core = {
        "SNP_rsID": base["SNP"],
        "CHR": base[chr_val] if chr_val in base.columns else pd.Series(pd.NA, index=base.index),
        "POS": base[bp_val] if bp_val in base.columns else pd.Series(pd.NA, index=base.index),
        "EA": base["A1"],
        "OA": base["A2"],
        "GWAS_pval": base["_P"],
        "GWAS_OR": base["_OR"],
        "GWAS_beta": base["_BETA"],
        "GWAS_log10p": base["_log10p"],
        "GWAS_chi2": base["_chi2"],
        "GWAS_SE": base["_SE"],
        "GWAS_CI95_lower": base["_ci95_low"],
        "GWAS_CI95_upper": base["_ci95_high"],
        "FINEMAP_PIP": base["_fm_pip"],
        "SuSiE_PIP": base["_su_pip"],
        "Composite_PIP": base["_composite_pip"],
        "CS_membership": base["_cs_membership"],
        "CS_size": base["_cs_size"],
        "LD_score": base["_ld_score"],
        "AF_EA": base["AF_EA"],
        "AF_OA": base["AF_OA"],
    }
    master = pd.DataFrame(master_core, index=base.index)

    # Append eQTL columns in dataset order in one operation.
    eqtl_cols_order = [col_name for col_name, _ in eqtl_col_defs if col_name in base.columns]
    if eqtl_cols_order:
        master = pd.concat([master, base[eqtl_cols_order]], axis=1)

    # ------------------------------------------------------------------
    # 10. Sort: fine-mapped rows by Composite_PIP desc, remainder by P asc
    # ------------------------------------------------------------------
    has_pip  = master["Composite_PIP"].notna()
    finemaped = master.loc[has_pip].sort_values("Composite_PIP", ascending=False)
    not_fmap  = master.loc[~has_pip].sort_values("GWAS_pval", ascending=True)
    master = pd.concat([finemaped, not_fmap], ignore_index=True)

    # Ensure NA (not empty string) for all missing cells
    master = master.where(master.notna(), other=pd.NA)

    # ------------------------------------------------------------------
    # 11. Write outputs
    # ------------------------------------------------------------------
    _write_tsv(master, args.out_master)
    _write_square_tsv(r2_df, args.out_ld_r2)
    _write_square_tsv(r_df,  args.out_ld_r)
    _write_square_tsv(dp_df, args.out_ld_dprime)

    with open(args.done_file, "w", encoding="utf-8") as fh:
        fh.write("ok\n")


if __name__ == "__main__":
    main()
