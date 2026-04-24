#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _find_column_case_insensitive(df: pd.DataFrame, name: str) -> str | None:
    want = name.lower()
    for col in df.columns:
        if str(col).lower() == want:
            return str(col)
    return None


def _compute_or(beta):
    """Compute odds ratio from beta (log-odds)."""
    return np.exp(beta)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export per-locus summary table and Excel")
    parser.add_argument("--target", required=True)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--matched-tsv", required=True, help="Matched GWAS TSV with SNP, CHR, BP, A1, A2, P, BETA, SE")
    parser.add_argument("--afreq", required=True, help="Allele frequency file with ID, REF, ALT, ALT_FREQS")
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

    # Load base data: matched GWAS and allele frequencies
    matched = pd.read_csv(args.matched_tsv, sep="\t")
    afreq = pd.read_csv(args.afreq, sep="\t")

    # Standardize allele columns
    matched["A1"] = matched["A1"].str.upper()
    matched["A2"] = matched["A2"].str.upper()
    afreq["REF"] = afreq["REF"].str.upper()
    afreq["ALT"] = afreq["ALT"].str.upper()

    # Merge allele frequency into matched
    matched = matched.merge(
        afreq[["ID", "REF", "ALT", "ALT_FREQS"]],
        left_on="SNP", right_on="ID", how="left"
    )

    # Compute EA (effect allele = A1) and OA (other allele)
    matched["EA"] = matched["A1"]
    matched["OA"] = matched["A2"]
    
    # Compute AF for effect allele (A1)
    matched["AF"] = matched.apply(
        lambda row: (
            row["ALT_FREQS"]
            if row["A1"] == row["ALT"]
            else (1 - row["ALT_FREQS"] if row["A1"] == row["REF"] else pd.NA)
        ),
        axis=1
    )

    # Compute OR from BETA
    matched["OR"] = matched["BETA"].apply(_compute_or)

    # Create base table with essential columns
    base = matched[["SNP", "CHR", "BP", "EA", "OA", "AF", "P", "BETA", "OR", "SE"]].copy()

    # Load FINEMAP
    finemap = pd.read_csv(args.finemap_tsv, sep="\t")
    finemap["SNP"] = finemap["SNP"].astype(str)
    
    # Find PIP column case-insensitively
    pip_col = _find_column_case_insensitive(finemap, "PIP")
    if pip_col:
        finemap["finemap_pip"] = pd.to_numeric(finemap[pip_col], errors="coerce")
    else:
        finemap["finemap_pip"] = pd.NA

    # Find CS column case-insensitively
    cs_col = _find_column_case_insensitive(finemap, "CS")
    if cs_col:
        finemap["finemap_cs"] = finemap[cs_col]
    else:
        finemap["finemap_cs"] = pd.NA

    finemap_keep = finemap[["SNP", "finemap_pip", "finemap_cs"]].drop_duplicates(subset=["SNP"], keep="first")

    # Load SuSiE and filter for credible set members
    susier_raw = pd.read_csv(args.susier_tsv, sep="\t")
    susier_raw["SNP"] = susier_raw["SNP"].astype(str)
    
    pip_col_su = _find_column_case_insensitive(susier_raw, "PIP")
    cs_col_su = _find_column_case_insensitive(susier_raw, "CS")

    if pip_col_su is None:
        raise ValueError("SuSiE summary must include a PIP column")
    if cs_col_su is None:
        raise ValueError("SuSiE summary must include a CS column for credible set membership")

    susier_raw[pip_col_su] = pd.to_numeric(susier_raw[pip_col_su], errors="coerce")
    cs_is_present = susier_raw[cs_col_su].notna() & (susier_raw[cs_col_su].astype(str).str.strip() != "")
    susier_filtered = susier_raw.loc[cs_is_present & (susier_raw[pip_col_su] > 0)].copy()

    susier_filtered["susie_pip"] = susier_filtered[pip_col_su]
    susier_filtered["susie_cs"] = susier_filtered[cs_col_su]
    susier_keep = susier_filtered[["SNP", "susie_pip", "susie_cs"]].drop_duplicates(subset=["SNP"], keep="first")

    # Load COJO
    cojo_raw = pd.read_csv(args.cojo_tsv, sep="\t")
    if len(cojo_raw) > 0:
        cojo_raw["SNP"] = cojo_raw["lead_snp"].astype(str)
        cojo_raw["cojo"] = "iter" + cojo_raw["iteration"].astype(str)
        cojo_keep = cojo_raw[["SNP", "cojo"]].drop_duplicates(subset=["SNP"], keep="first")
    else:
        cojo_keep = pd.DataFrame(columns=["SNP", "cojo"])

    # Merge all data on SNP
    summary = base.copy()
    summary = summary.merge(finemap_keep, on="SNP", how="left")
    summary = summary.merge(susier_keep, on="SNP", how="left")
    summary = summary.merge(cojo_keep, on="SNP", how="left")

    # Fill NAs with empty string for better display
    for col in ["finemap_pip", "finemap_cs", "susie_pip", "susie_cs", "cojo"]:
        if col in summary.columns:
            summary[col] = summary[col].fillna("")

    # Order by P value
    summary = summary.sort_values("P").reset_index(drop=True)

    # Write outputs
    summary.to_csv(args.out_tsv, sep="\t", index=False)
    with pd.ExcelWriter(args.out_xlsx, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
