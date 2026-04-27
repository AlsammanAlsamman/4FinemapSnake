#!/usr/bin/env python3
import argparse
import os
import glob
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


def _load_cojo_iterations(cojo_dir: str) -> dict:
    """
    Load all COJO iteration files and extract conditional p-values and ORs.
    Returns dict: {iter_num: {SNP: {p: pC, or: OR}, ...}, ...}
    """
    iterations = {}
    
    # Find all iter_*.cma.cojo files
    iter_files = sorted(glob.glob(os.path.join(cojo_dir, "iter_*.cma.cojo")))
    
    for iter_file in iter_files:
        # Extract iteration number from filename (e.g., iter_01.cma.cojo -> 01)
        basename = os.path.basename(iter_file)
        iter_num_str = basename.split("_")[1].split(".")[0]  # Get "01" from "iter_01.cma.cojo"
        iter_num = int(iter_num_str)
        
        # Read the COJO output file
        try:
            df = pd.read_csv(iter_file, sep="\t")
            df["SNP"] = df["SNP"].astype(str)
            
            # Extract conditional p-value and beta
            pc_col = _find_column_case_insensitive(df, "pC")
            bc_col = _find_column_case_insensitive(df, "bC")
            
            if pc_col and bc_col:
                # Compute OR from conditional beta
                df["or_cond"] = df[bc_col].apply(_compute_or)
                
                # Create dictionary for this iteration
                iter_data = {}
                for _, row in df.iterrows():
                    snp = row["SNP"]
                    pc = pd.to_numeric(row[pc_col], errors="coerce")
                    or_cond = pd.to_numeric(row["or_cond"], errors="coerce")
                    iter_data[snp] = {"p": pc, "or": or_cond}
                
                iterations[iter_num] = iter_data
        except Exception as e:
            print(f"Warning: Could not read {iter_file}: {e}")
    
    return iterations


def main() -> None:
    parser = argparse.ArgumentParser(description="Export per-locus summary table and Excel")
    parser.add_argument("--target", required=True)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--matched-tsv", required=True, help="Matched GWAS TSV with SNP, CHR, BP, A1, A2, P, BETA, SE")
    parser.add_argument("--afreq", required=True, help="Allele frequency file with ID, REF, ALT, ALT_FREQS")
    parser.add_argument("--finemap-tsv", required=True)
    parser.add_argument("--susier-tsv", required=True)
    parser.add_argument("--cojo-tsv", required=True)
    parser.add_argument("--cojo-dir", required=True, help="Directory with COJO iteration files (gcta_tmp)")
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-xlsx", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.out_xlsx)
    _mkdir_for(args.done_file)

    # Load COJO iterations
    cojo_iterations = _load_cojo_iterations(args.cojo_dir)

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

    # FINEMAP credible set TSV contains only SNPs in credible sets, so finemap_cs = 1 for all
    finemap["finemap_cs"] = 1

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

    # Add COJO conditional p-values and ORs for each iteration
    for iter_num in sorted(cojo_iterations.keys()):
        iter_data = cojo_iterations[iter_num]
        col_p = f"p_iter{iter_num}"
        col_or = f"or_iter{iter_num}"
        
        summary[col_p] = summary["SNP"].map(lambda snp: iter_data.get(snp, {}).get("p", np.nan))
        summary[col_or] = summary["SNP"].map(lambda snp: iter_data.get(snp, {}).get("or", np.nan))

    # Fill NAs with empty string for better display for specific columns
    for col in ["finemap_pip", "finemap_cs", "susie_pip", "susie_cs", "cojo"]:
        if col in summary.columns:
            summary[col] = summary[col].fillna("")

    # Reorder columns to match desired output: 
    # SNP, CHR, BP, EA, OA, AF, P, BETA, OR, SE, cojo, finemap_pip, finemap_cs, susie_pip, susie_cs, [p_iterN, or_iterN, ...]
    desired_cols = ["SNP", "CHR", "BP", "EA", "OA", "AF", "P", "BETA", "OR", "SE", "cojo", "finemap_pip", "finemap_cs", "susie_pip", "susie_cs"]
    
    # Add iteration columns at the end
    iter_cols = sorted([c for c in summary.columns if c.startswith("p_iter") or c.startswith("or_iter")])
    desired_cols.extend(iter_cols)
    
    available_cols = [c for c in desired_cols if c in summary.columns]
    summary = summary[available_cols]

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
