#!/usr/bin/env python3
import argparse
import json
import os

import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _load_ref_bim(prefix: str) -> pd.DataFrame:
    path = f"{prefix}.bim"
    ref = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "REF_A1", "REF_A2"],
    )
    ref["SNP"] = ref["SNP"].astype(str)
    ref["CHR"] = ref["CHR"].astype(str).str.replace("chr", "", regex=False)
    ref["BP"] = pd.to_numeric(ref["BP"], errors="coerce")
    ref["REF_A1"] = ref["REF_A1"].astype(str).str.upper()
    ref["REF_A2"] = ref["REF_A2"].astype(str).str.upper()
    return ref[["SNP", "CHR", "BP", "REF_A1", "REF_A2"]]


def _harmonize_row(row):
    a1 = str(row["A1"]).upper()
    a2 = str(row["A2"]).upper()
    r1 = str(row["REF_A1"]).upper()
    r2 = str(row["REF_A2"]).upper()

    if a1 == r1 and a2 == r2:
        return row, False, "match"
    if a1 == r2 and a2 == r1:
        row["A1"], row["A2"] = r1, r2
        if "BETA" in row and pd.notna(row["BETA"]):
            row["BETA"] = -float(row["BETA"])
        return row, True, "flipped"
    return row, False, "dropped_allele_mismatch"


def _harmonize_matches(matches: pd.DataFrame, source: str):
    kept_rows = []
    flipped = 0
    dropped = 0

    for _, row in matches.iterrows():
        row2, was_flipped, status = _harmonize_row(row.copy())
        if status == "dropped_allele_mismatch":
            dropped += 1
            continue
        if source == "chr_bp" and "SNP_REF" in row2:
            row2["SNP"] = str(row2["SNP_REF"])
        if was_flipped:
            flipped += 1
        row2["_map_source"] = source
        kept_rows.append(row2)

    if kept_rows:
        out = pd.DataFrame(kept_rows)
    else:
        out = matches.head(0).copy()
    return out, flipped, dropped


def main() -> None:
    parser = argparse.ArgumentParser(description="Harmonize GWAS alleles to reference panel BIM alleles")
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--ref-prefix", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    gwas = pd.read_csv(args.gwas, sep="\t")
    gwas["_row_id"] = range(gwas.shape[0])
    gwas["SNP"] = gwas["SNP"].astype(str)
    gwas["CHR"] = gwas["CHR"].astype(str).str.replace("chr", "", regex=False)
    gwas["BP"] = pd.to_numeric(gwas["BP"], errors="coerce")
    gwas["A1"] = gwas["A1"].astype(str).str.upper()
    gwas["A2"] = gwas["A2"].astype(str).str.upper()
    ref = _load_ref_bim(args.ref_prefix)

    # 1) Match by existing SNP identifier where possible.
    merged_snp = gwas.merge(ref, on="SNP", how="inner", suffixes=("", "_REF"))
    mapped_snp, flipped_snp, dropped_snp = _harmonize_matches(merged_snp, "snp")
    mapped_snp_ids = set(mapped_snp["_row_id"].tolist()) if not mapped_snp.empty else set()

    # 2) For remaining variants, match by CHR/BP and assign the reference SNP (RSID).
    gwas_unmapped = gwas[~gwas["_row_id"].isin(mapped_snp_ids)].copy()
    ref_by_pos = ref.rename(columns={"SNP": "SNP_REF"})
    merged_chr_bp = gwas_unmapped.merge(ref_by_pos, on=["CHR", "BP"], how="inner")
    mapped_chr_bp, flipped_chr_bp, dropped_chr_bp = _harmonize_matches(merged_chr_bp, "chr_bp")

    if not mapped_chr_bp.empty:
        # Resolve any ambiguous multi-hit positions deterministically.
        mapped_chr_bp = mapped_chr_bp.sort_values(["_row_id", "SNP"]).drop_duplicates(subset=["_row_id"], keep="first")

    out = pd.concat([mapped_snp, mapped_chr_bp], ignore_index=True)
    flipped = int(flipped_snp + flipped_chr_bp)
    dropped = int(dropped_snp + dropped_chr_bp)

    keep_cols = [c for c in ["SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE"] if c in out.columns]
    if keep_cols:
        out = out[keep_cols].sort_values(["CHR", "BP"]).drop_duplicates(subset=["SNP"], keep="first")
    else:
        out = pd.DataFrame(columns=["SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE"])

    if "BP" in out.columns:
        out["BP"] = pd.to_numeric(out["BP"], errors="coerce")
        out = out.dropna(subset=["BP"]).copy()
        out["BP"] = out["BP"].round().astype("int64")
    out.to_csv(args.out_tsv, sep="\t", index=False)

    n_mapped_snp = int(mapped_snp["_row_id"].nunique()) if not mapped_snp.empty else 0
    n_mapped_chr_bp = int(mapped_chr_bp["_row_id"].nunique()) if not mapped_chr_bp.empty else 0
    n_mapped_total = n_mapped_snp + n_mapped_chr_bp

    diag = {
        "ref_prefix": args.ref_prefix,
        "n_input_snps": int(gwas.shape[0]),
        "n_candidate_snpid_matches": int(merged_snp.shape[0]),
        "n_candidate_chr_bp_matches": int(merged_chr_bp.shape[0]),
        "n_mapped_by_snpid": n_mapped_snp,
        "n_mapped_by_chr_bp": n_mapped_chr_bp,
        "n_unmapped": int(gwas.shape[0] - n_mapped_total),
        "n_output_snps": int(out.shape[0]),
        "n_flipped": int(flipped),
        "n_dropped_allele_mismatch": int(dropped),
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
