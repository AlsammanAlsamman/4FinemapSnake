#!/usr/bin/env python3
import argparse
import json
import os

import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _load_ref_info(prefix: str) -> pd.DataFrame:
    ref = pd.read_csv(
        f"{prefix}.bim",
        sep="\t",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "REF_A1", "REF_A2"],
    )
    ref["REF_A1"] = ref["REF_A1"].astype(str).str.upper()
    ref["REF_A2"] = ref["REF_A2"].astype(str).str.upper()
    ref["ref_order"] = range(ref.shape[0])
    return ref[["SNP", "ref_order", "REF_A1", "REF_A2"]]


def _load_ref_maf(prefix: str) -> pd.DataFrame:
    afreq = f"{prefix}.afreq"
    if not os.path.exists(afreq):
        return pd.DataFrame(columns=["SNP", "REF_MAF"])

    freq = pd.read_csv(afreq, sep="\s+")
    if "ID" in freq.columns and "ALT_FREQS" in freq.columns:
        freq["ALT_FREQS"] = pd.to_numeric(freq["ALT_FREQS"], errors="coerce")
        freq["REF_MAF"] = freq["ALT_FREQS"].apply(lambda x: min(x, 1 - x) if pd.notna(x) else None)
        return freq[["ID", "REF_MAF"]].rename(columns={"ID": "SNP"})
    return pd.DataFrame(columns=["SNP", "REF_MAF"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter low-MAF and enforce GWAS/reference SNP order match")
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--ref-prefix", required=True)
    parser.add_argument("--maf-min", type=float, required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-snplist", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.out_snplist)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    gwas = pd.read_csv(args.harmonized, sep="\t")
    gwas["A1"] = gwas["A1"].astype(str).str.upper()
    gwas["A2"] = gwas["A2"].astype(str).str.upper()

    ref_info = _load_ref_info(args.ref_prefix)
    ref_maf = _load_ref_maf(args.ref_prefix)

    merged = gwas.merge(ref_info, on="SNP", how="inner")
    n_missing_ref = int(gwas.shape[0] - merged.shape[0])

    before_allele = merged.shape[0]
    merged = merged[(merged["A1"] == merged["REF_A1"]) & (merged["A2"] == merged["REF_A2"])].copy()
    n_allele_mismatch_removed = int(before_allele - merged.shape[0])

    if not ref_maf.empty:
        merged = merged.merge(ref_maf, on="SNP", how="left")
        before_maf = merged.shape[0]
        merged = merged[(merged["REF_MAF"].isna()) | (merged["REF_MAF"] >= args.maf_min)].copy()
        n_low_maf_removed = int(before_maf - merged.shape[0])
    else:
        n_low_maf_removed = 0

    merged = merged.sort_values("ref_order")

    # Persist SNP list in exact reference order so PLINK subset has the same SNP set/order.
    merged[["SNP"]].to_csv(args.out_snplist, sep="\t", index=False, header=False)

    merged = merged.drop(columns=["ref_order", "REF_A1", "REF_A2"])
    merged.to_csv(args.out_tsv, sep="\t", index=False)

    diag = {
        "ref_prefix": args.ref_prefix,
        "maf_min": args.maf_min,
        "n_input": int(gwas.shape[0]),
        "n_missing_in_reference": n_missing_ref,
        "n_removed_allele_mismatch": n_allele_mismatch_removed,
        "n_removed_low_maf": n_low_maf_removed,
        "n_output": int(merged.shape[0]),
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
