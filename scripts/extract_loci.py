#!/usr/bin/env python3
import argparse
import json
import os
from typing import Dict, Tuple

import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    col_map = {
        "chrom": "CHR",
        "chr": "CHR",
        "bp": "BP",
        "pos": "BP",
        "position": "BP",
        "p": "P",
        "pval": "P",
        "pvalue": "P",
        "snp": "SNP",
        "rsid": "SNP",
        "snpid": "SNP",
        "markername": "SNP",
        "varid": "SNP",
        "variant_id": "SNP",
        "id": "SNP",
        "a1": "A1",
        "ea": "A1",
        "effect_allele": "A1",
        "a2": "A2",
        "nea": "A2",
        "other_allele": "A2",
        "beta": "BETA",
        "se": "SE",
    }
    rename = {}
    for c in df.columns:
        key = str(c).strip().lower()
        if key in col_map:
            rename[c] = col_map[key]
    out = df.rename(columns=rename)

    # Some inputs include multiple SNP-like columns (e.g., snpid and markername).
    # Choose the most informative source (highest uniqueness / variant-like tokens),
    # then fall back to other columns only when values are missing.
    if list(out.columns).count("SNP") > 1:
        snp_frame = out.loc[:, out.columns == "SNP"].copy()
        snp_frame.columns = [f"SNP_{i}" for i in range(snp_frame.shape[1])]
        snp_clean = snp_frame.astype(str).replace({"nan": "", "None": ""})

        def _score_series(s: pd.Series) -> tuple:
            non_empty = s[s != ""]
            n_non_empty = int(non_empty.shape[0])
            n_unique = int(non_empty.nunique())
            has_colon = int(non_empty.str.contains(":", regex=False).sum())
            # Prefer more unique IDs; colon-rich IDs (chr:pos:ref:alt) often carry full variant identity.
            return (n_unique, has_colon, n_non_empty)

        best_col = sorted(snp_clean.columns, key=lambda c: _score_series(snp_clean[c]), reverse=True)[0]
        ordered_cols = [best_col] + [c for c in snp_clean.columns if c != best_col]
        out["SNP"] = snp_clean[ordered_cols].bfill(axis=1).iloc[:, 0]
        out = out.loc[:, ~out.columns.duplicated()]

    required = ["CHR", "BP", "P", "SNP", "A1", "A2"]
    missing = [x for x in required if x not in out.columns]
    if missing:
        raise ValueError(f"GWAS file missing required columns after standardization: {missing}")

    out["CHR"] = out["CHR"].astype(str).str.replace("chr", "", regex=False)
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce")
    out = out.dropna(subset=["BP", "P"]).copy()
    out["BP"] = out["BP"].round().astype("int64")

    if "BETA" in out.columns:
        out["BETA"] = pd.to_numeric(out["BETA"], errors="coerce")
    if "SE" in out.columns:
        out["SE"] = pd.to_numeric(out["SE"], errors="coerce")

    return out


def _read_loci(loci_path: str) -> Dict[str, Tuple[str, int, int]]:
    loci = {}
    with open(loci_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() in {"locus", "loci", "id"}:
                continue
            if len(parts) >= 4:
                try:
                    start = int(parts[2])
                    end = int(parts[3])
                except ValueError:
                    # Skip malformed or header-like rows.
                    continue
                loci[parts[0]] = (str(parts[1]).replace("chr", ""), start, end)
            else:
                loci[parts[0]] = ("", -1, -1)
    return loci


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract one locus GWAS and standardize columns")
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_tsv)
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    gwas = pd.read_csv(args.gwas, sep=None, engine="python")
    gwas = _standardize_columns(gwas)

    loci = _read_loci(args.loci)
    if args.locus in loci:
        chr_, start, end = loci[args.locus]
        if start >= 0:
            gwas = gwas[(gwas["CHR"] == chr_) & (gwas["BP"] >= start) & (gwas["BP"] <= end)].copy()
    else:
        chr_, start, end = "", -1, -1

    gwas = gwas.sort_values(["CHR", "BP"]).drop_duplicates(subset=["SNP"], keep="first")
    gwas.to_csv(args.out_tsv, sep="\t", index=False)

    diag = {
        "locus": args.locus,
        "requested_window": {"chr": chr_, "start": start, "end": end},
        "n_snps_extracted": int(gwas.shape[0]),
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
