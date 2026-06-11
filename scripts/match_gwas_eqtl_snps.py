#!/usr/bin/env python3
import argparse
import json
import os

import pandas as pd


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _find_eqtl_snp_col(df: pd.DataFrame) -> str:
    candidates = ["SNP", "rsID", "rsid", "snp", "matched_rsid"]
    lowered = {str(c).lower(): str(c) for c in df.columns}
    for col in candidates:
        key = col.lower()
        if key in lowered:
            return lowered[key]
    raise ValueError(
        "Could not find SNP column in eQTL table. Tried: "
        + ", ".join(candidates)
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Match GWAS and eQTL subsets to shared SNPs in identical order"
    )
    parser.add_argument("--gwas-tsv", required=True)
    parser.add_argument("--eqtl-tsv", required=True)
    parser.add_argument("--gwas-snp-col", default="SNP")
    parser.add_argument("--eqtl-snp-col", default=None)
    parser.add_argument("--out-gwas-tsv", required=True)
    parser.add_argument("--out-eqtl-tsv", required=True)
    parser.add_argument("--out-snplist", required=True)
    parser.add_argument("--out-diagnostics-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_gwas_tsv)
    _mkdir_for(args.out_eqtl_tsv)
    _mkdir_for(args.out_snplist)
    _mkdir_for(args.out_diagnostics_json)
    _mkdir_for(args.done_file)

    gwas = pd.read_csv(args.gwas_tsv, sep="\t", dtype="object")
    eqtl = pd.read_csv(args.eqtl_tsv, sep="\t", dtype="object")

    gwas_input_rows = int(gwas.shape[0])
    eqtl_input_rows = int(eqtl.shape[0])

    if args.gwas_snp_col not in gwas.columns:
        raise ValueError(f"GWAS SNP column not found: {args.gwas_snp_col}")

    eqtl_snp_col = args.eqtl_snp_col or _find_eqtl_snp_col(eqtl)
    if eqtl_snp_col not in eqtl.columns:
        raise ValueError(f"eQTL SNP column not found: {eqtl_snp_col}")

    gwas[args.gwas_snp_col] = gwas[args.gwas_snp_col].astype(str).str.strip()
    eqtl[eqtl_snp_col] = eqtl[eqtl_snp_col].astype(str).str.strip()

    gwas = gwas[gwas[args.gwas_snp_col] != ""].copy()
    eqtl[eqtl_snp_col] = eqtl[eqtl_snp_col].astype(str).str.strip()
    eqtl = eqtl[eqtl[eqtl_snp_col] != ""].copy()

    gwas["_gwas_order"] = range(gwas.shape[0])

    # Keep first occurrence for GWAS SNPs so we can enforce strict 1:1 order.
    gwas_unique = gwas.drop_duplicates(subset=[args.gwas_snp_col], keep="first").copy()

    shared = set(gwas_unique[args.gwas_snp_col]).intersection(set(eqtl[eqtl_snp_col]))

    gwas_out = gwas_unique[gwas_unique[args.gwas_snp_col].isin(shared)].copy()
    gwas_out = gwas_out.sort_values("_gwas_order")

    order_map = {
        snp: i for i, snp in enumerate(gwas_out[args.gwas_snp_col].tolist())
    }

    eqtl_out = eqtl[eqtl[eqtl_snp_col].isin(shared)].copy()
    eqtl_out["_gwas_order"] = eqtl_out[eqtl_snp_col].map(order_map)
    eqtl_out = eqtl_out.sort_values(["_gwas_order"], kind="mergesort")

    gwas_out = gwas_out.drop(columns=["_gwas_order"])
    eqtl_out = eqtl_out.drop(columns=["_gwas_order"])

    gwas_out.to_csv(args.out_gwas_tsv, sep="\t", index=False)
    eqtl_out.to_csv(args.out_eqtl_tsv, sep="\t", index=False)
    gwas_out[[args.gwas_snp_col]].to_csv(args.out_snplist, sep="\t", index=False, header=False)

    diagnostics = {
        "gwas_input_rows": gwas_input_rows,
        "eqtl_input_rows": eqtl_input_rows,
        "gwas_unique_snps": int(gwas_unique.shape[0]),
        "shared_snps": int(gwas_out.shape[0]),
        "gwas_output_rows": int(gwas_out.shape[0]),
        "eqtl_output_rows": int(eqtl_out.shape[0]),
        "gwas_snp_col": args.gwas_snp_col,
        "eqtl_snp_col": eqtl_snp_col,
    }

    with open(args.out_diagnostics_json, "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()