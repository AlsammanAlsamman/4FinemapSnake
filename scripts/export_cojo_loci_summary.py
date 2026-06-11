#!/usr/bin/env python3
"""Aggregate per-locus GCTA-COJO results into one genome-wide loci summary table.

Scans {results_dir}/{target}/09_cojo_gcta/{locus}/ for cojo_independent_signals.tsv
and diagnostics.json (written by run_cojo_iterative_gcta), and emits a single
TSV + Excel table with one row per (target, locus) describing the locus span,
the number of independent COJO signals found, the lead SNPs, and the key
analysis parameters/diagnostics.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def _mkdir_for(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _read_diagnostics(diag_path: Path) -> dict:
    if not diag_path.exists():
        return {}
    with diag_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_locus_span(matched_tsv: Path):
    if not matched_tsv.exists():
        return None, None, None, None
    matched = pd.read_csv(matched_tsv, sep="\t", usecols=["CHR", "BP"])
    chrom = matched["CHR"].iloc[0] if not matched.empty else None
    return chrom, int(matched["BP"].min()), int(matched["BP"].max()), int(matched.shape[0])


def _summarize_locus(signals_path: Path, results_dir: Path) -> dict:
    rel = signals_path.relative_to(results_dir)
    target, _step, locus = rel.parts[0], rel.parts[1], rel.parts[2]

    signals = pd.read_csv(signals_path, sep="\t")
    diag = _read_diagnostics(signals_path.parent / "diagnostics.json")
    chrom, locus_start, locus_end, n_snps_tested = _read_locus_span(
        results_dir / target / "03_match" / locus / "matched.tsv"
    )

    lead_snps = signals["lead_snp"].astype(str).tolist() if "lead_snp" in signals.columns else []
    top_row = signals.iloc[0] if not signals.empty else None

    return {
        "target": target,
        "locus": locus,
        "chr": chrom,
        "locus_start_bp": locus_start,
        "locus_end_bp": locus_end,
        "n_snps_tested": n_snps_tested,
        "n_independent_signals": int(signals.shape[0]),
        "lead_snps": ";".join(lead_snps),
        "top_lead_snp": top_row["lead_snp"] if top_row is not None else "",
        "top_p_original": top_row["p_original"] if top_row is not None else float("nan"),
        "min_p_conditional": pd.to_numeric(signals.get("p_conditional"), errors="coerce").min()
        if "p_conditional" in signals.columns
        else float("nan"),
        "sample_size": diag.get("sample_size"),
        "p_value_cutoff": diag.get("p_cutoff"),
        "max_iterations": diag.get("max_iterations"),
        "iterations_run": diag.get("iterations_run"),
        "stop_reason": diag.get("stop_reason"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a genome-wide COJO loci summary table (one row per locus)")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-xlsx", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_tsv = Path(args.out_tsv)
    out_xlsx = Path(args.out_xlsx)
    done_file = Path(args.done_file)

    _mkdir_for(out_tsv)
    _mkdir_for(out_xlsx)
    _mkdir_for(done_file)

    signal_files = sorted(results_dir.glob("*/09_cojo_gcta/*/cojo_independent_signals.tsv"))
    if not signal_files:
        raise FileNotFoundError(f"No cojo_independent_signals.tsv files found under {results_dir}")

    rows = []
    for signals_path in signal_files:
        try:
            rows.append(_summarize_locus(signals_path, results_dir))
        except Exception as exc:
            print(f"[WARN] Skipping {signals_path}: {exc}", file=sys.stderr)

    summary = pd.DataFrame(rows).sort_values(["target", "locus"]).reset_index(drop=True)
    summary.to_csv(out_tsv, sep="\t", index=False)
    summary.to_excel(out_xlsx, index=False)

    with open(done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
