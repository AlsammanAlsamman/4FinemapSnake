#!/usr/bin/env python3
"""
Merge per-locus colocalization results into genome-wide summary tables.

This script aggregates colocalization results from all loci across all targets,
produces genome-wide summary tables, and optionally creates Excel workbooks.

Outputs:
  - Genome-wide colocalization results (PP.H4 sorted)
  - 95% credible set variants across all loci
  - Summary statistics by dataset/tissue
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import json
from datetime import datetime


def find_coloc_results(results_dir):
    """
    Recursively find all coloc_results.tsv files in the results directory.
    
    Expected structure: {results_dir}/{target}/{step}/coloc_results.tsv
    """
    coloc_files = []
    results_path = Path(results_dir)
    
    for tsv_file in results_path.glob("**/18_coloc/*/coloc_results.tsv"):
        coloc_files.append(tsv_file)
    
    return sorted(coloc_files)


def find_credible_set_results(results_dir):
    """Find all coloc_credible_set95.tsv files."""
    credible_files = []
    results_path = Path(results_dir)
    
    for tsv_file in results_path.glob("**/18_coloc/*/coloc_credible_set95.tsv"):
        credible_files.append(tsv_file)
    
    return sorted(credible_files)


def read_coloc_results(filepath):
    """Safely read coloc results with error handling."""
    try:
        if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
            return None
        
        df = pd.read_csv(filepath, sep="\t", dtype={"SNP": str})
        
        # Add locus identifier from filepath
        parts = filepath.parts
        if "18_coloc" in parts:
            idx = parts.index("18_coloc")
            if idx > 0:
                target = parts[idx - 1]
                if idx + 1 < len(parts):
                    locus = parts[idx + 1]
                    df["target"] = target
                    df["locus"] = locus
        
        return df
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
        return None


def merge_genome_wide_results(results_dir):
    """
    Merge all per-locus colocalization results.
    """
    print("[INFO] Finding colocalization result files...")
    coloc_files = find_coloc_results(results_dir)
    print(f"[INFO] Found {len(coloc_files)} result files")
    
    if len(coloc_files) == 0:
        print("[WARNING] No colocalization results found", file=sys.stderr)
        return None
    
    # Read all result files
    all_results = []
    n_skipped = 0
    
    for filepath in coloc_files:
        df = read_coloc_results(filepath)
        if df is not None and len(df) > 0:
            all_results.append(df)
        else:
            n_skipped += 1
    
    if len(all_results) == 0:
        print("[ERROR] No valid result files could be read", file=sys.stderr)
        return None
    
    print(f"[INFO] Successfully read {len(all_results)} files ({n_skipped} skipped)")
    
    # Concatenate all results
    genome_wide_df = pd.concat(all_results, ignore_index=True)
    
    # Sort by PP.H4 descending
    if "pp_h4" in genome_wide_df.columns:
        genome_wide_df = genome_wide_df.sort_values("pp_h4", ascending=False)
    
    print(f"[INFO] Merged {len(genome_wide_df)} total colocalization results")
    
    return genome_wide_df


def merge_credible_sets(results_dir):
    """
    Merge all credible set variant records.
    """
    print("[INFO] Finding credible set files...")
    credible_files = find_credible_set_results(results_dir)
    print(f"[INFO] Found {len(credible_files)} credible set files")
    
    if len(credible_files) == 0:
        print("[WARNING] No credible set files found", file=sys.stderr)
        return None
    
    # Read all credible set files
    all_credible = []
    n_skipped = 0
    
    for filepath in credible_files:
        df = read_coloc_results(filepath)
        if df is not None and len(df) > 0:
            all_credible.append(df)
        else:
            n_skipped += 1
    
    if len(all_credible) == 0:
        print("[WARNING] No valid credible set files could be read", file=sys.stderr)
        return None
    
    print(f"[INFO] Successfully read {len(all_credible)} credible set files ({n_skipped} skipped)")
    
    # Concatenate all credible sets
    credible_df = pd.concat(all_credible, ignore_index=True)
    
    # Sort by PP.H4 and SNP.PP.H4 descending
    if "pp_h4" in credible_df.columns and "SNP.PP.H4" in credible_df.columns:
        credible_df = credible_df.sort_values(
            ["pp_h4", "SNP.PP.H4"], 
            ascending=[False, False]
        )
    
    print(f"[INFO] Merged {len(credible_df)} total credible set variants")
    
    return credible_df


def compute_summary_statistics(genome_wide_df):
    """
    Compute summary statistics by dataset/tissue.
    """
    if genome_wide_df is None or len(genome_wide_df) == 0:
        return None
    
    summary_stats = []
    
    # Group by dataset
    if "dataset" in genome_wide_df.columns:
        grouped = genome_wide_df.groupby("dataset").agg({
            "pp_h4": ["count", "mean", "median", "min", "max"],
            "gene_id": "nunique"
        }).reset_index()
        
        # Flatten column names
        grouped.columns = ["_".join(col).strip("_") for col in grouped.columns.values]
        
        return grouped
    
    return None


def write_excel_summary(genome_wide_df, summary_stats, output_xlsx):
    """
    Write results to Excel workbook with multiple sheets.
    """
    try:
        import openpyxl
        from openpyxl.utils.dataframe import dataframe_to_rows
    except ImportError:
        print("[WARNING] openpyxl not available; skipping Excel output", file=sys.stderr)
        return
    
    try:
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            # Genome-wide results
            if genome_wide_df is not None and len(genome_wide_df) > 0:
                # Limit to top 10000 rows for Excel
                df_write = genome_wide_df.head(10000)
                df_write.to_excel(writer, sheet_name="Genome-Wide Results", index=False)
                print(f"[INFO] Wrote genome-wide results sheet ({len(df_write)} rows)")
            
            # Summary statistics
            if summary_stats is not None and len(summary_stats) > 0:
                summary_stats.to_excel(writer, sheet_name="Summary Stats", index=False)
                print(f"[INFO] Wrote summary statistics sheet")
        
        print(f"[INFO] Excel workbook written to {output_xlsx}")
    except Exception as e:
        print(f"[ERROR] Could not write Excel file: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-locus colocalization results into genome-wide summaries"
    )
    parser.add_argument("--results-dir", required=True, 
                       help="Directory containing per-locus results")
    parser.add_argument("--out-genome-wide-tsv", required=True,
                       help="Output file for genome-wide results")
    parser.add_argument("--out-credible-tsv", required=True,
                       help="Output file for credible set variants")
    parser.add_argument("--out-summary-xlsx", default=None,
                       help="Optional: Output Excel workbook with summary")
    parser.add_argument("--done-file", default=None,
                       help="Optional: Marker file indicating completion")
    
    args = parser.parse_args()
    
    print(f"[INFO] Starting coloc results merge at {datetime.now().isoformat()}")
    print(f"[INFO] Results directory: {args.results_dir}")
    
    # Merge genome-wide results
    genome_wide_df = merge_genome_wide_results(args.results_dir)
    
    if genome_wide_df is not None and len(genome_wide_df) > 0:
        # Write genome-wide results
        os.makedirs(os.path.dirname(args.out_genome_wide_tsv), exist_ok=True)
        genome_wide_df.to_csv(args.out_genome_wide_tsv, sep="\t", index=False)
        print(f"[INFO] Wrote genome-wide results to {args.out_genome_wide_tsv}")
    
    # Merge credible sets
    credible_df = merge_credible_sets(args.results_dir)
    
    if credible_df is not None and len(credible_df) > 0:
        # Write credible set results
        os.makedirs(os.path.dirname(args.out_credible_tsv), exist_ok=True)
        credible_df.to_csv(args.out_credible_tsv, sep="\t", index=False)
        print(f"[INFO] Wrote credible set variants to {args.out_credible_tsv}")
    
    # Compute summary statistics
    summary_stats = compute_summary_statistics(genome_wide_df)
    
    # Write Excel summary if requested
    if args.out_summary_xlsx:
        write_excel_summary(genome_wide_df, summary_stats, args.out_summary_xlsx)
    
    # Write done marker
    if args.done_file:
        os.makedirs(os.path.dirname(args.done_file), exist_ok=True)
        with open(args.done_file, 'w') as f:
            f.write(f"Merge completed at {datetime.now().isoformat()}\n")
        print(f"[INFO] Marked completion with {args.done_file}")
    
    print(f"[INFO] Merge complete at {datetime.now().isoformat()}")


if __name__ == "__main__":
    main()
