#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import pandas as pd


RSID_CANDIDATES = ["rsid", "rsID", "SNP", "snp", "markername", "MarkerName"]
CHR_CANDIDATES = ["chr", "chrom", "chromosome", "CHR", "#CHR"]
POS_CANDIDATES = ["pos", "bp", "position", "POS", "BP"]


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _find_rsid_column(columns) -> str | None:
    lowered = {str(col).lower(): str(col) for col in columns}
    for candidate in RSID_CANDIDATES:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def _find_column(columns, candidates) -> str | None:
    lowered = {str(col).lower(): str(col) for col in columns}
    for candidate in candidates:
        if str(candidate).lower() in lowered:
            return lowered[str(candidate).lower()]
    return None


def _normalize_chr(value) -> str:
    text = str(value).strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        text = text[3:]
    return text.upper()


def _selected_snps(summary: pd.DataFrame, cojo_label: str) -> pd.DataFrame:
    finemap_hit = pd.to_numeric(summary.get("finemap_pip"), errors="coerce").notna()
    susie_hit = pd.to_numeric(summary.get("susie_pip"), errors="coerce").notna()
    cojo_hit = summary.get("cojo", pd.Series(index=summary.index, dtype="object")).fillna("").astype(str).eq(cojo_label)

    keep = finemap_hit | susie_hit | cojo_hit
    out = summary.loc[keep, ["SNP", "CHR", "BP", "EA", "OA", "AF", "P", "BETA", "OR", "SE", "cojo", "finemap_pip", "finemap_cs", "susie_pip", "susie_cs"]].copy()
    out["selected_by_cojo"] = cojo_hit.loc[keep].astype(int)
    out["selected_by_finemap"] = finemap_hit.loc[keep].astype(int)
    out["selected_by_susie"] = susie_hit.loc[keep].astype(int)
    return out.drop_duplicates(subset=["SNP"]).reset_index(drop=True)


def _extract_locus_bounds(summary: pd.DataFrame) -> tuple[str, int, int]:
    chr_col = _find_column(summary.columns, CHR_CANDIDATES)
    pos_col = _find_column(summary.columns, POS_CANDIDATES)
    if chr_col is None or pos_col is None:
        raise ValueError("summary.tsv must include chromosome and position columns (e.g., CHR/BP)")

    chr_series = summary[chr_col].map(_normalize_chr)
    chr_non_empty = chr_series[chr_series != ""]
    if chr_non_empty.empty:
        raise ValueError("summary.tsv chromosome column is empty after normalization")

    locus_chr = chr_non_empty.mode().iloc[0]
    pos_series = pd.to_numeric(summary[pos_col], errors="coerce")
    pos_series = pos_series[pos_series.notna()]
    if pos_series.empty:
        raise ValueError("summary.tsv position column has no valid numeric values")

    locus_start = int(pos_series.min())
    locus_end = int(pos_series.max())
    return locus_chr, locus_start, locus_end


def _empty_subset(header_cols) -> pd.DataFrame:
    return pd.DataFrame({col: pd.Series(dtype="object") for col in header_cols})


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset all eQTL resources to SNPs inside locus genomic range")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--eqtl-dir", required=True)
    parser.add_argument("--cojo-label", default="iter1")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--selected-snps-tsv", required=True)
    parser.add_argument("--manifest-tsv", required=True)
    parser.add_argument("--combined-tsv", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.selected_snps_tsv)
    _mkdir_for(args.manifest_tsv)
    _mkdir_for(args.combined_tsv)
    _mkdir_for(args.done_file)
    os.makedirs(args.out_dir, exist_ok=True)

    summary = pd.read_csv(args.summary_tsv, sep="\t")
    selected = _selected_snps(summary, args.cojo_label)

    locus_chr, locus_start, locus_end = _extract_locus_bounds(summary)
    selected["locus_chr"] = locus_chr
    selected["locus_start"] = locus_start
    selected["locus_end"] = locus_end
    selected.to_csv(args.selected_snps_tsv, sep="\t", index=False)

    eqtl_dir = Path(args.eqtl_dir)
    eqtl_files = sorted(eqtl_dir.glob("*.tsv"))

    manifest_rows = []
    combined_parts = []

    for eqtl_file in eqtl_files:
        header = pd.read_csv(eqtl_file, sep="\t", nrows=0)
        rsid_col = _find_rsid_column(header.columns)
        chr_col = _find_column(header.columns, CHR_CANDIDATES)
        pos_col = _find_column(header.columns, POS_CANDIDATES)
        out_path = Path(args.out_dir) / eqtl_file.name

        if chr_col is None or pos_col is None:
            _empty_subset(list(header.columns)).to_csv(out_path, sep="\t", index=False)
            manifest_rows.append({
                "dataset": eqtl_file.stem,
                "input_file": str(eqtl_file),
                "output_file": str(out_path),
                "rsid_column": rsid_col or "",
                "chr_column": chr_col or "",
                "pos_column": pos_col or "",
                "matched_rows": 0,
                "status": "missing_chr_or_pos_column",
            })
            continue

        matched_chunks = []
        for chunk in pd.read_csv(eqtl_file, sep="\t", chunksize=200000):
            chr_norm = chunk[chr_col].map(_normalize_chr)
            pos_num = pd.to_numeric(chunk[pos_col], errors="coerce")
            keep = (chr_norm == locus_chr) & pos_num.between(locus_start, locus_end)
            if keep.any():
                matched_chunks.append(chunk.loc[keep].copy())

        if matched_chunks:
            subset = pd.concat(matched_chunks, ignore_index=True, sort=False)
        else:
            subset = _empty_subset(list(header.columns))

        subset.to_csv(out_path, sep="\t", index=False)

        if not subset.empty:
            tagged = subset.copy()
            tagged.insert(0, "source_dataset", eqtl_file.stem)
            if rsid_col is not None:
                tagged.insert(1, "matched_rsid", tagged[rsid_col].astype(str))
            else:
                tagged.insert(1, "matched_rsid", pd.NA)
            combined_parts.append(tagged)

        manifest_rows.append({
            "dataset": eqtl_file.stem,
            "input_file": str(eqtl_file),
            "output_file": str(out_path),
            "rsid_column": rsid_col,
            "chr_column": chr_col,
            "pos_column": pos_col,
            "locus_chr": locus_chr,
            "locus_start": locus_start,
            "locus_end": locus_end,
            "matched_rows": int(len(subset)),
            "status": "ok",
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(args.manifest_tsv, sep="\t", index=False)

    if combined_parts:
        combined = pd.concat(combined_parts, ignore_index=True, sort=False)
    else:
        combined = pd.DataFrame(columns=["source_dataset", "matched_rsid"])
    combined.to_csv(args.combined_tsv, sep="\t", index=False)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()