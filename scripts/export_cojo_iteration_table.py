#!/usr/bin/env python3
import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd


def _mkdir_for(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def _or_ci(beta_series: pd.Series, se_series: pd.Series):
    beta = pd.to_numeric(beta_series, errors="coerce")
    se = pd.to_numeric(se_series, errors="coerce")
    or_v = beta.apply(lambda b: math.exp(b) if pd.notna(b) else float("nan"))
    lo = (beta - 1.96 * se).apply(lambda b: math.exp(b) if pd.notna(b) else float("nan"))
    hi = (beta + 1.96 * se).apply(lambda b: math.exp(b) if pd.notna(b) else float("nan"))
    return or_v, lo, hi


def _load_afreq(path: Path) -> pd.DataFrame:
    af = pd.read_csv(path, sep=r"\s+", engine="python")
    needed = {"ID", "REF", "ALT", "ALT_FREQS"}
    missing = needed.difference(set(af.columns))
    if missing:
        raise ValueError(f"Missing columns in afreq file {path}: {sorted(missing)}")
    af = af[["ID", "REF", "ALT", "ALT_FREQS"]].copy()
    af["ID"] = af["ID"].astype(str)
    af["REF"] = af["REF"].astype(str).str.upper()
    af["ALT"] = af["ALT"].astype(str).str.upper()
    af["ALT_FREQS"] = pd.to_numeric(af["ALT_FREQS"], errors="coerce")
    return af


def _compute_effect_allele_freq(df: pd.DataFrame) -> pd.Series:
    def _pick(row):
        effect = str(row.get("refA", "")).upper()
        ref = str(row.get("REF", "")).upper()
        alt = str(row.get("ALT", "")).upper()
        altf = _safe_float(row.get("ALT_FREQS"))
        if not math.isfinite(altf):
            return float("nan")
        if effect == alt:
            return altf
        if effect == ref:
            return 1.0 - altf
        return float("nan")

    return df.apply(_pick, axis=1)


def _compute_other_allele(df: pd.DataFrame) -> pd.Series:
    def _pick(row):
        effect = str(row.get("refA", "")).upper()
        ref = str(row.get("REF", "")).upper()
        alt = str(row.get("ALT", "")).upper()
        if effect == ref and alt:
            return alt
        if effect == alt and ref:
            return ref
        return row.get("A2", "")

    return df.apply(_pick, axis=1)


def _fmt_scientific(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    return s.apply(lambda x: f"{x:.3e}" if pd.notna(x) else "")


def _fmt_fixed(series: pd.Series, digits: int) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    return s.apply(lambda x: f"{x:.{digits}f}" if pd.notna(x) else "")


def _load_iteration_snps(cond_path: Path) -> set:
    if not cond_path.exists():
        return set()
    snps = set()
    with cond_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            snps.add(line.split()[0])
    return snps


def main() -> None:
    parser = argparse.ArgumentParser(description="Export COJO iteration SNP table with pre/post conditioning OR and CI")
    parser.add_argument("--matched-tsv", required=True)
    parser.add_argument("--afreq", required=True)
    parser.add_argument("--cojo-tmp-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-xlsx", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    out_tsv = Path(args.out_tsv)
    out_xlsx = Path(args.out_xlsx)
    diag_json = Path(args.diag_json)
    done_file = Path(args.done_file)
    cojo_tmp_dir = Path(args.cojo_tmp_dir)

    _mkdir_for(out_tsv)
    _mkdir_for(out_xlsx)
    _mkdir_for(diag_json)
    _mkdir_for(done_file)

    matched = pd.read_csv(args.matched_tsv, sep="\t")
    required_matched = {"SNP", "A1", "A2", "BETA", "SE", "P", "CHR", "BP"}
    missing_matched = required_matched.difference(set(matched.columns))
    if missing_matched:
        raise ValueError(f"Missing matched columns: {sorted(missing_matched)}")

    matched = matched[["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]].copy()
    matched["SNP"] = matched["SNP"].astype(str)
    matched["A1"] = matched["A1"].astype(str).str.upper()
    matched["A2"] = matched["A2"].astype(str).str.upper()

    afreq = _load_afreq(Path(args.afreq))

    cma_files = sorted(cojo_tmp_dir.glob("iter_*.cma.cojo"))
    if not cma_files:
        raise FileNotFoundError(f"No COJO iteration files found in {cojo_tmp_dir}: expected iter_*.cma.cojo")

    all_rows_tsv = []
    all_rows_xlsx = []
    diag_iters = []

    for cma_path in cma_files:
        m = re.match(r"iter_(\d+)\.cma\.cojo$", cma_path.name)
        if not m:
            continue
        iter_n = int(m.group(1))
        iter_label = f"iter{iter_n}"
        cond_path = cojo_tmp_dir / f"iter_{iter_n:02d}.cond.snps"
        given_path = cojo_tmp_dir / f"iter_{iter_n:02d}.given.cojo"
        iter_conditioned_snps = _load_iteration_snps(cond_path)

        cma = pd.read_csv(cma_path, sep=r"\s+", engine="python")
        required_cma = {"SNP", "refA", "freq", "b", "se", "p", "bC", "bC_se", "pC"}
        missing_cma = required_cma.difference(set(cma.columns))
        if missing_cma:
            raise ValueError(f"Missing COJO columns in {cma_path}: {sorted(missing_cma)}")

        cma = cma[["SNP", "Chr", "bp", "refA", "freq", "b", "se", "p", "bC", "bC_se", "pC"]].copy()
        cma["SNP"] = cma["SNP"].astype(str)
        cma["refA"] = cma["refA"].astype(str).str.upper()

        if given_path.exists() and iter_conditioned_snps:
            given = pd.read_csv(given_path, sep=r"\s+", engine="python")
            needed_given = {"SNP", "Chr", "bp", "refA", "freq", "b", "se", "p"}
            missing_given = needed_given.difference(set(given.columns))
            if missing_given:
                raise ValueError(f"Missing columns in {given_path}: {sorted(missing_given)}")

            given = given[["SNP", "Chr", "bp", "refA", "freq", "b", "se", "p"]].copy()
            given["SNP"] = given["SNP"].astype(str)
            given["refA"] = given["refA"].astype(str).str.upper()
            given["bC"] = float("nan")
            given["bC_se"] = float("nan")
            given["pC"] = float("nan")

            missing_cond_snps = iter_conditioned_snps.difference(set(cma["SNP"]))
            if missing_cond_snps:
                add_rows = given[given["SNP"].isin(missing_cond_snps)].copy()
                cma = pd.concat([cma, add_rows], ignore_index=True)

        merged = cma.merge(matched, on="SNP", how="left", suffixes=("_cojo", "_matched"))
        merged = merged.merge(afreq, left_on="SNP", right_on="ID", how="left")

        merged["effect_allele_frequency_from_afreq"] = _compute_effect_allele_freq(merged)
        merged["effect_allele_frequency"] = merged["effect_allele_frequency_from_afreq"].where(
            merged["effect_allele_frequency_from_afreq"].notna(), pd.to_numeric(merged["freq"], errors="coerce")
        )
        merged["other_allele_effect_pair"] = _compute_other_allele(merged)

        merged["allele_alignment_ok"] = merged["A1"].eq(merged["refA"])

        or_before, lo_before, hi_before = _or_ci(merged["b"], merged["se"])
        or_after, lo_after, hi_after = _or_ci(merged["bC"], merged["bC_se"])

        out_num = pd.DataFrame(
            {
                "af": merged["effect_allele_frequency"],
                "snp": merged["SNP"],
                "chr": merged["Chr"].where(merged["Chr"].notna(), merged["CHR"]),
                "bp": merged["bp"].where(merged["bp"].notna(), merged["BP"]),
                "ea": merged["refA"],
                "oa": merged["other_allele_effect_pair"],
                "iter": merged["SNP"].apply(lambda s: iter_label if s in iter_conditioned_snps else ""),
                "p_pre_cond": pd.to_numeric(merged["p"], errors="coerce"),
                "or_pre_cond": or_before,
                "or_l95_pre_cond": lo_before,
                "or_u95_pre_cond": hi_before,
                "p_post_cond": pd.to_numeric(merged["pC"], errors="coerce"),
                "or_post_cond": or_after,
                "or_l95_post_cond": lo_after,
                "or_u95_post_cond": hi_after,
                "allele_align_ok": merged["allele_alignment_ok"],
                "cojo_freq": pd.to_numeric(merged["freq"], errors="coerce"),
                "af_ea_from_afreq": merged["effect_allele_frequency_from_afreq"],
            }
        )

        out_num["af"] = pd.to_numeric(out_num["af"], errors="coerce").round(3)
        out_num["or_pre_cond"] = pd.to_numeric(out_num["or_pre_cond"], errors="coerce").round(2)
        out_num["or_l95_pre_cond"] = pd.to_numeric(out_num["or_l95_pre_cond"], errors="coerce").round(2)
        out_num["or_u95_pre_cond"] = pd.to_numeric(out_num["or_u95_pre_cond"], errors="coerce").round(2)
        out_num["or_post_cond"] = pd.to_numeric(out_num["or_post_cond"], errors="coerce").round(2)
        out_num["or_l95_post_cond"] = pd.to_numeric(out_num["or_l95_post_cond"], errors="coerce").round(2)
        out_num["or_u95_post_cond"] = pd.to_numeric(out_num["or_u95_post_cond"], errors="coerce").round(2)

        out_num = out_num.sort_values(["p_pre_cond", "bp"], na_position="last")

        out_tsv_df = out_num.copy()
        out_tsv_df["af"] = _fmt_fixed(out_tsv_df["af"], 3)
        out_tsv_df["p_pre_cond"] = _fmt_scientific(out_tsv_df["p_pre_cond"])
        out_tsv_df["p_post_cond"] = _fmt_scientific(out_tsv_df["p_post_cond"])
        out_tsv_df["or_pre_cond"] = _fmt_fixed(out_tsv_df["or_pre_cond"], 2)
        out_tsv_df["or_l95_pre_cond"] = _fmt_fixed(out_tsv_df["or_l95_pre_cond"], 2)
        out_tsv_df["or_u95_pre_cond"] = _fmt_fixed(out_tsv_df["or_u95_pre_cond"], 2)
        out_tsv_df["or_post_cond"] = _fmt_fixed(out_tsv_df["or_post_cond"], 2)
        out_tsv_df["or_l95_post_cond"] = _fmt_fixed(out_tsv_df["or_l95_post_cond"], 2)
        out_tsv_df["or_u95_post_cond"] = _fmt_fixed(out_tsv_df["or_u95_post_cond"], 2)

        all_rows_tsv.append(out_tsv_df)
        all_rows_xlsx.append(out_num)

        diag_iters.append(
            {
                "iteration": iter_label,
                "cma_file": str(cma_path),
                "cond_snps_file": str(cond_path),
                "n_iteration_tagged_snps": int((out_num["iter"] == iter_label).sum()),
                "n_rows": int(out_num.shape[0]),
                "n_allele_alignment_false": int((~out_num["allele_align_ok"]).fillna(False).sum()),
            }
        )

    final_tsv = pd.concat(all_rows_tsv, ignore_index=True) if all_rows_tsv else pd.DataFrame()
    final_xlsx = pd.concat(all_rows_xlsx, ignore_index=True) if all_rows_xlsx else pd.DataFrame()
    final_tsv.to_csv(out_tsv, sep="\t", index=False)
    final_xlsx.to_excel(out_xlsx, index=False)

    diag = {
        "cojo_tmp_dir": str(cojo_tmp_dir),
        "n_iterations": len(diag_iters),
        "iterations": diag_iters,
        "n_output_rows": int(final_xlsx.shape[0]),
    }
    with open(diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
