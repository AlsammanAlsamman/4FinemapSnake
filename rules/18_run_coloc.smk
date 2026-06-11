import re
import sys
import os
from pathlib import Path
import glob as glob_module

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_software_module

# Build list of all results_* directories to support multiple datasets
RESULTS_DIRS = sorted([d for d in glob_module.glob("results_*") if os.path.isdir(d)])
if not RESULTS_DIRS:
    RESULTS_DIRS = ["results_Asian_IL12AB"]


def _safe_locus(value):
    """Sanitize locus names for use in paths."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def _load_loci_for_target(target, results_dir):
    """Load loci for a given target from config file."""
    loci_path = Path(config["targets"][target]["loci"])
    if not loci_path.exists():
        return ["locus1"]

    loci = []
    with loci_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() in {"locus", "loci", "id"}:
                continue
            loci.append(_safe_locus(parts[0]))

    return sorted(set(loci)) if loci else ["locus1"]


# Build target-locus pairs for each results directory
TARGET_LOCUS_PAIRS_PER_DIR = {}
TARGETS = sorted((config.get("targets") or {}).keys())
for results_dir in RESULTS_DIRS:
    TARGET_LOCUS_PAIRS = []
    for _target in TARGETS:
        for _locus in _load_loci_for_target(_target, results_dir):
            TARGET_LOCUS_PAIRS.append((_target, _locus))
    TARGET_LOCUS_PAIRS_PER_DIR[results_dir] = TARGET_LOCUS_PAIRS


# Retrieve all eQTL dataset names from config
def get_eqtl_datasets():
    """Extract available eQTL dataset names from analysis config."""
    datasets = {}
    try:
        eqtl_config = get_analysis_value("resources.eqtl_datasets")
        if eqtl_config and isinstance(eqtl_config, dict):
            datasets = {name: info for name, info in eqtl_config.items()}
    except Exception as e:
        print(f"Warning: Could not load eQTL datasets from config: {e}")
    return datasets


# Get list of dataset names for expansion
EQTL_DATASETS = list(get_eqtl_datasets().keys())


rule match_coloc_snps:
    """
    Match GWAS and eQTL subsets to shared SNPs in the same SNP order.
    """
    input:
        summary="{results_dir}/{target}/10_summary/{locus}/summary.tsv",
        eqtl_combined="{results_dir}/{target}/14_coloc/{locus}/eqtl_subsets/all_eqtls.tsv",
    output:
        matched_gwas="{results_dir}/{target}/18_coloc/{locus}/inputs/summary_matched.tsv",
        matched_eqtl="{results_dir}/{target}/18_coloc/{locus}/inputs/eqtl_matched.tsv",
        matched_snplist="{results_dir}/{target}/18_coloc/{locus}/inputs/shared_snps.snplist",
        diagnostics_json="{results_dir}/{target}/18_coloc/{locus}/inputs/match_diagnostics.json",
        done="{results_dir}/{target}/18_coloc/{locus}/inputs/match.done",
    log:
        "{results_dir}/log/18_coloc_match_{target}_{locus}.log",
    resources:
        mem_mb=16000,
        time="00:30:00",
        cores=1,
    shell:
        """
        python scripts/match_gwas_eqtl_snps.py \
          --gwas-tsv {input.summary} \
          --eqtl-tsv {input.eqtl_combined} \
          --out-gwas-tsv {output.matched_gwas} \
          --out-eqtl-tsv {output.matched_eqtl} \
          --out-snplist {output.matched_snplist} \
          --out-diagnostics-json {output.diagnostics_json} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule run_coloc:
    """
    Run colocalization analysis for a locus using extracted eQTLs.
    
    Performs Bayesian colocalization (coloc.abf) between GWAS signals at a locus
    and eQTL signals from multiple tissue/cell-type datasets. 
    
    Uses only SNPs selected by finemapping/fine-structure analysis (finemap, susie, cojo).
    """
    input:
        matched_summary="{results_dir}/{target}/18_coloc/{locus}/inputs/summary_matched.tsv",
        matched_eqtl="{results_dir}/{target}/18_coloc/{locus}/inputs/eqtl_matched.tsv",
        match_done="{results_dir}/{target}/18_coloc/{locus}/inputs/match.done",
        selected_snps="{results_dir}/{target}/14_coloc/{locus}/selected_snps.tsv",
    output:
        results_tsv="{results_dir}/{target}/18_coloc/{locus}/coloc_results.tsv",
        summary_tsv="{results_dir}/{target}/18_coloc/{locus}/coloc_summary_pp.tsv",
        credible_tsv="{results_dir}/{target}/18_coloc/{locus}/coloc_credible_set95.tsv",
        diagnostics_json="{results_dir}/{target}/18_coloc/{locus}/coloc_diagnostics.json",
        log_file="{results_dir}/{target}/18_coloc/{locus}/coloc_run.log",
        done="{results_dir}/{target}/18_coloc/{locus}/coloc.done",
    params:
        target=lambda wc: wc.target,
        locus=lambda wc: wc.locus,
        sample_size_gwas=lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        pp_h4_threshold=lambda wc: float(get_analysis_value("coloc_params.pp_h4_threshold")),
        pp_h4_credible=lambda wc: float(get_analysis_value("coloc_params.pp_h4_credible_set")),
        analysis_type=lambda wc: str(get_analysis_value("coloc_params.analysis_type")),
        output_all_hits=lambda wc: str(get_analysis_value("coloc_params.output_all_hits")).lower(),
        r_module=get_software_module("r"),
    log:
        "{results_dir}/log/18_coloc_{target}_{locus}.log",
    resources:
        mem_mb=64000,
        time="02:00:00",
        cores=4,
    shell:
        """
        module load {params.r_module}
        
        Rscript scripts/run_coloc.R \
                    --summary-tsv {input.matched_summary} \
          --selected-snps-tsv {input.selected_snps} \
                    --eqtl-combined-tsv {input.matched_eqtl} \
          --target {params.target} \
          --locus {params.locus} \
          --sample-size-gwas {params.sample_size_gwas} \
          --pp-h4-threshold {params.pp_h4_threshold} \
          --pp-h4-credible {params.pp_h4_credible} \
          --analysis-type {params.analysis_type} \
          --output-all-hits {params.output_all_hits} \
          --out-tsv {output.results_tsv} \
          --out-summary-tsv {output.summary_tsv} \
          --out-credible-tsv {output.credible_tsv} \
          --out-diagnostics-json {output.diagnostics_json} \
          --out-log-file {output.log_file} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule run_coloc_all:
    """
    Aggregate all colocalization analyses across all loci for a target.
    Uses standard expand with defined target-locus pairs.
    """
    input:
        lambda wc: expand(
            "{results_dir}/{target}/18_coloc/{locus}/coloc.done",
            results_dir=wc.results_dir,
            target=[t for t, _ in TARGET_LOCUS_PAIRS_PER_DIR.get(wc.results_dir, [])],
            locus=[l for _, l in TARGET_LOCUS_PAIRS_PER_DIR.get(wc.results_dir, [])],
        ),
    output:
        "{results_dir}/18_coloc/all_coloc.done",
    shell:
        """
        mkdir -p $(dirname {output})
        echo ok > {output}
        """


rule coloc_merge_results:
    """
    Merge colocalization results across all loci into genome-wide summary tables.
    """
    input:
        "{results_dir}/18_coloc/all_coloc.done",
    output:
        results_tsv="{results_dir}/18_coloc/coloc_genome_wide_results.tsv",
        credible_tsv="{results_dir}/18_coloc/coloc_genome_wide_credible_set95.tsv",
        summary_xlsx="{results_dir}/18_coloc/coloc_summary_statistics.xlsx",
    params:
        results_dir=lambda wc: wc.results_dir,
    shell:
        """
        /usr/local/analysis/python/3.10.2/bin/python3.10 scripts/merge_coloc_results.py \
          --results-dir {params.results_dir} \
          --out-tsv {output.results_tsv} \
          --out-credible-tsv {output.credible_tsv} \
          --out-xlsx {output.summary_xlsx}
        """
