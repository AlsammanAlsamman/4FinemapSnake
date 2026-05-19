import sys
import os

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()


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


rule run_coloc:
    """
    Run colocalization analysis for a locus using extracted eQTLs.
    
    Performs Bayesian colocalization (coloc.abf) between GWAS signals at a locus
    and eQTL signals from multiple tissue/cell-type datasets. 
    
    Uses only SNPs selected by finemapping/fine-structure analysis (finemap, susie, cojo).
    """
    input:
        summary=f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.tsv",
        selected_snps=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/selected_snps.tsv",
        eqtl_combined=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/eqtl_subsets/all_eqtls.tsv",
    output:
        results_tsv=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc_results.tsv",
        summary_tsv=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc_summary_pp.tsv",
        credible_tsv=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc_credible_set95.tsv",
        diagnostics_json=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc_diagnostics.json",
        log_file=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc_run.log",
        done=f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc.done",
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
        f"{RESULTS_DIR}/log/18_coloc_{{target}}_{{locus}}.log",
    resources:
        mem_mb=64000,
        time="02:00:00",
        cores=4,
    shell:
        """
        module load {params.r_module}
        
        Rscript scripts/run_coloc.R \
          --summary-tsv {input.summary} \
          --selected-snps-tsv {input.selected_snps} \
          --eqtl-combined-tsv {input.eqtl_combined} \
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
    """
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/18_coloc/{{locus}}/coloc.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/18_coloc/all_coloc.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "All colocalization analyses completed at $(date)" > {output.done}
        """


rule coloc_merge_results:
    """
    Merge per-locus colocalization results into genome-wide summary tables.
    """
    input:
        all_done=f"{RESULTS_DIR}/18_coloc/all_coloc.done",
    output:
        genome_wide_tsv=f"{RESULTS_DIR}/18_coloc/coloc_genome_wide_results.tsv",
        genome_wide_credible_tsv=f"{RESULTS_DIR}/18_coloc/coloc_genome_wide_credible_set95.tsv",
        summary_xlsx=f"{RESULTS_DIR}/18_coloc/coloc_summary_statistics.xlsx",
        done=f"{RESULTS_DIR}/18_coloc/coloc_merge.done",
    resources:
        mem_mb=32000,
        time="01:00:00",
        cores=2,
    shell:
        """
        python scripts/merge_coloc_results.py \
          --results-dir {RESULTS_DIR} \
          --out-genome-wide-tsv {output.genome_wide_tsv} \
          --out-credible-tsv {output.genome_wide_credible_tsv} \
          --out-summary-xlsx {output.summary_xlsx} \
          --done-file {output.done}
        """
