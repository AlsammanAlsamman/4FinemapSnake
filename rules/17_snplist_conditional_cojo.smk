import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_field, get_software_module

RESULTS_DIR = get_results_dir()


rule snplist_conditional_cojo:
    """
    For each target/locus, condition the locus summary stats on every SNP from
    the user-supplied SNP list (evidances.snplist) that falls within the locus.
    Uses GCTA --cojo-cond for conditional analysis.

    Outputs per target/locus:
      plots/locus_unconditional.png          — Manhattan without conditioning
      plots/locus_conditional_<SNP>.png      — Manhattan conditioned on each SNP
      plots/impact_plot.png                  — % AUC decrease bar chart
      snplist_conditional_results.xlsx       — SNP table + impact summary (two sheets)
      conditional.done                       — marker file
    """
    input:
        matched  = f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ref_bed  = f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bed",
        ref_bim  = f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
        ref_fam  = f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.fam",
        snplist  = get_analysis_value("evidances.snplist"),
    output:
        done = f"{RESULTS_DIR}/{{target}}/17_snplist_conditional/{{locus}}/conditional.done",
        xlsx = f"{RESULTS_DIR}/{{target}}/17_snplist_conditional/{{locus}}/snplist_conditional_results.xlsx",
    params:
        sample_size  = lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        bfile_prefix = lambda wc: f"{RESULTS_DIR}/{wc.target}/03_match/{wc.locus}/refpanel_matched",
        out_dir      = lambda wc: f"{RESULTS_DIR}/{wc.target}/17_snplist_conditional/{wc.locus}",
        loci_file    = lambda wc: str(config["targets"][wc.target]["loci"]),
        gcta_module  = lambda wc: str(get_software_module("gcta")),
        gcta_bin     = lambda wc: str(get_software_field("gcta", "extension", "gcta")),
        r_libs_user  = "~/R/packages:/usr/local/analysis/bioconductor/3.19/R/library:/usr/local/analysis/R/4.4.1-mkl/lib/R/library",
    log:
        f"{RESULTS_DIR}/log/17_snplist_cojo_{{target}}_{{locus}}.log",
    resources:
        mem_mb = 32000,
        time   = "01:00:00",
        cores  = 2,
    shell:
        """
        module load {params.gcta_module}
        module load R/4.5.1-mkl
        export R_LIBS_USER="{params.r_libs_user}"
        mkdir -p {params.out_dir}
        Rscript scripts/run_snplist_cojo_conditional.R \
          --input-tsv   {input.matched} \
          --bfile-prefix {params.bfile_prefix} \
          --snp-list    {input.snplist} \
          --loci-file   {params.loci_file} \
          --locus       {wildcards.locus} \
          --out-dir     {params.out_dir} \
          --sample-size {params.sample_size} \
          --gcta-bin    {params.gcta_bin} \
          --done-file   {output.done} \
          > {log} 2>&1
        """


rule snplist_conditional_cojo_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/17_snplist_conditional/{{locus}}/conditional.done",
            zip,
            target = PAIR_TARGETS,
            locus  = PAIR_LOCI,
        ),
    output:
        done = f"{RESULTS_DIR}/17_snplist_conditional/all.done",
    resources:
        mem_mb = 1000,
        time   = "00:05:00",
        cores  = 1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
