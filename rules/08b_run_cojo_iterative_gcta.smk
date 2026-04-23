import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_field, get_software_module

RESULTS_DIR = get_results_dir()

rule run_cojo_iterative_gcta:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ref_bed=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bed",
        ref_bim=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
        ref_fam=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.fam",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/09_cojo_gcta/{{locus}}/cojo_independent_signals.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/09_cojo_gcta/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/09_cojo_gcta/{{locus}}/cojo.done",
    params:
        p_cutoff=lambda wc: float(get_analysis_value("cojo_params.pvalue_cutoff")),
        max_iter=lambda wc: int(get_analysis_value("cojo_params.max_iterations")),
        sample_size=lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        bfile_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/03_match/{wc.locus}/refpanel_matched",
        plot_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/09_cojo_gcta/{wc.locus}/iterations",
        gcta_module=lambda wc: str(get_software_module("gcta")),
        gcta_bin=lambda wc: str(get_software_field("gcta", "extension", "gcta")),
    log:
        f"{RESULTS_DIR}/log/08b_cojo_gcta_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load {params.gcta_module}
        module load R/4.5.1-mkl
        Rscript scripts/run_cojo_iterative_gcta.R \
          --input-tsv {input.matched} \
          --bfile-prefix {params.bfile_prefix} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          --sample-size {params.sample_size} \
          --p-cutoff {params.p_cutoff} \
          --max-iterations {params.max_iter} \
          --plot-dir {params.plot_dir} \
          --gcta-bin {params.gcta_bin} \
          > {log} 2>&1
        """
