import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()

rule run_cojo_iterative:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_matrix_qc.tsv",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/09_cojo/{{locus}}/cojo_independent_signals.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/09_cojo/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/09_cojo/{{locus}}/cojo.done",
    params:
        p_cutoff=lambda wc: float(get_analysis_value("cojo_params.pvalue_cutoff")),
        max_iter=lambda wc: int(get_analysis_value("cojo_params.max_iterations")),
        sample_size=lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        plot_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/09_cojo/{wc.locus}/iterations",
    log:
        f"{RESULTS_DIR}/log/08_cojo_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
                module load R/4.5.1-mkl
                Rscript scripts/run_cojo_iterative.R \
          --input-tsv {input.matched} \
          --ld-matrix {input.ld} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
                --sample-size {params.sample_size} \
          --p-cutoff {params.p_cutoff} \
          --max-iterations {params.max_iter} \
                    --plot-dir {params.plot_dir} \
          > {log} 2>&1
        """
