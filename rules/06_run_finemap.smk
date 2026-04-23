import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()

rule run_finemap:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_matrix_qc.tsv",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap_credible_set.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap.done",
    params:
        max_causal=lambda wc: int(get_analysis_value("finemap_params.max_causal")),
    log:
        f"{RESULTS_DIR}/log/06_finemap_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        bash scripts/run_finemap.sh \
          --input-tsv {input.matched} \
          --ld-matrix {input.ld} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          --max-causal {params.max_causal} \
          > {log} 2>&1
        """
