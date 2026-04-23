import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()

rule qc_fix_ld_matrix:
    input:
        matrix=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.tsv",
    output:
        matrix=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_matrix_qc.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_qc.done",
    params:
        tolerance=lambda wc: float(get_analysis_value("filters.matrix_symmetry_tolerance")),
    log:
        f"{RESULTS_DIR}/log/05_ld_qc_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
                module load R/4.5.1-mkl
                Rscript scripts/qc_fix_ld_matrix.R \
          --input-matrix {input.matrix} \
          --out-matrix {output.matrix} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          --tolerance {params.tolerance} \
          > {log} 2>&1
        """
