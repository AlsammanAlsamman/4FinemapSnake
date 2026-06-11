import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()

rule build_ld_matrix:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ref_bim=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
    output:
        matrix=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.done",
    log:
        f"{RESULTS_DIR}/log/04_ld_matrix_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
                module load R/4.5.1-mkl
                Rscript scripts/build_ld_matrix.R \
          --input-tsv {input.matched} \
                    --ref-bim {input.ref_bim} \
          --out-matrix {output.matrix} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          > {log} 2>&1
        """
