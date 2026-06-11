import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule calculate_ldscores_all_snps:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld_matrix=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.tsv",
    output:
        table=f"{RESULTS_DIR}/{{target}}/04_ldscore_table/{{locus}}/ldscores_all_snps.tsv",
        done=f"{RESULTS_DIR}/{{target}}/04_ldscore_table/{{locus}}/ldscores.done",
    log:
        f"{RESULTS_DIR}/log/03b_ldscore_table_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/calculate_ldscores_table.py \
          --matched-tsv {input.matched} \
          --ld-matrix {input.ld_matrix} \
          --output-tsv {output.table} \
          --done-file {output.done} \
          > {log} 2>&1
        """