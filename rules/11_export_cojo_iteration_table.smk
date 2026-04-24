import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule export_cojo_iteration_table:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        afreq=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.afreq",
        cojo_done=f"{RESULTS_DIR}/{{target}}/09_cojo_gcta/{{locus}}/cojo.done",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/11_cojo_iteration_table/{{locus}}/cojo_iteration_snps.tsv",
        xlsx=f"{RESULTS_DIR}/{{target}}/11_cojo_iteration_table/{{locus}}/cojo_iteration_snps.xlsx",
        diag=f"{RESULTS_DIR}/{{target}}/11_cojo_iteration_table/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/11_cojo_iteration_table/{{locus}}/cojo_iteration_table.done",
    params:
        cojo_tmp_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/09_cojo_gcta/{wc.locus}/gcta_tmp",
    log:
        f"{RESULTS_DIR}/log/11_cojo_iteration_table_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/export_cojo_iteration_table.py \
          --matched-tsv {input.matched} \
          --afreq {input.afreq} \
          --cojo-tmp-dir {params.cojo_tmp_dir} \
          --out-tsv {output.tsv} \
          --out-xlsx {output.xlsx} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule export_cojo_iteration_table_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/11_cojo_iteration_table/{{locus}}/cojo_iteration_table.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/11_cojo_iteration_table/all_cojo_iteration_tables.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
