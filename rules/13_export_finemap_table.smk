import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule export_finemap_table:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        afreq=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.afreq",
        finemap=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap_credible_set.tsv",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/13_finemap_table/{{locus}}/finemap_all_snps.tsv",
        xlsx=f"{RESULTS_DIR}/{{target}}/13_finemap_table/{{locus}}/finemap_all_snps.xlsx",
        diag=f"{RESULTS_DIR}/{{target}}/13_finemap_table/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/13_finemap_table/{{locus}}/finemap_table.done",
    log:
        f"{RESULTS_DIR}/log/13_finemap_table_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/export_finemap_table.R \
          --matched-tsv {input.matched} \
          --afreq {input.afreq} \
          --finemap-tsv {input.finemap} \
          --out-tsv {output.tsv} \
          --out-xlsx {output.xlsx} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule export_finemap_table_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/13_finemap_table/{{locus}}/finemap_table.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/13_finemap_table/all_finemap_tables.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
