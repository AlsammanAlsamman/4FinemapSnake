import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()

rule export_summary:
    input:
        finemap=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap_credible_set.tsv",
        susier=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_credible_set.tsv",
        cojo=f"{RESULTS_DIR}/{{target}}/09_cojo_gcta/{{locus}}/cojo_independent_signals.tsv",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.tsv",
        xlsx=f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.xlsx",
        done=f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.done",
    params:
        target=lambda wc: wc.target,
        locus=lambda wc: wc.locus,
    log:
        f"{RESULTS_DIR}/log/10_summary_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/export_summary.py \
          --target {params.target} \
          --locus {params.locus} \
          --finemap-tsv {input.finemap} \
          --susier-tsv {input.susier} \
          --cojo-tsv {input.cojo} \
          --out-tsv {output.tsv} \
          --out-xlsx {output.xlsx} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule export_summary_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/10_summary/all_summaries.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
