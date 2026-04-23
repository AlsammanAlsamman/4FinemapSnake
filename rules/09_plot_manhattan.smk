import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()

rule plot_manhattan:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        finemap=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap_credible_set.tsv",
        susier=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_credible_set.tsv",
        cojo=f"{RESULTS_DIR}/{{target}}/09_cojo/{{locus}}/cojo_independent_signals.tsv",
    output:
        png=f"{RESULTS_DIR}/{{target}}/10_manhattan/{{locus}}/manhattan.png",
        done=f"{RESULTS_DIR}/{{target}}/10_manhattan/{{locus}}/manhattan.done",
    log:
        f"{RESULTS_DIR}/log/09_manhattan_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/plot_manhattan.py \
          --matched-tsv {input.matched} \
          --finemap-tsv {input.finemap} \
          --susier-tsv {input.susier} \
          --cojo-tsv {input.cojo} \
          --plot-png {output.png} \
          --done-file {output.done} \
          > {log} 2>&1
        """
