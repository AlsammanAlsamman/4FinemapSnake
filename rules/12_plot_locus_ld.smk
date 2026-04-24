import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule plot_locus_ld:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld_matrix=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_matrix_qc.tsv",
    output:
        pdf=f"{RESULTS_DIR}/{{target}}/12_locus_ld_plot/{{locus}}/locus_ld.pdf",
        png=f"{RESULTS_DIR}/{{target}}/12_locus_ld_plot/{{locus}}/locus_ld.png",
        done=f"{RESULTS_DIR}/{{target}}/12_locus_ld_plot/{{locus}}/locus_ld.done",
    params:
        window_bp=config.get("locus_ld_window_bp", 100000),
    log:
        f"{RESULTS_DIR}/log/12_locus_ld_plot_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/plot_locus_ld.R \
          --matched-tsv {input.matched} \
          --ld-matrix   {input.ld_matrix} \
          --out-pdf     {output.pdf} \
          --out-png     {output.png} \
          --window-bp   {params.window_bp} \
          --done-file   {output.done} \
          > {log} 2>&1
        """


rule plot_locus_ld_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/12_locus_ld_plot/{{locus}}/locus_ld.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/12_locus_ld_plot/all_locus_ld_plots.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
