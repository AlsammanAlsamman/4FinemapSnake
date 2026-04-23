import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()

rule ldscore_diagnostics:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld_matrix=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.tsv",
    output:
        plot=f"{RESULTS_DIR}/{{target}}/04_ldscore_plot/{{locus}}/zscore_vs_ldscore.png",
        diag_tsv=f"{RESULTS_DIR}/{{target}}/04_ldscore_plot/{{locus}}/correlation.tsv",
        done=f"{RESULTS_DIR}/{{target}}/04_ldscore_plot/{{locus}}/ldscore_plot.done",
    params:
        require_corr=lambda wc: str(bool(get_analysis_value("diagnostics.require_ldscore_correlation"))).lower(),
        min_corr=lambda wc: float(get_analysis_value("diagnostics.min_ldscore_p_correlation")),
    log:
        f"{RESULTS_DIR}/log/03_ldscore_plot_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/plot_ldscore_vs_p.R \
          --input-tsv {input.matched} \
          --ld-matrix {input.ld_matrix} \
          --plot-png {output.plot} \
          --diag-tsv {output.diag_tsv} \
          --done-file {output.done} \
          --require-correlation {params.require_corr} \
          --min-correlation {params.min_corr} \
          > {log} 2>&1
        """
