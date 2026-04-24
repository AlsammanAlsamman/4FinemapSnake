import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module, get_software_param

RESULTS_DIR = get_results_dir()

rule run_susier:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        bed=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bed",
        bim=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
        fam=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.fam",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_credible_set.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/diagnostics.json",
        diag_log=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/diagnostics.log",
        plot_overview=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_overview.png",
        plot_pip=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_pip.png",
        plot_ld=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_ld.png",
        plot_zscore=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_zscore.png",
        done=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier.done",
    params:
        L=lambda wc: int(get_analysis_value("susie_params.L")),
        coverage=lambda wc: float(get_analysis_value("susie_params.coverage")),
        window_kb=lambda wc: float(get_analysis_value("susie_params.window_kb")),
        sample_size=lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        plink_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/03_match/{wc.locus}/refpanel_matched",
        r_module=get_software_module("r"),
        r_libs_user=get_software_param("r", "r_libs_user", ""),
    log:
        f"{RESULTS_DIR}/log/07_susier_{{target}}_{{locus}}.log",
    resources:
        mem_mb=128000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load {params.r_module}
        export R_LIBS_USER="{params.r_libs_user}"
        Rscript scripts/run_susier.R \
          --input-tsv {input.matched} \
                    --plink-prefix {params.plink_prefix} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
                    --diag-log {output.diag_log} \
                    --plot-overview {output.plot_overview} \
                    --plot-pip {output.plot_pip} \
                    --plot-ld {output.plot_ld} \
                    --plot-zscore {output.plot_zscore} \
          --done-file {output.done} \
                    --sample-size {params.sample_size} \
          --L {params.L} \
          --coverage {params.coverage} \
                    --window-kb {params.window_kb} \
          > {log} 2>&1
        """
