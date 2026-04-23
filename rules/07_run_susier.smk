import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module, get_software_param

RESULTS_DIR = get_results_dir()

rule run_susier:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld=f"{RESULTS_DIR}/{{target}}/06_ld_qc/{{locus}}/ld_matrix_qc.tsv",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_credible_set.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier.done",
    params:
        L=lambda wc: int(get_analysis_value("susie_params.L")),
        coverage=lambda wc: float(get_analysis_value("susie_params.coverage")),
        sample_size=lambda wc: int(get_analysis_value(f"targets.{wc.target}.samplesize")),
        r_module=get_software_module("r"),
        r_libs_user=get_software_param("r", "r_libs_user", ""),
    log:
        f"{RESULTS_DIR}/log/07_susier_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load {params.r_module}
        export R_LIBS_USER="{params.r_libs_user}"
        Rscript scripts/run_susier.R \
          --input-tsv {input.matched} \
          --ld-matrix {input.ld} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
                    --sample-size {params.sample_size} \
          --L {params.L} \
          --coverage {params.coverage} \
          > {log} 2>&1
        """
