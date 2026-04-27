import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule causal_snp_analysis:
    input:
        master=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/snp_master_table.tsv",
        ld_r2=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_r2.tsv",
        ld_r=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_r.tsv",
        ld_dprime=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_Dprime.tsv",
    output:
        pdf=f"{RESULTS_DIR}/{{target}}/16_causal_snp_analysis/{{locus}}/all_plots.pdf",
        xlsx=f"{RESULTS_DIR}/{{target}}/16_causal_snp_analysis/{{locus}}/causal_snp_summary.xlsx",
        done=f"{RESULTS_DIR}/{{target}}/16_causal_snp_analysis/{{locus}}/causal_snp_analysis.done",
    log:
        f"{RESULTS_DIR}/log/16_causal_snp_analysis_{{target}}_{{locus}}.log",
    resources:
        mem_mb=64000,
        time="02:00:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/causal_snp_analysis.R \
          --snp-master {input.master} \
          --ld-r2 {input.ld_r2} \
          --ld-r {input.ld_r} \
          --ld-dprime {input.ld_dprime} \
          --out-pdf {output.pdf} \
          --out-xlsx {output.xlsx} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule causal_snp_analysis_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/16_causal_snp_analysis/{{locus}}/causal_snp_analysis.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/16_causal_snp_analysis/all_causal_snp_analysis.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
