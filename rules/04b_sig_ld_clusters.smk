import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()


rule cluster_significant_ld:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        ld_matrix=f"{RESULTS_DIR}/{{target}}/05_ld_matrix/{{locus}}/ld_matrix.tsv",
    output:
        summary=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/{{locus}}/cluster_summary.tsv",
        html=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/{{locus}}/clusters_manhattan_interactive.html",
        diag=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/{{locus}}/cluster.done",
    params:
        loci_file=lambda wc: str(config["targets"][wc.target]["loci"]),
        p_threshold=lambda wc: float(get_analysis_value("filters.gwas_pvalue_threshold")),
        r2_threshold=lambda wc: float(get_analysis_value("clustering_LD.sig_ld_cluster_r2_threshold")),
        distance_kb=lambda wc: float(get_analysis_value("clustering_LD.sig_distance_cluster_r2_threshold")),
        out_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/04b_sig_ld_clusters/{wc.locus}",
    log:
        f"{RESULTS_DIR}/log/04b_sig_ld_clusters_{{target}}_{{locus}}.log",
    resources:
        mem_mb=16000,
        time="00:30:00",
        cores=2,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/cluster_significant_ld.R \
          --matched-tsv {input.matched} \
          --ld-matrix {input.ld_matrix} \
          --loci-file {params.loci_file} \
          --locus {wildcards.locus} \
          --p-threshold {params.p_threshold} \
          --r2-threshold {params.r2_threshold} \
          --distance-kb {params.distance_kb} \
          --out-summary {output.summary} \
          --out-html {output.html} \
          --out-diag {output.diag} \
          --out-dir {params.out_dir} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule cluster_significant_ld_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/{{locus}}/cluster.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/04b_sig_ld_clusters/all.done",
    resources:
        mem_mb=1000,
        time="00:05:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
