import sys
from pathlib import Path

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()


def _load_loci_for_target(target):
    """Helper to load locus IDs for a target."""
    loci_path = Path(config["targets"][target]["loci"])
    if not loci_path.exists():
        return []
    loci = []
    with loci_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() in {"locus", "loci", "id"}:
                continue
            loci.append(parts[0])
    return sorted(set(loci)) if loci else []


# Pre-compute split_loci_by_clusters inputs per target
SPLIT_INPUTS = {}
for target in TARGETS:
    loci_list = _load_loci_for_target(target)
    SPLIT_INPUTS[target] = [
        f"{RESULTS_DIR}/{target}/04b_sig_ld_clusters/{loc}/cluster.done"
        for loc in loci_list
    ]


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


rule split_loci_by_clusters:
    input:
        cluster_done=lambda wc: SPLIT_INPUTS.get(wc.target, []),
    output:
        loci=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/loci_split_by_clusters.tsv",
        done=f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/split.done",
    params:
        loci_file=lambda wc: str(config["targets"][wc.target]["loci"]),
        cluster_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/04b_sig_ld_clusters",
    log:
        f"{RESULTS_DIR}/log/split_loci_by_clusters_{{target}}.log",
    resources:
        mem_mb=4000,
        time="00:15:00",
        cores=1,
    shell:
        """
        module load R/4.5.1-mkl
        Rscript scripts/split_loci_by_clusters.R \
          --loci-file {params.loci_file} \
          --cluster-dir {params.cluster_dir} \
          --target {wildcards.target} \
          --out-loci {output.loci} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule split_loci_by_clusters_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/04b_sig_ld_clusters/split.done",
            target=TARGETS,
        ),
    output:
        done=f"{RESULTS_DIR}/04b_sig_ld_clusters/split_all.done",
    resources:
        mem_mb=1000,
        time="00:05:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
