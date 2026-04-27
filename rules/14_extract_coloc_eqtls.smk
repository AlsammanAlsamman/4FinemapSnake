import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()


rule extract_coloc_eqtls:
    input:
        summary=f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.tsv",
    output:
        selected_snps=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/selected_snps.tsv",
        manifest=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/eqtl_subset_manifest.tsv",
        combined=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/eqtl_subsets/all_eqtls.tsv",
        done=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/coloc_extract.done",
    params:
        eqtl_dir=lambda wc: str(get_analysis_value("resources.eqtls_dir")),
        out_dir=lambda wc: f"{RESULTS_DIR}/{wc.target}/14_coloc/{wc.locus}/eqtl_subsets",
        cojo_label=lambda wc: str(get_analysis_value("coloc_params.cojo_keep_label")),
    log:
        f"{RESULTS_DIR}/log/14_coloc_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="01:00:00",
        cores=2,
    shell:
        """
        python scripts/extract_coloc_eqtls.py \
          --summary-tsv {input.summary} \
          --eqtl-dir {params.eqtl_dir} \
          --cojo-label {params.cojo_label} \
          --out-dir {params.out_dir} \
          --selected-snps-tsv {output.selected_snps} \
          --manifest-tsv {output.manifest} \
          --combined-tsv {output.combined} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule extract_coloc_eqtls_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/coloc_extract.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/14_coloc/all_coloc_extract.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """