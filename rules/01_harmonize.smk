import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()

rule harmonize_with_refpanel:
    input:
        gwas=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/extracted.tsv",
        ref_bim=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.bim",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/02_harmonize/{{locus}}/harmonized.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/02_harmonize/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/02_harmonize/{{locus}}/harmonize.done",
    params:
        ref_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/01_extract/{wc.locus}/refpanel",
    log:
        f"{RESULTS_DIR}/log/01_harmonize_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/harmonize_gwas.py \
          --gwas {input.gwas} \
          --ref-prefix {params.ref_prefix} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          > {log} 2>&1
        """
