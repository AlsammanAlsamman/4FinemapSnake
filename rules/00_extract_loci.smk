import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()

rule extract_locus:
    input:
        gwas=lambda wc: target_cfg(wc.target)["gwaspath"],
        loci=lambda wc: target_cfg(wc.target)["loci"],
        ref_bed=lambda wc: f"{refpanel_prefix_for_target(wc.target)}.bed",
        ref_bim=lambda wc: f"{refpanel_prefix_for_target(wc.target)}.bim",
        ref_fam=lambda wc: f"{refpanel_prefix_for_target(wc.target)}.fam",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/extracted.tsv",
        diag=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/diagnostics.json",
        ref_bed=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.bed",
        ref_bim=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.bim",
        ref_fam=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.fam",
        ref_afreq=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.afreq",
        ref_diag=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel_diagnostics.json",
        ref_done=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel_extract.done",
        done=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/extract.done",
    params:
        locus=lambda wc: wc.locus,
        ref_prefix=lambda wc: refpanel_prefix_for_target(wc.target),
        out_ref_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/01_extract/{wc.locus}/refpanel",
    log:
        f"{RESULTS_DIR}/log/00_extract_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/extract_loci.py \
          --gwas {input.gwas} \
          --loci {input.loci} \
          --locus {params.locus} \
          --out-tsv {output.tsv} \
          --diag-json {output.diag} \
          --done-file {output.done} \
                    > {log} 2>&1

                module load plink2/2.00a3.3lm
                python scripts/extract_refpanel_locus.py \
                    --ref-prefix {params.ref_prefix} \
                    --loci {input.loci} \
                    --locus {params.locus} \
                    --out-prefix {params.out_ref_prefix} \
                    --diag-json {output.ref_diag} \
                    --done-file {output.ref_done} \
                    >> {log} 2>&1
        """
