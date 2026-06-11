import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir

RESULTS_DIR = get_results_dir()

rule filter_and_match_snps:
    input:
        harmonized=f"{RESULTS_DIR}/{{target}}/02_harmonize/{{locus}}/harmonized.tsv",
        ref_bed=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.bed",
        ref_bim=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.bim",
        ref_fam=f"{RESULTS_DIR}/{{target}}/01_extract/{{locus}}/refpanel.fam",
    output:
        tsv=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        snplist=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.snplist",
        ref_bed=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bed",
        ref_bim=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
        ref_fam=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.fam",
        ref_afreq=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.afreq",
        diag=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/diagnostics.json",
        done=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/match.done",
    params:
        ref_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/01_extract/{wc.locus}/refpanel",
        out_ref_prefix=lambda wc: f"{RESULTS_DIR}/{wc.target}/03_match/{wc.locus}/refpanel_matched",
        maf_min=lambda wc: float(get_analysis_value("filters.refpanel_maf_min")),
    log:
        f"{RESULTS_DIR}/log/02_match_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/filter_and_match_snps.py \
          --harmonized {input.harmonized} \
          --ref-prefix {params.ref_prefix} \
          --maf-min {params.maf_min} \
          --out-tsv {output.tsv} \
                    --out-snplist {output.snplist} \
          --diag-json {output.diag} \
          --done-file {output.done} \
          > {log} 2>&1

                module load plink2/2.00a3.3lm
                plink2 \
                    --bfile {params.ref_prefix} \
                    --extract {output.snplist} \
                    --make-bed \
                    --out {params.out_ref_prefix} \
                    >> {log} 2>&1

                plink2 \
                    --bfile {params.out_ref_prefix} \
                    --freq \
                    --out {params.out_ref_prefix} \
                    >> {log} 2>&1
        """
