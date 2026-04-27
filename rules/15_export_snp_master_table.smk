import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

RESULTS_DIR = get_results_dir()


rule export_snp_master_table:
    input:
        matched=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/matched.tsv",
        afreq=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.afreq",
        bed=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bed",
        bim=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.bim",
        fam=f"{RESULTS_DIR}/{{target}}/03_match/{{locus}}/refpanel_matched.fam",
        finemap=f"{RESULTS_DIR}/{{target}}/07_finemap/{{locus}}/finemap_credible_set.tsv",
        susier=f"{RESULTS_DIR}/{{target}}/08_susier/{{locus}}/susier_credible_set.tsv",
        eqtl_combined=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/eqtl_subsets/all_eqtls.tsv",
        eqtl_manifest=f"{RESULTS_DIR}/{{target}}/14_coloc/{{locus}}/eqtl_subset_manifest.tsv",
    output:
        master=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/snp_master_table.tsv",
        ld_r2=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_r2.tsv",
        ld_r=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_r.tsv",
        ld_dprime=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/ld_matrix_Dprime.tsv",
        done=f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/snp_master_table.done",
    log:
        f"{RESULTS_DIR}/log/15_snp_master_{{target}}_{{locus}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/export_snp_master_table.py \
          --matched-tsv {input.matched} \
          --afreq {input.afreq} \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --finemap-tsv {input.finemap} \
          --susier-tsv {input.susier} \
          --eqtl-combined {input.eqtl_combined} \
          --eqtl-manifest {input.eqtl_manifest} \
          --out-master {output.master} \
          --out-ld-r2 {output.ld_r2} \
          --out-ld-r {output.ld_r} \
          --out-ld-dprime {output.ld_dprime} \
          --done-file {output.done} \
          > {log} 2>&1
        """


rule export_snp_master_table_all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/15_snp_master/{{locus}}/snp_master_table.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
    output:
        done=f"{RESULTS_DIR}/15_snp_master/all_snp_master_tables.done",
    resources:
        mem_mb=1000,
        time="00:10:00",
        cores=1,
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo ok > {output.done}
        """
