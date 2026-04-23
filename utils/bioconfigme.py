from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import yaml

_BASE_DIR = Path(__file__).resolve().parents[1]
_ANALYSIS_PATH = _BASE_DIR / "configs" / "analysis.yml"
_SOFTWARE_PATH = _BASE_DIR / "configs" / "software.yml"

_analysis_cache: dict[str, Any] | None = None
_software_cache: dict[str, Any] | None = None


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing config file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config root must be a mapping: {path}")
    return data


def _get_analysis_config() -> dict[str, Any]:
    global _analysis_cache
    if _analysis_cache is None:
        _analysis_cache = _load_yaml(_ANALYSIS_PATH)
    return _analysis_cache


def _get_software_config() -> dict[str, Any]:
    global _software_cache
    if _software_cache is None:
        _software_cache = _load_yaml(_SOFTWARE_PATH)
    return _software_cache


def get_results_dir() -> str:
    """Return results directory from analysis config or inferred default."""
    cfg = _get_analysis_config()
    explicit = cfg.get("results_dir")
    if explicit:
        return str(explicit)

    project_name = str((cfg.get("metadata") or {}).get("project_name", "")).strip()
    if "analysis-" in project_name:
        suffix = project_name.split("analysis-", 1)[1]
        suffix = re.sub(r"[^a-zA-Z0-9]+", "_", suffix).strip("_").lower()
        if suffix:
            return f"results_{suffix}"

    cohort_name = str((cfg.get("cohort") or {}).get("name", "cohort")).strip().lower()
    cohort_name = re.sub(r"[^a-z0-9]+", "_", cohort_name).strip("_") or "cohort"
    return f"results_{cohort_name}"


def get_software_module(tool: str) -> str:
    """Return module name for a given tool in software.yml."""
    software = _get_software_config()
    tool_config = software.get(tool)
    if not isinstance(tool_config, dict):
        raise KeyError(f"software.yml missing tool: {tool}")
    module = tool_config.get("module")
    if not module:
        raise KeyError(f"software.yml missing module for tool: {tool}")
    return str(module)


def get_software_param(tool: str, param: str, default: Any = None) -> Any:
    """Return a parameter value from software.yml for a given tool.

    Example: get_software_param("r", "r_libs_user", "")
    """
    software = _get_software_config()
    tool_config = software.get(tool)
    if not isinstance(tool_config, dict):
        return default
    params = tool_config.get("params")
    if not isinstance(params, dict):
        return default
    return params.get(param, default)


def get_software_field(tool: str, field: str, default: Any = None) -> Any:
    """Return a top-level field from software.yml for a given tool.

    Example: get_software_field("gcta", "extension", "gcta")
    """
    software = _get_software_config()
    tool_config = software.get(tool)
    if not isinstance(tool_config, dict):
        return default
    return tool_config.get(field, default)


def get_gwas_validation_files() -> dict[str, str]:
    """Return {ancestry_label: file_path} from analysis.yml gwas_validation section."""
    cfg = _get_analysis_config()
    val_cfg = cfg.get("gwas_validation")
    if not isinstance(val_cfg, dict):
        return {}
    result: dict[str, str] = {}
    for anc, entry in val_cfg.items():
        if isinstance(entry, dict) and "gwas_file" in entry:
            result[str(anc)] = str(entry["gwas_file"])
        elif isinstance(entry, str):
            result[str(anc)] = entry
    return result


def get_analysis_value(path: str) -> Any:
    """Return nested value from analysis.yml using dot-separated path.

    Supports legacy keys used by existing rules by mapping them onto the new
    minimal analysis.yml structure.
    """
    cfg: Any = _get_analysis_config()

    # Direct lookup first.
    node: Any = cfg
    found = True
    for key in path.split("."):
        if not isinstance(node, dict) or key not in node:
            found = False
            break
        node = node[key]
    if found:
        return node

    results_dir = get_results_dir()
    cohort_cfg = cfg.get("cohort", {}) if isinstance(cfg.get("cohort"), dict) else {}
    ref_cfg = cfg.get("reference", {}) if isinstance(cfg.get("reference"), dict) else {}
    anc_cfg = cfg.get("ancestry", {}) if isinstance(cfg.get("ancestry"), dict) else {}
    data_cfg = cfg.get("data", {}) if isinstance(cfg.get("data"), dict) else {}
    resources_cfg = cfg.get("resources", {}) if isinstance(cfg.get("resources"), dict) else {}
    default_res_cfg = resources_cfg.get("default", {}) if isinstance(resources_cfg.get("default"), dict) else {}

    populations = [str(p).strip() for p in ref_cfg.get("populations", ["YRI", "PEL", "CEU", "CHB"])]
    populations = [p for p in populations if p]
    if not populations:
        raise ValueError("analysis.yml: reference.populations must contain at least one population")
    if len(set(populations)) != len(populations):
        raise ValueError("analysis.yml: reference.populations contains duplicates; order must be unique")

    pop_tag = "_".join(populations)
    pop_to_label = anc_cfg.get("population_to_label") or {}
    if not isinstance(pop_to_label, dict) or not pop_to_label:
        pop_to_label = {p: p for p in populations}
    pop_to_label = {str(k).strip(): str(v).strip() for k, v in pop_to_label.items()}

    missing_pop_labels = [p for p in populations if p not in pop_to_label or not pop_to_label[p]]
    if missing_pop_labels:
        raise ValueError(
            "analysis.yml: ancestry.population_to_label missing entries for populations: "
            + ", ".join(missing_pop_labels)
        )

    labels_ordered: list[str] = []
    for pop in populations:
        labels_ordered.append(pop_to_label.get(pop, pop))
    if len(set(labels_ordered)) != len(labels_ordered):
        raise ValueError(
            "analysis.yml: ancestry.population_to_label maps multiple populations to the same label; "
            "labels must be unique to keep ancestry codes deterministic"
        )

    label_to_code = {lbl: str(i) for i, lbl in enumerate(labels_ordered)}
    code_to_label = {str(i): lbl for i, lbl in enumerate(labels_ordered)}

    label_to_color = anc_cfg.get("label_to_color") or {}
    if not isinstance(label_to_color, dict):
        label_to_color = {}
    label_to_color = {str(k).strip(): str(v).strip() for k, v in label_to_color.items()}
    missing_label_colors = [lbl for lbl in labels_ordered if lbl not in label_to_color or not label_to_color[lbl]]
    if missing_label_colors:
        raise ValueError(
            "analysis.yml: ancestry.label_to_color missing colors for labels: "
            + ", ".join(missing_label_colors)
        )

    ref_ancestry = str(anc_cfg.get("reference_ancestry", labels_ordered[0])).strip()
    if ref_ancestry not in labels_ordered:
        raise ValueError(
            "analysis.yml: ancestry.reference_ancestry must match one ancestry label from population_to_label"
        )

    plink_prefix_name = str(cohort_cfg.get("plink_prefix", "Hispanic"))
    full_plink_prefix = plink_prefix_name if "/" in plink_prefix_name else f"inputs/{plink_prefix_name}"
    phenotype_format = str(cohort_cfg.get("phenotype_format", "fam"))
    case_code = str(cohort_cfg.get("case_code", 2))
    control_code = str(cohort_cfg.get("control_code", 1))

    loci_file = str(data_cfg.get("loci_file", "loci.txt"))
    loci_path = loci_file if "/" in loci_file else f"input/{loci_file}"
    prs_file = str(data_cfg.get("prs_file", "LAMR_PRS_using_MEX.txt"))
    prs_path = prs_file if "/" in prs_file else f"input/{prs_file}"
    cov_file = str(data_cfg.get("covariates_file", "covariates_Hispanic_5PC.txt"))
    cov_path = cov_file if "/" in cov_file else f"inputs/{cov_file}"
    cov_n_pcs = int(data_cfg.get("covariates_n_pcs", 5))

    prs_columns = data_cfg.get("prs_columns") if isinstance(data_cfg.get("prs_columns"), dict) else {}
    prs_iid = str(prs_columns.get("iid", "IID"))
    prs_score = str(prs_columns.get("score", "PRS"))
    prs_in_reg = str(prs_columns.get("in_regression", "In_Regression"))
    prs_in_reg_keep = str(prs_columns.get("in_regression_value", "Yes"))

    # Root-level minimal parameters.
    bin_size_bp = int(cfg.get("bin_size_bp", 250000))
    max_genes_per_bin = int(cfg.get("max_genes_per_bin", 3))
    ancestry_bins = int(cfg.get("ancestry_bins", 10))

    legacy_map: dict[str, Any] = {
        # first_step
        "first_step.panel_file": "resources/1000g_phase3/panel/integrated_call_samples_v3.20130502.ALL.panel",
        "first_step.output_dir": f"{results_dir}/reference_panel_prep",
        "first_step.populations": populations,
        "first_step.ancestry_labels": pop_to_label,

        # step2
        "step2.reference_release_dir": "resources/1000g_phase3/release_20130502",
        "step2.sample_list": f"{results_dir}/reference_panel_prep/{pop_tag}.samples.txt",
        "step2.output_dir": f"{results_dir}/reference_panel_subset",
        "step2.output_suffix": pop_tag,

        # step3
        "step3.plink_prefix": full_plink_prefix,
        "step3.loci_file": loci_path,
        "step3.output_dir": f"{results_dir}/cohort_region_subset",
        "step3.output_suffix": "loci",

        # step3c
        "step3c.output_dir": f"{results_dir}/cohort_region_subset_harmonized",
        "step3c.output_suffix": "loci.harmonized",

        # step4
        "step4.cohort_input_pattern": f"{results_dir}/cohort_region_subset_harmonized/chr{{chr}}.loci.harmonized.vcf.gz",
        "step4.reference_vcf_pattern": f"{results_dir}/reference_panel_subset/chr{{chr}}.{pop_tag}.vcf.gz",
        "step4.genetic_map_pattern": "resources/geneticmaps37/chr{chr}.b37.gmap.gz",
        "step4.output_dir": f"{results_dir}/cohort_phasing",

        # step5
        "step5.cohort_phased_pattern": f"{results_dir}/cohort_phasing/chr{{chr}}.phased.vcf.gz",
        "step5.reference_vcf_pattern": f"{results_dir}/reference_panel_subset/chr{{chr}}.{pop_tag}.vcf.gz",
        "step5.genetic_map_pattern": "resources/geneticmaps37/chr{chr}.b37.gmap.gz",
        "step5.sample_map_file": f"{results_dir}/reference_panel_prep/{pop_tag}.sample_map.txt",
        "step5.sample_map_label_codes": label_to_code,
        "step5.min_common_snps": int(cfg.get("step5", {}).get("min_common_snps", 2)) if isinstance(cfg.get("step5"), dict) else 2,
        "step5.output_dir": f"{results_dir}/rfmix",

        # step6
        "step6.rfmix_msp_pattern": f"{results_dir}/rfmix/chr{{chr}}.deconvoluted.msp.tsv",
        "step6.output_dir": f"{results_dir}/rfmix_plots",
        "step6.output_prefix_pattern": "chr{chr}.rfmix_ancestry_tracts",
        "step6.ancestry_labels": code_to_label,
        "step6.ancestry_colors": label_to_color,

        # step7
        "step7.rfmix_msp_pattern": f"{results_dir}/rfmix/chr{{chr}}.deconvoluted.msp.tsv",
        "step7.output_dir": f"{results_dir}/rfmix_bin_case_control",
        "step7.output_prefix_pattern": "chr{chr}.rfmix_case_control_bins",
        "step7.bin_size_bp": bin_size_bp,
        "step7.gene_annotation_file": "resources/NCBI37.3.gene.loc",
        "step7.max_genes_per_bin_label": max_genes_per_bin,
        "step7.phenotype_file": f"{full_plink_prefix}.fam",
        "step7.phenotype_format": phenotype_format,
        "step7.fam_case_code": case_code,
        "step7.fam_control_code": control_code,
        "step7.table_sample_column": "IID",
        "step7.table_status_column": "PHENO",
        "step7.table_case_value": "case",
        "step7.table_control_value": "control",
        "step7.ancestry_labels": code_to_label,
        "step7.ancestry_colors": label_to_color,

        # step8
        "step8.rfmix_q_pattern": f"{results_dir}/rfmix/chr{{chr}}.deconvoluted.rfmix.Q",
        "step8.ancestry_contribution_bins": ancestry_bins,
        "step8.phenotype_file": f"{full_plink_prefix}.fam",
        "step8.phenotype_format": phenotype_format,
        "step8.fam_case_code": case_code,
        "step8.fam_control_code": control_code,
        "step8.prs_file": prs_path,
        "step8.prs_iid_column": prs_iid,
        "step8.prs_score_column": prs_score,
        "step8.prs_in_regression_column": prs_in_reg,
        "step8.prs_in_regression_keep_value": prs_in_reg_keep,
        "step8.output_dir": f"{results_dir}/rfmix_prs_correlation",
        "step8.per_chr_prefix_pattern": "chr{chr}.rfmix_prs_correlation",
        "step8.matrix_prefix": "rfmix_prs_r2_matrix",
        "step8.ancestry_labels": code_to_label,
        "step8.ancestry_colors": label_to_color,

        # tractor
        "tractor.vcf_pattern": f"{results_dir}/cohort_phasing/chr{{chr}}.phased.vcf.gz",
        "tractor.msp_pattern": f"{results_dir}/rfmix/chr{{chr}}.deconvoluted.msp.tsv",
        "tractor.num_ancs": len(labels_ordered),
        "tractor.output_dir": f"{results_dir}/tractor",
        "tractor.output_vcf": False,
        "tractor.compress_output": False,

        # tractor_gwas
        "tractor_gwas.extract_output_dir": f"{results_dir}/tractor",
        "tractor_gwas.phenotype_file": f"{full_plink_prefix}.fam",
        "tractor_gwas.phenotype_format": phenotype_format,
        "tractor_gwas.fam_case_code": case_code,
        "tractor_gwas.fam_control_code": control_code,
        "tractor_gwas.table_iid_column": "IID",
        "tractor_gwas.table_phenotype_column": "PHENO",
        "tractor_gwas.output_dir": f"{results_dir}/tractor_gwas",
        "tractor_gwas.per_chr_prefix_pattern": "chr{chr}.model_{pair}",
        "tractor_gwas.merged_prefix_pattern": "merged.model_{pair}",
        "tractor_gwas.ancestry_labels": code_to_label,
        "tractor_gwas.min_partitions": 32,
        "tractor_gwas.covariates.enabled": True,
        "tractor_gwas.covariates.pca_file": cov_path,
        "tractor_gwas.covariates.format": "plink_eigenvec",
        "tractor_gwas.covariates.has_header": False,
        "tractor_gwas.covariates.fid_column": "FID",
        "tractor_gwas.covariates.iid_column": "IID",
        "tractor_gwas.covariates.n_pcs": cov_n_pcs,

        # analysis groups
        "annotation_meta.annotation": "resources/NCBI37.3.gene.loc",
        "maf_summary.output_dir": f"{results_dir}/maf_summary",
        "maf_gwas_summary.output_dir": f"{results_dir}/maf_gwas_summary",

        # locus ancestry report
        "locus_ancestry_report.input_tsv": f"{results_dir}/maf_gwas_summary/merged.maf_gwas_summary.tsv",
        "locus_ancestry_report.output_dir": f"{results_dir}/locus_ancestry_report",
        "locus_ancestry_report.cohort_maf_min": float(cfg.get("cohort_maf_min", 0.01)),
        "locus_ancestry_report.ancestry_maf_min": float(cfg.get("ancestry_maf_min", 0.01)),
        "locus_ancestry_report.local_ancestry_min_logp": float(cfg.get("local_ancestry_min_logp", 1.0)),
        "locus_ancestry_report.dosage_gwas_max_p": float(cfg.get("dosage_gwas_max_p", 5e-5)),
        "locus_ancestry_report.max_p_dosage": float(cfg.get("dosage_gwas_max_p", cfg.get("max_p_dosage", 0.01))),
        "locus_ancestry_report.min_or_dosage_deviation": float(cfg.get("min_or_dosage_deviation", 0.10)),
        "locus_ancestry_report.max_p_local_ancestry": float(cfg.get("max_p_local_ancestry", 10 ** (-float(cfg.get("local_ancestry_min_logp", 1.0))))),
        "locus_ancestry_report.min_or_local_deviation": float(cfg.get("min_or_local_deviation", 0.05)),
        "locus_ancestry_report.min_supporting_snps": int(cfg.get("min_supporting_snps", 1)),
        "locus_ancestry_report.clustering_method": str(cfg.get("clustering_method", "ward")),
        "locus_ancestry_report.clustering_metric": str(cfg.get("clustering_metric", "euclidean")),

        # admixture validation
        "admixture_validation.summary_tsv": f"{results_dir}/maf_gwas_summary/merged.maf_gwas_summary.tsv",
        "admixture_validation.rfmix_q_pattern": f"{results_dir}/rfmix/chr{{chr}}.deconvoluted.rfmix.Q",
        "admixture_validation.tractor_dosage_pattern": f"{results_dir}/tractor/chr{{chr}}.phased.anc{{anc}}.dosage.txt",
        "admixture_validation.output_dir": f"{results_dir}/admixture_validation",
        "admixture_validation.include_sex": _as_bool(cfg.get("include_sex", True), True),
        "admixture_validation.include_age": _as_bool(cfg.get("include_age", False), False),
        "admixture_validation.age_file": str(cfg.get("age_file", "")),
        "admixture_validation.age_iid_column": str(cfg.get("age_iid_column", "IID")),
        "admixture_validation.age_column": str(cfg.get("age_column", "AGE")),
        "admixture_validation.global_ancestry_reference": ref_ancestry,
        "admixture_validation.snp_pvalue_threshold": float(cfg.get("snp_pvalue_threshold", 0.01)),
        "admixture_validation.max_snps_per_ancestry": int(cfg.get("max_snps_per_ancestry", 50)),
        "admixture_validation.use_top_snp_per_locus": _as_bool(cfg.get("use_top_snp_per_locus", True), True),
        "admixture_validation.cv_folds": int(cfg.get("cv_folds", 5)),
        "admixture_validation.random_state": int(cfg.get("random_state", 42)),
        "admixture_validation.genome_build": str(cfg.get("genome_build", "GRCh37")),

        # admixture importance
        "admixture_importance.sample_features_tsv": f"{results_dir}/admixture_validation/sample_prediction_features.tsv",
        "admixture_importance.selected_snps_tsv": f"{results_dir}/admixture_validation/selected_admixture_snps.tsv",
        "admixture_importance.tractor_dosage_pattern": f"{results_dir}/tractor/chr{{chr}}.phased.anc{{anc}}.dosage.txt",
        "admixture_importance.output_dir": f"{results_dir}/admixture_importance",
        "admixture_importance.top_n_loci_plot": int(cfg.get("top_n_loci_plot", 20)),
        "admixture_importance.global_ancestry_reference": ref_ancestry,

        # resources compatibility
        "default_resources.mem_mb": int(default_res_cfg.get("mem_mb", 32000)),
        "default_resources.cores": int(default_res_cfg.get("cores", 2)),
        "default_resources.time": str(default_res_cfg.get("time", "00:30:00")),
    }

    if path in legacy_map:
        return legacy_map[path]

    raise KeyError(f"analysis.yml missing key path: {path}")


def get_default_resource(name: str, fallback: Any) -> Any:
    """Return default resources from either legacy or minimal config."""
    cfg = _get_analysis_config()
    legacy = cfg.get("default_resources", {})
    if isinstance(legacy, dict) and name in legacy:
        return legacy.get(name, fallback)
    modern = cfg.get("resources", {})
    if isinstance(modern, dict):
        modern_default = modern.get("default", {})
        if isinstance(modern_default, dict):
            return modern_default.get(name, fallback)
    return fallback


def _chromosome_sort_key(value: str) -> tuple[int, Any]:
    value = str(value).replace("chr", "")
    return (0, int(value)) if value.isdigit() else (1, value)


def get_loci_chromosomes(loci_file: str | None = None) -> list[str]:
    """Return sorted unique chromosomes from the configured loci file."""
    if not loci_file:
        loci_file = str(get_analysis_value("step3.loci_file"))

    path = Path(loci_file)
    if not path.is_absolute():
        path = _BASE_DIR / path

    if not path.exists():
        raise FileNotFoundError(f"Loci file not found: {path}")

    chromosomes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            if str(parts[0]).strip().lower() == "locus" and str(parts[1]).strip().lower() in {"chr", "chrom", "chromosome"}:
                continue
            chromosomes.add(str(parts[1]).replace("chr", ""))

    if not chromosomes:
        raise ValueError(f"No chromosomes found in loci file: {path}")

    return sorted(chromosomes, key=_chromosome_sort_key)


def get_step5_successful_chromosomes(fallback_to_loci: bool = True) -> list[str]:
    """Return chromosomes that produced usable RFMix outputs.

    When the Step 5 manifest does not exist yet, optionally fall back to the
    configured loci chromosomes so dry-runs and early-step parsing still work.
    """
    manifest = _BASE_DIR / get_results_dir() / "rfmix" / "successful_chromosomes.txt"

    if not manifest.exists():
        if fallback_to_loci:
            return get_loci_chromosomes()
        raise FileNotFoundError(f"RFMix successful chromosome manifest not found: {manifest}")

    chromosomes: list[str] = []
    with manifest.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            chromosomes.append(str(line).replace("chr", ""))

    chromosomes = sorted(set(chromosomes), key=_chromosome_sort_key)
    if not chromosomes:
        raise ValueError(f"RFMix successful chromosome manifest is empty: {manifest}")

    return chromosomes
