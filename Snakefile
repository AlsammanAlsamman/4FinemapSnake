import re
import sys
from pathlib import Path

configfile: "configs/analysis.yml"

sys.path.append("utils")
from bioconfigme import get_default_resource, get_results_dir

RESULTS_DIR = get_results_dir()
TARGETS = sorted((config.get("targets") or {}).keys())
if not TARGETS:
    raise ValueError("analysis.yml must define at least one target under 'targets'")


def _safe_locus(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def _load_loci_for_target(target):
    loci_path = Path(config["targets"][target]["loci"])
    if not loci_path.exists():
        return ["locus1"]

    loci = []
    with loci_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() in {"locus", "loci", "id"}:
                continue
            loci.append(_safe_locus(parts[0]))

    return sorted(set(loci)) if loci else ["locus1"]


TARGET_LOCUS_PAIRS = []
for _target in TARGETS:
    for _locus in _load_loci_for_target(_target):
        TARGET_LOCUS_PAIRS.append((_target, _locus))

PAIR_TARGETS = [x[0] for x in TARGET_LOCUS_PAIRS]
PAIR_LOCI = [x[1] for x in TARGET_LOCUS_PAIRS]

DEFAULT_MEM_MB = int(get_default_resource("mem_mb", 32000))
DEFAULT_CORES = int(get_default_resource("cores", 2))
DEFAULT_TIME = str(get_default_resource("time", "00:30:00"))

wildcard_constraints:
    target = r"[A-Za-z0-9_.-]+",
    locus = r"[A-Za-z0-9_.-]+"

include: "rules/common.smk"
include: "rules/00_extract_loci.smk"
include: "rules/01_harmonize.smk"
include: "rules/02_filter_match_snps.smk"
include: "rules/03_ldscore_diagnostics.smk"
include: "rules/04_make_ld_matrix.smk"
include: "rules/05_ld_matrix_qc_fix.smk"
include: "rules/06_run_finemap.smk"
include: "rules/07_run_susier.smk"
include: "rules/08_run_cojo_iterative.smk"
include: "rules/08b_run_cojo_iterative_gcta.smk"
include: "rules/09_plot_manhattan.smk"
include: "rules/10_export_summary.smk"

rule all:
    input:
        expand(
            f"{RESULTS_DIR}/{{target}}/10_summary/{{locus}}/summary.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        ),
        expand(
            f"{RESULTS_DIR}/{{target}}/04_ldscore_plot/{{locus}}/ldscore_plot.done",
            zip,
            target=PAIR_TARGETS,
            locus=PAIR_LOCI,
        )
