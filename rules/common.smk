import os
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir


def target_cfg(target):
    return get_analysis_value(f"targets.{target}")


def refpanel_prefix_for_target(target):
    cfg = target_cfg(target)
    refpanel_name = cfg["refpanel"]
    return get_analysis_value(f"refpanels.{refpanel_name}.path")


def ensure_parent(path_value):
    os.makedirs(os.path.dirname(path_value), exist_ok=True)


RESULTS_DIR = get_results_dir()
