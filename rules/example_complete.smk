import sys

sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module

RESULTS_DIR = get_results_dir()
PYTHON_MODULE = get_software_module("python")

rule example_complete:
    input:
        "inputs/snps.txt",
    output:
        done=f"{RESULTS_DIR}/example/sample_{{sample}}/example.done",
    params:
        module=PYTHON_MODULE,
        sample=lambda wc: wc.sample,
    log:
        f"{RESULTS_DIR}/log/example_complete_{{sample}}.log",
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2,
    shell:
        """
        python scripts/example_helper.py \
          --sample {params.sample} \
          --module {params.module} \
          --input-file {input} \
          --done-file {output.done} \
          > {log} 2>&1
        """
