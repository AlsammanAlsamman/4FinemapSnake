#!/usr/bin/env bash
set -euo pipefail

TARGET=${1:-}
if [[ -z "$TARGET" ]]; then
  echo "Usage: $0 <done-target-path>" >&2
  exit 1
fi

snakemake \
  --snakefile rules/example_complete.smk \
  --configfile configs/analysis.yml \
  --cores 2 \
  "$TARGET"
