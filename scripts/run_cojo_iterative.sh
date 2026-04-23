#!/usr/bin/env bash
set -euo pipefail

INPUT_TSV=""
LD_MATRIX=""
OUT_TSV=""
DIAG_JSON=""
DONE_FILE=""
P_CUTOFF="5e-8"
MAX_ITER="25"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-tsv) INPUT_TSV="$2"; shift 2 ;;
    --ld-matrix) LD_MATRIX="$2"; shift 2 ;;
    --out-tsv) OUT_TSV="$2"; shift 2 ;;
    --diag-json) DIAG_JSON="$2"; shift 2 ;;
    --done-file) DONE_FILE="$2"; shift 2 ;;
    --p-cutoff) P_CUTOFF="$2"; shift 2 ;;
    --max-iterations) MAX_ITER="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$(dirname "$OUT_TSV")" "$(dirname "$DIAG_JSON")" "$(dirname "$DONE_FILE")"

python - "$INPUT_TSV" "$OUT_TSV" "$DIAG_JSON" "$P_CUTOFF" "$MAX_ITER" <<'PY'
import json
import sys
import pandas as pd

input_tsv, out_tsv, diag_json, p_cutoff, max_iter = sys.argv[1:6]
p_cutoff = float(p_cutoff)
max_iter = int(max_iter)

df = pd.read_csv(input_tsv, sep="\t")
if "SNP" not in df.columns or "P" not in df.columns:
    raise ValueError("Matched table must include SNP and P")

work = df.copy()
work["P"] = pd.to_numeric(work["P"], errors="coerce")
work = work.dropna(subset=["P"]).sort_values("P")

signals = []
for i in range(max_iter):
    sig = work[work["P"] <= p_cutoff]
    if sig.empty:
        break
    top = sig.iloc[0].copy()
    top["iteration"] = i + 1
    signals.append(top)
    work = work[work["SNP"] != top["SNP"]]

if signals:
    out = pd.DataFrame(signals)
else:
    out = pd.DataFrame(columns=["SNP", "P", "iteration"])

out.to_csv(out_tsv, sep="\t", index=False)

with open(diag_json, "w", encoding="utf-8") as h:
    json.dump({"method": "cojo_placeholder", "p_cutoff": p_cutoff, "iterations_run": len(signals)}, h, indent=2)
PY

echo "ok" > "$DONE_FILE"
