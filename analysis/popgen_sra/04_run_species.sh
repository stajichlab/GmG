#!/bin/bash
# usage: 04_run_species.sh SLUG   (run from the popsra folder; BIN = popsketch binary)
# Full run on a folder of sketches, then a --min-floors sweep (no PCA) and a one-line summary per setting.
set -euo pipefail
SLUG=${1:?slug}; BIN=${BIN:-../popsketch_build/target/release/popsketch}; T=${THREADS:-8}
$BIN run --from-sketches $SLUG/sketches -o ${SLUG}_out --threads $T -q
for m in 1 2 3 5 8; do
  $BIN run --from-sketches $SLUG/sketches -o ${SLUG}_sweep/mf$m --threads $T --min-floors $m --no-pca -q
done
python3 - "$SLUG" <<'PY'
import json, sys
s = sys.argv[1]
for m in (1, 2, 3, 5, 8):
    d = json.load(open(f"{s}_sweep/mf{m}/popsketch.json"))
    g = d["clusters"]["gap"]
    print(s, "min_floors", m, g["verdict"], "top_jump_floors %.2f" % g["gap_over_floor"], "k_used", d["clusters"]["k_used"], "in_analysis", d["samples"]["in_analysis"])
PY
