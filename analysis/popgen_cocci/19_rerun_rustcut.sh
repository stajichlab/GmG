#!/bin/bash
# rerun the Cocci and Bd analyses with gmg-core's abundance cutoff (no decontamination), suffix _rc
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
cd $P/popgen_cocci
EX=TX13,NM_7898,CA10,CA12
RL="485B-1_L_OLD_CPA0023=Coccidioides posadasii,B3476=Coccidioides immitis"
python3 02_popdist.py sk200 manifest.tsv final_s200_rc --exclude $EX --relabel "$RL" --cutoffs rust_cutoff_cocci_s200.tsv 2>&1 | grep -v -i warn
for sp in immitis posadasii; do
  python3 02_popdist.py sk200 manifest.tsv final_s200_${sp}_rc --exclude $EX --relabel "$RL" --only-label "Coccidioides $sp" --cutoffs rust_cutoff_cocci_s200.tsv --cut 0.004 2>&1 | grep -v -i warn
done
cd $P/popgen_bd
# half-content samples: estimated size < 17 Mb under the Rust cutoff
EXB=$(python3 - <<PY
import pandas as pd
r = pd.read_csv("rust_cutoff_bd_s200.tsv", sep="\t"); print(",".join(r[r.est_mb < 17].id))
PY
)
python3 02_popdist.py sk200 manifest.tsv out_s200_qc_rc --exclude "$EXB" --cutoffs rust_cutoff_bd_s200.tsv --cut 0.01 2>&1 | grep -v -i warn
