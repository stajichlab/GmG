#!/bin/bash
# pick 12 Cocci and 12 Bd samples (seeded), run the split-half sketching 8 at a time
set -uo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
cd $P; mkdir -p splithalf
python3 - <<PY
import pandas as pd
rows = []
for ds, qc, man in (("cocci", "popgen_cocci/final_s200_rc/qc.tsv", "popgen_cocci/manifest.tsv"), ("bd", "popgen_bd/out_s200_qc_rc/qc.tsv", "popgen_bd/manifest.tsv")):
    q = pd.read_csv(qc, sep="\t"); q = q[q.status == "ok"]
    m = pd.read_csv(man, sep="\t").set_index("id")
    for i in q.sample(12, random_state=1).id: rows.append((ds, i, m.r1[i]))
pd.DataFrame(rows).to_csv("splithalf/tasks.tsv", sep="\t", header=False, index=False)
PY
cat splithalf/tasks.tsv | tr '\t' ' ' | xargs -P 8 -L1 bash popgen_cocci/20_splithalf_task.sh
echo "split-half sketches: $(ls splithalf/*/sk/*.gsk | wc -l)"
