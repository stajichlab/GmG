#!/bin/bash
# re-checks with the Rust cutoff: (1) are the two label swaps still found when NOT relabelled? (2) SNP-tree clade tests
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci
cd $P
python3 02_popdist.py sk200 manifest.tsv norelabel_rc --exclude TX13,NM_7898,CA10,CA12 --cutoffs rust_cutoff_cocci_s200.tsv 2>&1 | grep -v -i warn
python3 - <<PY
import pandas as pd
lc = pd.read_csv("norelabel_rc/label_check.tsv", sep="\t")
print("label conflicts (no relabelling, Rust cutoff):"); print(lc[lc.conflict == 1].round(4).to_string())
PY
export SUF=_rc
python3 08_snp_concordance.py 3 2>&1 | grep -v -i warn
