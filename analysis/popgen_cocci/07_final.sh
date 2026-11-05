#!/bin/bash
# outliers dropped (TX13 Fusarium; NM_7898 Histoplasma; CA10, CA12 Lecanicillium), the two swaps relabelled
# to the species their k-mers and nearest neighbours say, then: all-Coccidioides + each species on its own
set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci
EX=TX13,NM_7898,CA10,CA12
RL="485B-1_L_OLD_CPA0023=Coccidioides posadasii,B3476=Coccidioides immitis"
for s in 1000 200; do
  python3 02_popdist.py sk$s manifest.tsv final_s$s --exclude $EX --relabel "$RL" --flag flagged.txt --ref-hashes refs_cocci_s$s.npy 2>&1 | grep -v -i warn
  python3 03_plots.py final_s$s || true
done
for sp in immitis posadasii; do
  python3 02_popdist.py sk200 manifest.tsv final_s200_$sp --exclude $EX --relabel "$RL" --only-label "Coccidioides $sp" \
     --flag flagged.txt --ref-hashes refs_cocci_s200.npy --cut 0.004 2>&1 | grep -v -i warn
done
