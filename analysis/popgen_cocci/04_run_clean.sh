#!/bin/bash
# decontaminated, Fusarium-free runs at both resolutions
set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci
for s in 1000 200; do
  python3 02_popdist.py sk$s manifest.tsv clean_s$s --exclude TX13 --flag flagged.txt --ref-hashes refs_cocci_s$s.npy 2>&1 | grep -v -i warn
  python3 03_plots.py clean_s$s || true
done
