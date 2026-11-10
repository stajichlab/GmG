#!/bin/bash -l
#SBATCH -p epyc -N1 -c4 --mem=128G -t 4:00:00 -J f1mkana -o %x_%j.out
# F1 step 1b: re-run only marker_resolution.py on existing blast tables (after run_markers.sbatch).
# usage: REPO=<checkout> W=<build dir> sbatch run_marker_analysis.sbatch
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}
module load miniconda3
python3 -I "$REPO/bench/F1_fusariumid/marker_resolution.py" "$W/fid" "$W/blast" "$W/markers"
