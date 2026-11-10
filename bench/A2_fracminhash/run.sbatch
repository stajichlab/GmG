#!/bin/bash
#SBATCH -p short -N1 -c4 --mem=16G -t 1:00:00 -J a2nest -o %x_%j.out
# A2: prediction table (predict.py) and the nested scaled-200 -> scaled-1000 check on real sketches (nested_check.py).
# usage: REPO=<worktree root> sbatch run.sbatch   (existing sketches only)
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/A2_fracminhash"
module load miniconda3
python3 "$OUTD/predict.py" > "$OUTD/prediction.tsv"
printf "slug\tn_pairs\tmedian_one_minus_ani\trms_z\tfrac_abs_z_gt_2\tmedian_hashes_200\tmedian_hashes_1000\tbytes_ratio\n" > "$OUTD/nested_check.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  python3 "$OUTD/nested_check.py" "$slug" "$P/$slug/sketches" 80 3 >> "$OUTD/nested_check.tsv" || echo "FAILED $slug" >&2
done
