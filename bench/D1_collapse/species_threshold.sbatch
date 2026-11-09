#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1sp -o %x_%j.out
# D1c: per-species collapse threshold from the split-half distances of ALL samples (--floor-samples 1000).
# usage: REPO=<worktree root> sbatch species_threshold.sbatch   (existing sketches only)
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\tn\tfloor_n\tmed\tp90\tnn_p10\tnn_p50\tnn_p90\treps_half_med\treps_med\treps_p90\n" > "$OUTD/species_threshold.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1c_$slug; rm -rf "$W"
  "$BIN" run --from-sketches "$P/$slug/sketches" -o "$W" --no-pca --threads 8 --floor-samples 1000 -q || { echo "run failed: $slug" >&2; continue; }
  ( cd "$OUTD" && python3 species_threshold.py "$slug" "$W/ani.tsv" "$W/popsketch.json" ) >> "$OUTD/species_threshold.tsv"
done
