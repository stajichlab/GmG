#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1pair -o %x_%j.out
# D1e: per-pair collapse threshold from each sample's own split-half distance (split_half.tsv), seven species.
# usage: REPO=<worktree root> sbatch per_pair.sbatch   (existing sketches only)
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\tn\tmedian\trule\treps\tlargest_cluster\n" > "$OUTD/per_pair.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1e_$slug; rm -rf "$W"
  "$BIN" run --from-sketches "$P/$slug/sketches" -o "$W" --no-pca --threads 8 --floor-samples 1000 -q || { echo "FAILED: $slug" >&2; continue; }
  [ -s "$W/split_half.tsv" ] || { echo "FAILED: no split_half.tsv for $slug" >&2; continue; }
  cp "$W/split_half.tsv" "$OUTD/split_half_$slug.tsv"
  ( cd "$OUTD" && python3 per_pair.py "$slug" "$W/ani.tsv" "$W/split_half.tsv" ) >> "$OUTD/per_pair.tsv"
done
