#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 2:00:00 -J d1order -o %x_%j.out
# D1g: stability of --collapse across visiting orders (the default noise order and 10 random seeds), seven species.
# For each: reps, cluster count k, agreement of the clusters (via representatives) with the uncollapsed clusters, and the
# label tests (study, country) on the representatives. Existing sketches only.
# usage: REPO=<worktree root> sbatch order_check.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\tseed\tn\treps\tk_clusters\tari_vs_baseline\tstudy_r2\tstudy_p\tcountry_r2\tcountry_p\n" > "$OUTD/order_check.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1g_$slug; rm -rf "$W"; mkdir -p "$W"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "id","study","country"} NR>1{print $1,$2,$10}' "$P/$slug/runs.tsv" > "$W/labels.tsv"
  COMMON=(run --from-sketches "$P/$slug/sketches" --labels "$W/labels.tsv" --permutations 199 --no-pca --threads 8 --floor-samples 1000 -q)
  "$BIN" "${COMMON[@]}" -o "$W/base" 2> "$W/base.log" || { echo "FAILED baseline $slug" >&2; cat "$W/base.log" >&2; continue; }
  for seed in noise 1 2 3 4 5 6 7 8 9 10; do
    SEEDARG=(); [ "$seed" != noise ] && SEEDARG=(--collapse-seed "$seed")
    "$BIN" "${COMMON[@]}" --collapse mean "${SEEDARG[@]}" -o "$W/s$seed" 2> "$W/s$seed.log" || { echo "FAILED $slug $seed" >&2; cat "$W/s$seed.log" >&2; continue; }
    ( cd "$OUTD" && python3 order_check.py "$slug" "$seed" "$W/base" "$W/s$seed" ) >> "$OUTD/order_check.tsv"
  done
done
