#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1purity -o %x_%j.out
# D1h: label purity of collapsed clusters (study, country), seven species, default collapse order.
# usage: REPO=<worktree root> sbatch label_purity.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\tcolumn\tclusters\tmulti_member\tmixed_clusters\tmixed_share_of_multi\tsamples_differing_from_rep\n" > "$OUTD/label_purity.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1h_$slug; rm -rf "$W"; mkdir -p "$W"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "id","study","country"} NR>1{print $1,$2,$10}' "$P/$slug/runs.tsv" > "$W/labels.tsv"
  "$BIN" run --from-sketches "$P/$slug/sketches" --no-pca --threads 8 --floor-samples 1000 --collapse mean -o "$W/out" -q 2> "$W/log" || { echo "FAILED $slug" >&2; cat "$W/log" >&2; continue; }
  ( cd "$OUTD" && python3 label_purity.py "$slug" "$W/out/collapse.tsv" "$W/labels.tsv" ) >> "$OUTD/label_purity.tsv"
done
