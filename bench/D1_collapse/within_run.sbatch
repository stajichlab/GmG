#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1within -o %x_%j.out
# D1i: what does a label-aware collapse cost? Representatives under --collapse mean with no restriction, within study,
# within country, and within study+country; plus label purity of the clusters (should be 100% by construction).
# usage: REPO=<worktree root> sbatch within_run.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\twithin\tn\treps\tlargest_cluster\tsingletons\tmixed_study_clusters\n" > "$OUTD/within_run.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1i_$slug; rm -rf "$W"; mkdir -p "$W"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "id","study","country"} NR>1{print $1,$2,$10}' "$P/$slug/runs.tsv" > "$W/labels.tsv"
  for within in none study country study,country; do
    W2="$W/$(echo $within | tr , _)"
    ARGS=(); [ "$within" != none ] && ARGS=(--collapse-within "$within")
    "$BIN" run --from-sketches "$P/$slug/sketches" --labels "$W/labels.tsv" --permutations 0 --no-pca --threads 8 --floor-samples 1000 --collapse mean "${ARGS[@]}" -o "$W2" -q 2> "$W2.log" || { echo "FAILED $slug $within" >&2; cat "$W2.log" >&2; continue; }
    mixed=$(cd "$OUTD" && python3 label_purity.py "$slug" "$W2/collapse.tsv" "$W/labels.tsv" | awk -F'\t' '$2=="study"{print $5}')
    awk -F'\t' -v s=$slug -v w=$within -v m="$mixed" 'NR>1 { n++; if ($5=="true") { reps++; if ($3>mx) mx=$3; if ($3==1) one++ } } END { printf "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", s, w, n, reps, mx, one, m }' "$W2/collapse.tsv" >> "$OUTD/within_run.tsv"
  done
done
