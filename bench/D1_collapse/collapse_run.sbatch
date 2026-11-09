#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1run -o %x_%j.out
# D1f: `popsketch run --collapse` on the seven species (existing sketches). Representatives per rule, and the largest cluster.
# usage: REPO=<worktree root> sbatch collapse_run.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\trule\tn\treps\tlargest_cluster\tsingletons\n" > "$OUTD/collapse_run.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  for rule in mean max median; do
    W=$SCRATCH/d1f_${slug}_$rule; rm -rf "$W"
    "$BIN" run --from-sketches "$P/$slug/sketches" -o "$W" --no-pca --threads 8 --floor-samples 1000 --collapse $rule -q 2> "$W.log" || { echo "FAILED: $slug $rule" >&2; cat "$W.log" >&2; continue; }
    [ -s "$W/collapse.tsv" ] || { echo "FAILED: no collapse.tsv $slug $rule" >&2; continue; }
    awk -F'\t' -v s=$slug -v r=$rule 'NR>1 { n++; if ($5=="true") { reps++; if ($3>mx) mx=$3; if ($3==1) one++ } } END { printf "%s\t%s\t%d\t%d\t%d\t%d\n", s, r, n, reps, mx, one }' "$W/collapse.tsv" >> "$OUTD/collapse_run.tsv"
  done
done
