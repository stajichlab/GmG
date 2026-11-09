#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:30:00 -J d1collapse -o %x_%j.out
# D1: representative counts after collapsing near-identical samples, per species (existing sketches, no downloads).
# usage: REPO=<worktree root> sbatch run.sbatch      Only results.tsv is kept; ani.tsv files live in $SCRATCH.
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"   # not dirname $0: Slurm runs a spooled copy
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\tn\tthreshold\tleaders\tcomponents\tlargest_leader\n" > "$OUTD/results.tsv"
for slug in z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae; do
  W=$SCRATCH/d1_$slug; rm -rf "$W"
  "$BIN" run --from-sketches "$P/$slug/sketches" -o "$W" --no-pca --threads 8 -q || { echo "run failed: $slug" >&2; continue; }
  python3 "$OUTD/collapse_counts.py" "$slug" "$W/ani.tsv" >> "$OUTD/results.tsv"
done
