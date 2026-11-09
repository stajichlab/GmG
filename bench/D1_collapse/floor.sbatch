#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=64G -t 1:30:00 -J d1floor -o %x_%j.out
# D1b: collapse thresholds taken from the split-half noise floor (1 - m * floor), then the pooled run of the
# representatives vs the pooled run of everything, to see the distance method, variable hashes and time of each.
# usage: REPO=<worktree root> sbatch floor.sbatch     (existing sketches only; no downloads)
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/D1_collapse"   # not dirname $0: Slurm runs a spooled copy
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
ALL="z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae"
printf "slug\tn\tfloor\tmultiple\tthreshold\tleaders\tcomponents\tlargest_leader\n" > "$OUTD/results_floor.tsv"
W=$SCRATCH/d1b; rm -rf "$W"; mkdir -p "$W/all" "$W/reps1" "$W/reps2" "$W/reps3" "$OUTD/logs"
for d in all reps1 reps2 reps3; do cp "$P/z_tritici/sketches/params.txt" "$W/$d/params.txt"; done
link() { # dir slug id
  for part in "" .half_a .half_b; do ln -sf "$P/$2/sketches/$3$part.gsk" "$W/$1/$3$part.gsk"; done
  ln -sf "$P/$2/sketches/$3.json" "$W/$1/$3.json"
}
for slug in $ALL; do
  "$BIN" run --from-sketches "$P/$slug/sketches" -o "$W/$slug" --no-pca --threads 8 2> "$W/$slug.log" || { echo "run failed: $slug" >&2; continue; }
  grep -h "distance method" "$W/$slug.log" | sed "s/^/$slug: /"
  F=$(python3 -c "import json;print(json.load(open('$W/$slug/popsketch.json'))['noise_floor']['value'])")
  python3 "$OUTD/collapse_counts.py" "$slug" "$W/$slug/ani.tsv" --floor "$F" --reps-out "$W/$slug.reps" >> "$OUTD/results_floor.tsv"
  for m in 1 2 3; do for id in $(cat "$W/$slug.reps.$m"); do link "reps$m" "$slug" "$id"; done; done
  for id in $(cut -f1 "$W/$slug/ani.tsv" | tail -n +2); do link all "$slug" "$id"; done
done
printf "set\tn_in\tdistance_line\twall_s\n" > "$OUTD/pooled.tsv"
for set in all reps1 reps2 reps3; do
  n=$(ls "$W/$set"/*.json | wc -l)
  t0=$SECONDS
  timeout 1500 "$BIN" run --from-sketches "$W/$set" -o "$W/out_$set" --no-pca --threads 8 2> "$W/$set.log"
  echo "exit $?" >> "$W/$set.log"; cp "$W/$set.log" "$OUTD/logs/pooled_$set.log"
  printf "%s\t%s\t%s\t%s\n" "$set" "$n" "$(grep -h 'distance method' "$W/$set.log" | tr '\n' ' ')" "$((SECONDS-t0))" >> "$OUTD/pooled.tsv"
done
