#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=64G -t 2:00:00 -J a1time -o %x_%j.out
# Per-stage timings of `popsketch run --from-sketches` (needs the stage-timings build) to find which stage makes n = 4,004 take over
# 25 minutes when n = 2,002 takes 2. Same symlink sets as measure.sbatch. Keeps only timings.txt.
# usage: REPO=<worktree root> sbatch timings.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD="$REPO/bench/A1_measure"
export CARGO_TARGET_DIR=$SCRATCH/a1t_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
W=$SCRATCH/a1t_work; rm -rf "$W"; mkdir -p "$W"
ALL="z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae"
mkset() { # name reps slug...
  local name=$1 reps=$2; shift 2
  mkdir -p "$W/$name"; cp "$P/$1/sketches/params.txt" "$W/$name/params.txt"
  for slug in "$@"; do
    for j in "$P/$slug"/sketches/*.json; do
      id=$(basename "$j" .json)
      for ((r=0; r<reps; r++)); do
        nid=$id; ((r>0)) && nid=${id}_r$r
        for part in "" .half_a .half_b; do ln -s "$P/$slug/sketches/$id$part.gsk" "$W/$name/$nid$part.gsk"; done
        sed "s/\"id\":\"$id\"/\"id\":\"$nid\"/" "$j" > "$W/$name/$nid.json"
      done
    done
  done
}
mkset n1001 1 c_neoformans c_beticola p_nodorum b_cinerea   # about 1,200 distinct real samples of 4 sets
mkset n2002 1 $ALL
mkset n4004 2 $ALL
: > "$OUTD/timings.txt"
for set in n1001 n2002 n4004; do
  echo "== $set ($(ls "$W/$set"/*.json | wc -l) samples), --no-pca" >> "$OUTD/timings.txt"
  timeout 2400 "$BIN" run --from-sketches "$W/$set" -o "$W/out_$set" --threads 8 --max-samples 7000 --flagged mark --no-pca > "$W/log_$set.txt" 2>&1
  echo "exit $?" >> "$W/log_$set.txt"
  grep -E "^stage|^exit|samples, " "$W/log_$set.txt" >> "$OUTD/timings.txt"
  rm -rf "$W/out_$set"
done
cat "$OUTD/timings.txt"
