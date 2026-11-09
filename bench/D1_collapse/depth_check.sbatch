#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=32G -t 1:00:00 -J d1depth -o %x_%j.out
# D1d: how far is a sample from its own lower-depth subsample? (12 runs x 1, 2, 4 Gb, z_tritici and c_neoformans),
# beside the split-half floor at each depth. Existing sketches only; no downloads.
# usage: REPO=<worktree root> sbatch depth_check.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
D=/bigdata/stajichlab/jstajich/projects/GmG_popgen/depth
OUTD="$REPO/bench/D1_collapse"
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/d1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
printf "slug\trun\tdepth_a\tdepth_b\tdistance\n" > "$OUTD/depth_self.tsv"
printf "slug\tdepth\tn\tfloor_median\tfloor_p90\n" > "$OUTD/depth_floor.tsv"
for slug in c_neoformans z_tritici; do
  W=$SCRATCH/d1d_$slug; rm -rf "$W"; mkdir -p "$W/pool"
  cp "$D/$slug/s1g/params.txt" "$W/pool/params.txt"
  for dep in s1g s2g s4g; do
    "$BIN" run --from-sketches "$D/$slug/$dep" -o "$W/out_$dep" --flagged mark --min-cov 0 --no-pca --no-tree --floor-samples 100 --threads 8 -q || continue
    python3 -c "import json;f=json.load(open('$W/out_$dep/popsketch.json'))['noise_floor'];print('$slug','$dep',f['n'],f['median'],f['p90'],sep='\t')" >> "$OUTD/depth_floor.tsv"
    for j in "$D/$slug/$dep"/*.json; do
      id=$(basename "$j" .json); nid=${id}__$dep
      for part in "" .half_a .half_b; do ln -sf "$D/$slug/$dep/$id$part.gsk" "$W/pool/$nid$part.gsk"; done
      sed "s/\"id\":\"$id\"/\"id\":\"$nid\"/" "$j" > "$W/pool/$nid.json"
    done
  done
  "$BIN" run --from-sketches "$W/pool" -o "$W/out_pool" --flagged mark --min-cov 0 --no-pca --no-tree --floor-samples 100 --threads 8 -q || continue
  ( cd "$OUTD" && python3 depth_check.py "$slug" "$W/out_pool/ani.tsv" ) >> "$OUTD/depth_self.tsv"
done
