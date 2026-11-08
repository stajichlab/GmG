#!/bin/bash
#SBATCH -p short -N1 -c8 --mem=64G -t 2:00:00 -J a1mem -o %x_%j.out
# A1 measurement: peak RSS and wall time of `popsketch run --from-sketches` against the table in README.md.
# Sets are symlink folders in $SCRATCH (replicas get new ids; they are exact copies, so distances between a sample and its
# replica are 0: this measures memory and the shape-driven cost, not the cost on structureless real data).
# Only results.tsv is kept. usage: REPO=<worktree root> sbatch measure.sbatch
set -uo pipefail
REPO=${REPO:?set REPO to the repository root}
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra
OUTD=$(cd "$(dirname "$0")" && pwd)
export CARGO_TARGET_DIR=$SCRATCH/a1_target
( cd "$REPO" && cargo build --release -p popsketch ) || exit 1
BIN=$CARGO_TARGET_DIR/release/popsketch
W=$SCRATCH/a1_work; rm -rf "$W"; mkdir -p "$W"

mkset() { # name reps slug...
  local name=$1 reps=$2; shift 2
  mkdir -p "$W/$name"
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
ALL="z_tritici c_neoformans c_beticola p_nodorum b_cinerea n_discreta p_oryzae"
mkset n264  1 c_neoformans
mkset n2002 1 $ALL
mkset n4004 2 $ALL
mkset n6006 3 $ALL

printf "set\tvariant\tn_input\tin_analysis\tmax_rss_gb\twall_s\texit\n" > "$OUTD/results.tsv"
run() { # set variant extra-args...
  local set=$1 var=$2; shift 2
  local n=$(ls "$W/$set"/*.json | wc -l)
  /usr/bin/time -v -o "$W/time.txt" timeout 1500 "$BIN" run --from-sketches "$W/$set" -o "$W/out_$set$var" --threads 8 --max-samples 7000 --flagged mark -q "$@" > "$W/log_$set$var.txt" 2>&1
  local ec=$?
  local rss=$(awk '/Maximum resident/{printf "%.2f",$6/1048576}' "$W/time.txt")
  local wall=$(awk '/Elapsed/{n=split($NF,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; print s}' "$W/time.txt")
  local ina=$(python3 -c "import json;print(json.load(open('$W/out_$set$var/popsketch.json'))['samples']['in_analysis'])" 2>/dev/null || echo NA)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$set" "$var" "$n" "$ina" "$rss" "$wall" "$ec" >> "$OUTD/results.tsv"
  rm -rf "$W/out_$set$var"
}
run n264  _nopca --no-pca
run n264  _pca
run n2002 _nopca --no-pca
run n2002 _pca
run n4004 _nopca --no-pca
run n6006 _nopca --no-pca
run n4004 _pca
cat "$OUTD/results.tsv"
