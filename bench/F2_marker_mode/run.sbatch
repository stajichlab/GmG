#!/bin/bash -l
#SBATCH -p epyc -N1 -c16 --mem=64G -t 6:00:00 -J f2markers -o %x_%j.out
# F2: GmG marker mode on FusariumID. Held-out evaluation (20% of strains as queries, one file per marker sequence) and
# the full marker index. Needs the parsed FusariumID records of F1 ($FID = .../F1_fusariumid_build/fid).
# usage: REPO=<checkout> W=<build dir> FID=<fid dir> sbatch run.sbatch
# Keeps in $W: holdout/{meta.tsv,queries_truth.tsv}, markers_holdout.gmgi, run_holdout/results.tsv, full/meta.tsv,
# markers_full.gmgi; sketches and the gmg build live in $SCRATCH. Results tables are copied to the repo by hand.
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}; FID=${FID:?set FID}
module load miniconda3
export CARGO_TARGET_DIR=$SCRATCH/f2_target
( cd "$REPO" && cargo build -q --release -p gmg-cli ) || exit 1
GMG=$CARGO_TARGET_DIR/release/gmg
T=$SLURM_CPUS_PER_TASK
S="$REPO/scripts/fusariumid_markers.py"
mkdir -p "$W"
cat > "$W/permissive.json" <<'EOF'
{"min_shared": [5, 5, 5, 5, 5, 5], "min_ani": [0, 0, 0, 0, 0, 0], "tentative_genus_min_ani": 2}
EOF
for set in holdout full; do
  frac=0.2; [ $set = full ] && frac=0
  rm -rf "$W/$set"; python3 -I "$S" prepare "$FID/records.tsv" "$FID/loci" "$W/$set" --holdout $frac --seed 1
  rm -rf "$SCRATCH/sk_$set"
  /usr/bin/time -v $GMG --threads $T -q sketch -k 21 -s 1 -o "$SCRATCH/sk_$set" "$W/$set"/refs/*.fa 2> "$W/${set}_sketch.time" || { tail "$W/${set}_sketch.time"; exit 1; }
  $GMG --threads $T build --meta "$W/$set/meta.tsv" --sketches "$SCRATCH/sk_$set" -o "$W/markers_$set.gmgi" \
    --taxonomy-version "FusariumID 2026-06-23 (strains)" --note "F2 marker mode, $set"
  rm -rf "$W/$set/refs"   # the reference FASTAs are a copy of FusariumID; the index and meta.tsv are enough
done
ls -la "$W"/markers_*.gmgi
rm -rf "$W/run_holdout"
/usr/bin/time -v $GMG --threads $T -q run "$W/holdout/queries" --index "$W/markers_holdout.gmgi" --thresholds "$W/permissive.json" \
  -o "$W/run_holdout" --no-report 2> "$W/run_holdout.time" || { tail -20 "$W/run_holdout.time"; exit 1; }
grep -E "Elapsed|Maximum resident" "$W/run_holdout.time"
python3 -I "$S" score "$W/run_holdout/results.tsv" "$W/holdout/queries_truth.tsv" "$REPO/crates/gmg-core/data/species_complexes.tsv" \
  "$W/calibration.tsv" --meta "$W/holdout/meta.tsv"
# the shipped threshold set (`--thresholds markers`): the held-out numbers as a user gets them
rm -rf "$W/run_holdout_markers"
$GMG --threads $T -q run "$W/holdout/queries" --index "$W/markers_holdout.gmgi" --thresholds markers -o "$W/run_holdout_markers" --no-report
python3 -I "$S" score "$W/run_holdout_markers/results.tsv" "$W/holdout/queries_truth.tsv" "$REPO/crates/gmg-core/data/species_complexes.tsv" \
  "$W/shipped.tsv" --meta "$W/holdout/meta.tsv"
