#!/bin/bash
#SBATCH -p stajichlab -N1 -n1 -c8 --mem=16G -t 3:00:00 -J thrsweep -o logs/%x_%A.out
# #159 step 3: how many peaks change with the peak-window share (PEAK_WINDOW_SHARE in crates/gmg-core/src/coverage.rs, 0.3 now)?
# Builds one popsketch per threshold from a copy of the source (SRC, a GmG checkout at or after PR #158) in $SCRATCH, runs each over every species' finished sketches,
# writes thr_<value>/<slug>/qc.tsv under the working directory. Evaluate with 14_threshold_eval.py.
# usage: SRC=<GmG checkout> sbatch 13_threshold_sweep.sbatch      (run from popsra)
set -uo pipefail
: "${SRC:?GmG checkout}"
SLUGS="b_cinerea c_beticola c_neoformans n_discreta p_nodorum p_oryzae z_tritici"
cp -r "$SRC" $SCRATCH/gmg_src && cd $SCRATCH/gmg_src && rm -rf target
for t in 0.20 0.25 0.30 0.35 0.40; do
  sed -i "s/^const PEAK_WINDOW_SHARE: f64 = .*;/const PEAK_WINDOW_SHARE: f64 = $t;/" crates/gmg-core/src/coverage.rs
  grep -q "PEAK_WINDOW_SHARE: f64 = $t;" crates/gmg-core/src/coverage.rs || { echo "FAIL sed $t"; exit 1; }
  cargo build --release -q -p popsketch --target-dir $SCRATCH/target_$t || { echo "FAIL build $t"; exit 1; }
done
cd "$SLURM_SUBMIT_DIR"
for t in 0.20 0.25 0.30 0.35 0.40; do for s in $SLUGS; do
  $SCRATCH/target_$t/release/popsketch run --from-sketches $s/sketches -o thr_$t/$s --threads 8 -q >/dev/null 2>>thr_$t.err || echo "FAIL run $t $s"
done; done
echo done
