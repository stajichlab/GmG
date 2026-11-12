#!/bin/bash
# P6: where do the time and memory go in a chunk job? Re-runs the two largest B. bruxellensis runs of the phase 1 chunk
# that peaked at 20.0 of 20 GiB (task 0), plus two smaller runs, with `popseq fetch --profile`.
# Run on the HPCC from a clone of GmG that has scripts/popsra_fetch_chunks.slurm with --profile (GmG #186):
#   POPSEQ_ROOT=<PopSeq checkout at or after #45> bash bench/P6_profile_big/run.sh
# Everything it writes goes under $OUT (default ~/projects/GmG/bench/P6_profile_big/out): purgeable. Temp files: $SCRATCH.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GMG_ROOT=$(cd "$HERE/../.." && pwd)
: "${POPSEQ_ROOT:?}"
OUT=${OUT:-$HOME/projects/GmG/bench/P6_profile_big/out}
mkdir -p "$OUT"
# Build popsketch inside the job (node-local $SCRATCH), then run the unchanged chunk script on 2 chunks of 2 runs.
cat > "$OUT/job.sh" <<JOB
#!/bin/bash
set -euo pipefail
cd "$GMG_ROOT" && cargo build --release -q -p popsketch --target-dir "\$SCRATCH/target"
export POPSKETCH="\$SCRATCH/target/release/popsketch"
exec bash "$GMG_ROOT/scripts/popsra_fetch_chunks.slurm"
JOB
sbatch --parsable --array=0-1 --mem 32G -J p6profile -o "$OUT/slurm-%a.out" \
  --export=ALL,POPSEQ_ROOT="$POPSEQ_ROOT",GMG_ROOT="$GMG_ROOT",WORK="$OUT/work",MANIFEST="$HERE/manifest.csv",CHUNK=2,PAR=2,CONCURRENCY=2 \
  "$OUT/job.sh"
