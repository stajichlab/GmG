#!/bin/bash
# P7: phase 2, one full species (Pyricularia oryzae, 300 runs of panel pop_sra_ge10_cap300), with popseq fetch --profile.
# Run on the HPCC from a clone of GmG:
#   POPSEQ_ROOT=<PopSeq checkout at v0.2.0 or later> ARRAY=0 bash bench/P7_poryzae/run.sh     # the first chunk of 32 runs
#   POPSEQ_ROOT=... ARRAY=1-9 bash bench/P7_poryzae/run.sh                                     # the rest, once the first is read
# Sampling is the popseq default (auto); coverage 25. Output under $OUT (purgeable); temp files in node-local $SCRATCH.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GMG_ROOT=$(cd "$HERE/../.." && pwd)
: "${POPSEQ_ROOT:?}" "${ARRAY:?array range, e.g. 0 or 1-9}"
OUT=${OUT:-$HOME/projects/GmG/bench/P7_poryzae/out}
MEM=${MEM:-12G}
mkdir -p "$OUT"
cat > "$OUT/job.sh" <<JOB
#!/bin/bash
set -euo pipefail
cd "$GMG_ROOT" && cargo build --release -q -p popsketch --target-dir "\$SCRATCH/target"
export POPSKETCH="\$SCRATCH/target/release/popsketch"
exec bash "$GMG_ROOT/scripts/popsra_fetch_chunks.slurm"
JOB
sbatch --parsable --array="$ARRAY%4" -p intel,batch,epyc -c 16 -t 08:00:00 --mem "$MEM" -J p7poryzae -o "$OUT/slurm-%a.out" \
  --export=ALL,POPSEQ_ROOT="$POPSEQ_ROOT",GMG_ROOT="$GMG_ROOT",WORK="$OUT/work",MANIFEST="$HERE/manifest.csv",CHUNK=32,PAR=4,CONCURRENCY=16,SAMPLING=auto \
  "$OUT/job.sh"
