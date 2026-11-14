#!/bin/bash
# P9: background sets, first three (the Coccidioides pair (100 runs) and the A. fumigatus pool (150), analysis/popsets), with popseq fetch --profile. Runs over 30 Gb of bases are in manifest_deferred.csv, not here (no size ceiling in popseq fetch yet, PopSeq #48).
# Run on the HPCC from a clone of GmG:
#   POPSEQ_ROOT=<PopSeq checkout at v0.2.0 or later> ARRAY=0-7 bash bench/P9_background_sets/run.sh     # 8 chunks of 32 runs (249 runs)
# Sampling is the popseq default (auto); coverage 25. Output under $OUT (purgeable); temp files in node-local $SCRATCH.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GMG_ROOT=$(cd "$HERE/../.." && pwd)
: "${POPSEQ_ROOT:?}" "${ARRAY:?array range, e.g. 0 or 1-9}"
OUT=${OUT:-$HOME/projects/GmG/bench/P9_background_sets/out}
MEM=${MEM:-8G} # the process peaks at about 0.2 GiB per fetch (P6, P7); the rest is page cache
mkdir -p "$OUT"
cat > "$OUT/job.sh" <<JOB
#!/bin/bash
set -euo pipefail
cd "$GMG_ROOT" && cargo build --release -q -p popsketch --target-dir "\$SCRATCH/target"
export POPSKETCH="\$SCRATCH/target/release/popsketch"
exec bash "$GMG_ROOT/scripts/popsra_fetch_chunks.slurm"
JOB
sbatch --parsable --array="$ARRAY%4" -p "${PARTITION:-stajichlab}" -c 16 -t 08:00:00 --mem "$MEM" -J p9background -o "$OUT/slurm-%a.out" \
  --export=ALL,POPSEQ_ROOT="$POPSEQ_ROOT",GMG_ROOT="$GMG_ROOT",WORK="$OUT/work",MANIFEST="$HERE/manifest.csv",CHUNK=32,PAR=4,CONCURRENCY=16,SAMPLING=auto \
  "$OUT/job.sh"
