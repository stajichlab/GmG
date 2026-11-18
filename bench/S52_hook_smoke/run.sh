#!/bin/bash
# S52: smoke test of the species check in the fetch hook (PopSeq#52) with the REAL popsketch species-match and the small index.
# Two runs, one array task: SRR30534096 (filed C. neoformans, reads A. fumigatus: must be refused, exit 14, sketch in refused/) and
# DRR829882 (clean C. neoformans: must be kept, verdict match or genus_only). Purgeable: everything is under $OUT.
#   POPSEQ_ROOT=<PopSeq checkout at v0.2.0 or later> bash bench/S52_hook_smoke/run.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GMG_ROOT=$(cd "$HERE/../.." && pwd)
: "${POPSEQ_ROOT:?}"
OUT=${OUT:-$HOME/projects/GmG/bench/S52_hook_smoke/out}
INDEX=${SPECIES_INDEX:-/bigdata/stajichlab/shared/projects/GmG_runs/index/bfd_s10000.gmgi}
mkdir -p "$OUT"
cat > "$OUT/job.sh" <<JOB
#!/bin/bash
set -euo pipefail
cd "$GMG_ROOT" && cargo build --release -q -p popsketch --target-dir "\$SCRATCH/target"
export POPSKETCH="\$SCRATCH/target/release/popsketch"
exec bash "$GMG_ROOT/scripts/popsra_fetch_chunks.slurm"
JOB
sbatch --parsable --array=0 -p "${PARTITION:-stajichlab}" -c 8 -t 02:00:00 --mem 12G -J s52smoke -o "$OUT/slurm-%a.out" \
  --export=ALL,POPSEQ_ROOT="$POPSEQ_ROOT",GMG_ROOT="$GMG_ROOT",WORK="$OUT/work",MANIFEST="$HERE/manifest.csv",CHUNK=32,PAR=2,CONCURRENCY=2,SAMPLING=auto,SPECIES_INDEX="$INDEX",SPECIES_MATCH=refuse \
  "$OUT/job.sh"
