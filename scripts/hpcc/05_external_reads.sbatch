#!/bin/bash -l
#SBATCH -J gmg_ext
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/05_external_reads.%j.log
# Stage 5 (genomes not in BFD) and stage 6 (reads) scoring against the BFD index.
#   $GR/external/sketch/*.gsk   sketches of non-BFD assemblies (made locally with `gmg sketch`)
#   $GR/reads/sim/*.gsk         sketches of simulated reads (made locally with `gmg sketch --abund`)
#   $GR/reads/cocci_r1.list     real Coccidioides R1 files, sketched here
set -uo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
mkdir -p $GR/external $GR/reads
for S in 10000 1000; do
  TH=$GR/bench/summary_s${S}.thresholds.json; [ -s $TH ] || TH=$GR/bench/summary_s10000.thresholds.json
  THARG=""; [ -s $TH ] && THARG="--thresholds $TH"
  IDX=$GR/index/bfd_s${S}.gmgi
  [ -s $IDX ] || { echo "no index $IDX"; continue; }
  # stage 5
  ls $GR/external/sketch/*.gsk > /dev/null 2>&1 && \
    $GMG guess --index $IDX $THARG $GR/external/sketch/*.gsk > $GR/external/guess_ext_s${S}.jsonl
  # stage 6a: simulated reads (abundance-filtered at query time)
  if ls $GR/reads/sim/*.gsk > /dev/null 2>&1; then
    for A in 1 2; do
      $GMG guess --index $IDX $THARG --min-abund $A $GR/reads/sim/*.gsk > $GR/reads/guess_sim_s${S}_a${A}.jsonl
    done
  fi
done
# stage 6b: real reads (first N reads of R1), sketched here then guessed
for N in 100000 500000 2000000; do
  $GMG --threads 4 sketch --abund --max-reads $N -k 21 -s 1000 -o $GR/reads/cocci_$N --list $GR/reads/cocci_r1.list 2> $GR/reads/cocci_$N.sketch.log
  for S in 10000 1000; do
    TH=$GR/bench/summary_s${S}.thresholds.json; [ -s $TH ] || TH=$GR/bench/summary_s10000.thresholds.json
    THARG=""; [ -s $TH ] && THARG="--thresholds $TH"
    for A in 1 2; do
      $GMG guess --index $GR/index/bfd_s${S}.gmgi $THARG --min-abund $A $GR/reads/cocci_$N/*.gsk > $GR/reads/guess_cocci_${N}_s${S}_a${A}.jsonl
    done
  done
done
ls -la $GR/external $GR/reads | head -40
