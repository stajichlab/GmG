#!/bin/bash
#SBATCH -p short -N1 -n1 -c4 --mem=12G -t 2:00:00 -J depth -o logs/%x_%A_%a.out
# How does resolution change with depth? One deep run (7 Gb of bases or more) per array task (row of depth/tasks.tsv: slug, run, base_count).
# Fetches the whole run once and sketches it at 1, 2, 4 and 6 Gb of bases (every pair kept with probability target/bases, evenly spaced, so the
# amounts are exact) and the whole run, each into depth/<slug>/<variant>/. Compare with 07_subsample_compare.py (VARIANTS=full,s6g,s4g,s2g,s1g).
# usage: sbatch --array=1-N%8 08_depth_series.sbatch
set -uo pipefail
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 2>/dev/null || true
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=$P/popsketch_build/target/release/popsketch
cd $P
read -r slug run bases < <(awk -F'\t' -v n=$SLURM_ARRAY_TASK_ID 'NR==n' depth/tasks.tsv)
T=depth/$slug
[[ -s $T/s1g/$run.json ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/depth.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
parallel-fastq-dump --sra-id $run --threads 4 --split-files --gzip --tmpdir $work -O $work/raw >$work/dump.log 2>&1 || { echo "FAIL dump $run"; tail -5 $work/dump.log; exit 1; }
r1=$work/raw/${run}_1.fastq.gz; r2=$work/raw/${run}_2.fastq.gz
[[ -s $r1 && -s $r2 ]] || { echo "FAIL mates $run"; exit 1; }
sk() { $BIN sketch --id $run -o $T/$1 --max-bases 0 "${@:2}" || { echo "FAIL sketch $1 $run"; exit 1; }; }
for g in 1 2 4 6; do
  f=$(awk -v t=${g}000000000 -v b=$bases 'BEGIN{printf "%.9f", (t>=b)?1:t/b}')
  for m in 1 2; do
    zcat $work/raw/${run}_$m.fastq.gz | awk -v f=$f 'NR%4==1{i++; k=(int(i*f)!=int((i-1)*f))} k' | gzip > $work/s${g}_$m.fastq.gz
  done
  sk s${g}g $work/s${g}_1.fastq.gz $work/s${g}_2.fastq.gz
  rm -f $work/s${g}_?.fastq.gz
done
sk full $r1 $r2
echo "done $run bases=$bases"
