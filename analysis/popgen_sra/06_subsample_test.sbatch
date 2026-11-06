#!/bin/bash
#SBATCH -p short -N1 -n1 -c4 --mem=12G -t 2:00:00 -J subtest -o logs/%x_%A_%a.out
# Does a subsample give the same sketch-level answer as the whole run? One run per array task (row of subtest/tasks.tsv: slug, run, base_count).
# Fetches the full run once and sketches it five ways into subtest/<slug>/<variant>/ (one folder of sketches per variant):
#   full      all reads
#   stride1g  every Nth pair to about 1 Gb of bases (what the population sets use; copied from popsra/<slug>/sketches)
#   stride300 every Nth pair to about 300 Mb
#   head300   the first 300 Mb of reads
#   head100   the first 100 Mb of reads
# usage: sbatch --array=1-N%8 06_subsample_test.sbatch
set -uo pipefail
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 2>/dev/null || true
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=$P/popsketch_build/target/release/popsketch
cd $P
read -r slug run bases < <(awk -F'\t' -v n=$SLURM_ARRAY_TASK_ID 'NR==n' subtest/tasks.tsv)
T=subtest/$slug; mkdir -p $T/stride1g
[[ -s $T/stride1g/$run.json ]] || cp popsra/$slug/sketches/$run.* popsra/$slug/sketches/params.txt $T/stride1g/ 2>/dev/null || true
[[ -s $T/head100/$run.json ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/subtest.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
parallel-fastq-dump --sra-id $run --threads 4 --split-files --gzip --tmpdir $work -O $work/raw >$work/dump.log 2>&1 || { echo "FAIL dump $run"; tail -5 $work/dump.log; exit 1; }
r1=$work/raw/${run}_1.fastq.gz; r2=$work/raw/${run}_2.fastq.gz
[[ -s $r1 && -s $r2 ]] || { echo "FAIL mates $run"; exit 1; }
sk() { $BIN sketch --id $run -o $T/$1 --max-bases $2 "${@:3}" || { echo "FAIL sketch $1 $run"; exit 1; }; }
sk full 0 $r1 $r2
sk head300 300000000 $r1 $r2
sk head100 100000000 $r1 $r2
step=$(( (bases + 299999999) / 300000000 )); [[ $step -lt 1 ]] && step=1
for m in 1 2; do zcat $work/raw/${run}_$m.fastq.gz | awk -v s=$step 'NR%4==1{k=((NR-1)/4)%s==0} k' | gzip > $work/s_$m.fastq.gz; done
sk stride300 300000000 $work/s_1.fastq.gz $work/s_2.fastq.gz
echo "done $run step=$step"
