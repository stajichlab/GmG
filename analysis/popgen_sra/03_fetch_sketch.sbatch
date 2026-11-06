#!/bin/bash
#SBATCH -p short -N1 -n1 -c4 --mem=12G -t 2:00:00 -J popsra -o logs/%x_%A_%a.out
# One SRA run per array task: parallel-fastq-dump to scratch, keep every Nth read pair (about 1 Gb of bases in all: old SRA FASTQs
# are position-sorted, so the first reads would come from one region), sketch it into the species' folder of sketches, delete the reads.
# Only the sketch (about 8 MB with its split halves) lands under /bigdata.
# usage: N=$(( $(wc -l < <SLUG>/runs_ftp.tsv) - 1 )); sbatch --array=1-$N%12 --export=ALL,SLUG=<slug> 03_fetch_sketch.sbatch
set -uo pipefail
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 2>/dev/null || true
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=$P/popsketch_build/target/release/popsketch
cd $P/popsra
tsv=$SLUG/runs_ftp.tsv
run=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' $tsv)
bases=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==1{for(i=1;i<=NF;i++)if($i=="base_count")c=i} NR==n{print $c}' $tsv)
sk=$SLUG/sketches
mkdir -p $sk
[[ -s $sk/$run.json ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/popsra.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
parallel-fastq-dump --sra-id $run --threads 4 --split-files --gzip --tmpdir $work -O $work/raw >$work/dump.log 2>&1 || { echo "FAIL dump $run"; tail -5 $work/dump.log; exit 1; }
step=$(( (bases + 999999999) / 1000000000 )); [[ $step -lt 1 ]] && step=1
for m in 1 2; do
  f=$work/raw/${run}_$m.fastq.gz
  [[ -s $f ]] || { echo "FAIL missing mate $m of $run"; exit 1; }
  zcat $f | awk -v s=$step 'NR%4==1{k=((NR-1)/4)%s==0} k' | gzip > $work/sub_$m.fastq.gz
done
$BIN sketch --id $run -o $sk $work/sub_1.fastq.gz $work/sub_2.fastq.gz || { echo "FAIL sketch $run"; exit 1; }
