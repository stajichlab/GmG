#!/bin/bash -l
#SBATCH -J lr_pilot
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=0-79%12
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_lr_pilot.%A_%a.log
# Long reads (#44), part 2, step 1: measure the real error rate of every candidate run (data/validation/long_reads_candidates.tsv) from the first 4 MB
# of its FASTQ (HTTPS Range request to the ENA file server, ~100-600 reads) aligned to the run's own chromosome-level assembly (mappy).
# The chemistry in the candidate list is only a guess from the metadata; this gives the class used to choose the 20 runs. 4 MB x 80 = 0.3 GB downloaded.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${LR_SCRIPTS:-$GR/long_reads/scripts}
W=$GR/long_reads/pilot; mkdir -p $W/raw $W/out; cd $W
CAND=${CAND:-$GR/long_reads/long_reads_candidates.tsv}
row=$(awk -F'\t' 'NR>1' $CAND | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
run=$(echo "$row" | cut -f1); gid=$(echo "$row" | cut -f11); url=$(echo "$row" | cut -f22)
f=raw/$run.head.fastq.gz
[ -s $f ] || curl -s -f -r 0-4194303 -o $f https://$url || { echo "$run download failed"; exit 1; }
g=$BFD_GENOMES/clean/$gid.purge.fasta.gz   # the cleaned (contaminant-purged) assembly of the reference genome
python3 $S/long_reads_measure.py $run $f $g 300 > out/$run.tsv 2> out/$run.err
cat out/$run.tsv
