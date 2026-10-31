#!/bin/bash -l
#SBATCH -J lr_flag
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=28G
#SBATCH -t 02:00:00
#SBATCH --array=0-5
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/39_lr_flag.%A_%a.log
# Long reads (#44): the `high_error_reads` flag as built into gmg (private build 38), on re-simulated reads of 6 of the 24 genomes (every 4th), 8 error rates,
# depths 1 ... 30x, own genome hidden, min-abund 2: flag rate and error_est per (error, depth). Output $GR/long_reads/flag/<genome>/results.{tsv,jsonl}.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${LR_SCRIPTS:-$GR/long_reads/scripts}
GMGBIN=${GMGBIN:-$GR/long_reads/bin/gmg}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
LIST=${LIST:-$GR/long_reads/long_reads_sim_genomes.tsv}
ERRS="0.005 0.01 0.02 0.05 0.08 0.10 0.12 0.15"; DEPTHS="1,2,5,10,30"
G=$((SLURM_ARRAY_TASK_ID * 4))
row=$(awk -F'\t' 'NR>1' $LIST | sed -n "$((G + 1))p"); gid=$(echo "$row" | cut -f1); sp=$(echo "$row" | cut -f2)
W=$GR/long_reads/flag/$gid; mkdir -p $W/reads; cd $W
echo "gmg: $($GMGBIN --version) genome $gid"
g=$BFD_GENOMES/clean/$gid.purge.fasta.gz
printf "id\tpath\texpected\tref\n" > sheet.tsv
for e in $ERRS; do
  python3 $S/sim_long_reads.py $g reads/e$e --err $e --depths $DEPTHS --seed $((777 + G)) || exit 1
  for d in ${DEPTHS//,/ }; do printf "e%s_d%s\t%s/reads/e%s_d%s.fastq\t%s\t%s\n" $e $d $W $e $d "$sp" $gid >> sheet.tsv; done
done
$GMGBIN --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o out --no-report -q --kind reads --exclude-level self --min-abund 2
rm -rf reads; echo done
