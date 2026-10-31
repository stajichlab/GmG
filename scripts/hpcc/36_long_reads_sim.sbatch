#!/bin/bash -l
#SBATCH -J lr_sim
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=28G
#SBATCH -t 02:00:00
#SBATCH --array=0-95
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/36_lr_sim.%A_%a.log
# Long reads (#44), part 1: simulated ONT/PacBio-like reads at 8 error rates x 7 depths for the 24 genomes of data/validation/long_reads_sim_genomes.tsv.
# Array task = genome (24) x pair of error rates (4): 24 x 4 = 96 tasks. For every (genome, error) one 30x read stream is simulated (scripts/sim_long_reads.py,
# nested prefixes = the 7 depths 0.25 ... 30x), its realised error is measured by minimap2 alignment (scripts/long_reads_measure.py), and the 14 read sets
# are classified with the index of the release (full index, s10000) with the query's own genome hidden (--exclude-level self), and with its species or its genus
# hidden (default min-abund 2), min-abund 1 and min-abund auto (self only). The reads are deleted when the task is done; the results stay.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${LR_SCRIPTS:-$GR/long_reads/scripts}
GMGBIN=${GMGBIN:-$GR/bin/gmg-readscal}      # needs `--min-abund auto` and coverage_est (#32)
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
LIST=${LIST:-$GR/long_reads/long_reads_sim_genomes.tsv}
ERRS=(0.005 0.01 0.02 0.05 0.08 0.10 0.12 0.15)
DEPTHS="0.25,0.5,1,2,5,10,30"
i=$SLURM_ARRAY_TASK_ID; G=$((i / 4)); P=$((i % 4))
row=$(awk -F'\t' 'NR>1' $LIST | sed -n "$((G + 1))p")
gid=$(echo "$row" | cut -f1); sp=$(echo "$row" | cut -f2)
W=$GR/long_reads/sim/$gid/p$P; mkdir -p $W/reads; cd $W
echo "gmg: $($GMGBIN --version)  index: $INDEX  genome: $gid ($sp)  errors: ${ERRS[$((2*P))]} ${ERRS[$((2*P+1))]}"
g=$BFD_GENOMES/clean/$gid.purge.fasta.gz
printf "id\tpath\texpected\tref\n" > sheet.tsv
printf "label\tn_reads\tn_aligned\tmedian_len\tmean_len\terror_pct\tmismatch_pct\tindel_pct\n" > measured.tsv
for e in ${ERRS[$((2*P))]} ${ERRS[$((2*P+1))]}; do
  pre=reads/e$e
  [ -s ${pre}_d30.fastq ] || python3 $S/sim_long_reads.py $g $pre --err $e --depths $DEPTHS --seed $((1000 * G + i)) || { echo "simulation failed"; exit 1; }
  python3 $S/long_reads_measure.py e$e ${pre}_d5.fastq $g 300 >> measured.tsv
  for d in ${DEPTHS//,/ }; do
    printf "e%s_d%s\t%s/%s_d%s.fastq\t%s\t%s\n" $e $d $W ${pre} $d "$sp" $gid >> sheet.tsv
  done
done
run() {  # tag, extra args
  tag=$1; shift
  [ -s $tag/results.tsv ] || $GMGBIN --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o $tag --no-report -q --kind reads "$@"
}
run self_a2 --exclude-level self --min-abund 2
run self_a1 --exclude-level self --min-abund 1
run self_aauto --exclude-level self --min-abund auto
run species_a2 --exclude-level species --min-abund 2
run genus_a2 --exclude-level genus --min-abund 2
rm -rf reads
echo done
