#!/bin/bash -l
#SBATCH -J cmp_refs
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --array=0-2
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/cmp_refs.%A_%a.log
# Tool comparison (#46): build each tool's reference of the representatives set (the genomes of the GmG reps index, decoys left out).
# usage: META=.../meta_reps.tsv sbatch hpcc_00_refs.sbatch      array: 0 skani, 1 mash (s=10000), 2 sourmash (k=21, scaled=10000)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
C=$GR/compare; E=$GR/envs/skani/bin; mkdir -p $C/ref; cd $C/ref
META=${META:-$GR/release/2026.10.3/meta_reps.tsv}
if [ ! -s reflist.txt ]; then
  tail -n +2 $META | cut -f1 | while read id; do f=$BFD_GENOMES/$id.fa.gz; [ -s "$f" ] && echo "$f"; done > reflist.txt
fi
wc -l reflist.txt
case $SLURM_ARRAY_TASK_ID in
0) /usr/bin/time -v -o skani_build.time $E/skani sketch -l reflist.txt -o skani_db -t 8 ;;
1) module load mash/2.3; /usr/bin/time -v -o mash_build.time mash sketch -p 8 -s 10000 -k 21 -l reflist.txt -o mash_refs ;;
2) module load sourmash/4.9
   (echo "name,genome_filename,protein_filename"; while read f; do b=$(basename $f .fa.gz); echo "$b,$f,"; done < reflist.txt) > sm_manifest.csv
   /usr/bin/time -v -o sourmash_build.time sourmash sketch fromfile sm_manifest.csv -p dna,k=21,scaled=10000 -o sm_refs.zip ;;
esac
ls -la
