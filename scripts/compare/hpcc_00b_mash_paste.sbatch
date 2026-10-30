#!/bin/bash -l
#SBATCH -J cmp_mashpaste
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/cmp_mashpaste.%j.log
# Tool comparison (#46): the mash reference from the BFD pipeline's per-genome sketches (k=21, 10000 hashes: the parameters wanted here),
# for the same 7.4k reference genomes. `mash sketch -l` over all 7.4k FASTA files did its work (11.5k CPU-s) but exited with status 1 when
# writing the combined sketch, so the existing sketches are pasted instead.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load mash/2.3
B=$BFD_GENOMES/../results/mash_sketch; cd $GR/compare/ref
while read id; do f=$B/$id.fa.msh; [ -s "$f" ] && echo "$f"; done < refids.txt > mashlist.txt
wc -l refids.txt mashlist.txt
rm -f mash_refs.msh
/usr/bin/time -v -o mash_paste.time mash paste mash_refs mashlist.txt
mash info -H mash_refs.msh | head -8
ls -la mash_refs.msh
