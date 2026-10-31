#!/bin/bash -l
#SBATCH -J cmp_mashpaste
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/cmp_mashpaste.%j.log
# Tool comparison (#46): the mash reference = one combined sketch of the same 7.4k reference genomes. The BFD pipeline already has per-genome
# sketches with the parameters wanted here (k=21, 10000 hashes) for most of them; the others are sketched here the same way (per file, in parallel).
# (`mash sketch -l` over all 7.4k FASTA files did its work, 11.5k CPU-s, but exited with status 1 at the final write.)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load mash/2.3
B=$BFD_GENOMES/../results/mash_sketch; cd $GR/compare/ref
mkdir -p mashsk; : > mashlist.txt; : > mashmissing.txt
while read id; do
  if [ -s "$B/$id.fa.msh" ]; then echo "$B/$id.fa.msh" >> mashlist.txt
  elif [ -s "mashsk/$id.msh" ]; then echo "$PWD/mashsk/$id.msh" >> mashlist.txt
  else echo "$id" >> mashmissing.txt; fi
done < refids.txt
wc -l refids.txt mashlist.txt mashmissing.txt
/usr/bin/time -v -o mash_extra.time bash -c "cat mashmissing.txt | xargs -P 8 -I{} sh -c 'mash sketch -s 10000 -k 21 -o mashsk/{} $BFD_GENOMES/{}.fa.gz >/dev/null 2>&1'"
: > mashlist.txt
while read id; do
  if [ -s "$B/$id.fa.msh" ]; then echo "$B/$id.fa.msh" >> mashlist.txt
  elif [ -s "mashsk/$id.msh" ]; then echo "$PWD/mashsk/$id.msh" >> mashlist.txt; fi
done < refids.txt
wc -l mashlist.txt
rm -f mash_refs.msh
/usr/bin/time -v -o mash_paste.time mash paste mash_refs $(cat mashlist.txt)
mash info -H mash_refs.msh | head -8
ls -la mash_refs.msh
