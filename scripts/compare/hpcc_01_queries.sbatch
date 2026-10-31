#!/bin/bash -l
#SBATCH -J cmp_queries
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/cmp_queries.%j.log
# Tool comparison (#46): the query files. s1 = post-freeze assemblies (downloaded), s2 = BFD copies of RefSeq truth-set genomes (linked),
# s3 = the stress inputs of #41 (stress_sim.py: completeness, fragmentation, contamination, negative controls).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
C=$GR/compare; Q=$C/q; mkdir -p $Q/s1 $Q/s2 $Q/s3 $C/sets $C/work $C/s3_files; cd $C/work
python3 $GMG_REPO/scripts/compare/make_sets.py $GMG_REPO $C/sets 5
tail -n +2 $C/sets/s2_truth.tsv | while IFS=$'\t' read id bfd rest; do ln -sf $BFD_GENOMES/$bfd.fa.gz $Q/s2/$id.fa.gz; done
tail -n +2 $C/sets/s1_post_freeze.tsv | cut -f1 > s1.acc
datasets download genome accession --inputfile s1.acc --include genome --filename s1.zip --no-progressbar && unzip -oq s1.zip -d s1x
for f in s1x/ncbi_dataset/data/GC*/*.fna; do gzip -c "$f" > $Q/s1/$(basename $(dirname "$f")).fa.gz; done
datasets download genome accession --inputfile $GMG_REPO/data/validation/stress_hosts.txt --include genome --filename h.zip --no-progressbar
unzip -oq h.zip -d hx; mkdir -p hosts; for f in hx/ncbi_dataset/data/GC*/*.fna; do mv "$f" hosts/$(basename $(dirname "$f")).fna; done
python3 $GMG_REPO/scripts/stress_sim.py hosts $GMG_REPO/data/validation/stress_hosts.meta.tsv sim
cp sim/genomes/*.fa.gz $C/s3_files/
for f in $C/s3_files/*.fa.gz; do ln -sf $f $Q/s3/$(basename $f); done
cp sim/truth.tsv $C/sets/s3_truth.tsv
echo "s1 $(ls $Q/s1 | wc -l) s2 $(ls $Q/s2 | wc -l) s3 $(ls $Q/s3 | wc -l)"
