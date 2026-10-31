#!/bin/bash -l
#SBATCH -J mixdet_contam
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_mixdet_contam.%j.log
# #31 negative set: contamination spike-in of #41 (bacteria or human DNA at 1-30% in a fungal host), classified with the mixture-detector binary. The 24 stress hosts + E. coli/B. subtilis or human chr20 DNA at 1-30% (as 20_contam.sbatch).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
GMG=${MG:-$GR/mixdet/bin/gmg}; INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
# HOSTS/META/NAME select the host set (default: the 24 stress hosts, species not in the reference; contam_hosts_species: species in the reference)
HOSTS=${HOSTS:-$GMG_REPO/data/validation/stress_hosts.txt}; META=${META:-$GMG_REPO/data/validation/stress_hosts.meta.tsv}
R=$GR/mixdet/out/${NAME:-contam}; W=$R/work; mkdir -p $W/hosts $W/donors; cd $W
datasets download genome accession --inputfile $HOSTS --include genome --filename h.zip --no-progressbar || exit 1
unzip -oq h.zip -d hx; for f in hx/ncbi_dataset/data/GC*/*.fna; do mv "$f" hosts/$(basename $(dirname "$f")).fna; done
printf 'GCF_000005845.2\nGCF_000009045.1\n' > bac.txt
datasets download genome accession --inputfile bac.txt --include genome --filename b.zip --no-progressbar || exit 1
unzip -oq b.zip -d bx; cat bx/ncbi_dataset/data/GC*/*.fna > donors/bacteria.fna
datasets download genome accession GCF_000001405.40 --chromosomes 20 --include genome --filename hs.zip --no-progressbar || exit 1
unzip -oq hs.zip -d hsx; cat hsx/ncbi_dataset/data/GC*/*.fna > donors/human.fna
python3 $GMG_REPO/scripts/contam_sim.py hosts $META donors sim || exit 1
$GMG --threads ${SLURM_CPUS_PER_TASK} run sim/genomes --index $INDEX -o $R/out --thresholds known --no-report
cp sim/truth.tsv $R/truth.tsv
rm -rf $W
