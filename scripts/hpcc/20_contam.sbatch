#!/bin/bash -l
#SBATCH -J gmg_contam
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/20_contam.%j.log
# Contamination spike-in (#41): the 24 stress hosts + bacterial (E. coli, B. subtilis) or human (chr20) DNA at 1-30%, classified with
# the decoy release index, scored against the host lineage with `gmg validate`. Per-input output stays here; aggregates go in the repo.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
INDEX=${INDEX:-$GR/release/2026.10.2/bfd-2026.10.2-s10000.gmgi}
R=$GR/contam; W=$R/work; mkdir -p $W/hosts $W/donors; cd $W
datasets download genome accession --inputfile $GMG_REPO/data/validation/stress_hosts.txt --include genome --filename h.zip --no-progressbar || exit 1
unzip -oq h.zip -d hx; for f in hx/ncbi_dataset/data/GC*/*.fna; do mv "$f" hosts/$(basename $(dirname "$f")).fna; done
printf 'GCF_000005845.2\nGCF_000009045.1\n' > bac.txt
datasets download genome accession --inputfile bac.txt --include genome --filename b.zip --no-progressbar || exit 1
unzip -oq b.zip -d bx; cat bx/ncbi_dataset/data/GC*/*.fna > donors/bacteria.fna
datasets download genome accession GCF_000001405.40 --chromosomes 20 --include genome --filename hs.zip --no-progressbar || exit 1
unzip -oq hs.zip -d hsx; cat hsx/ncbi_dataset/data/GC*/*.fna > donors/human.fna
python3 $GMG_REPO/scripts/contam_sim.py hosts $GMG_REPO/data/validation/stress_hosts.meta.tsv donors sim || exit 1
$GMG --threads ${SLURM_CPUS_PER_TASK} run sim/genomes --index $INDEX -o $R/out --no-report
cp sim/truth.tsv $R/truth.tsv
$GMG validate $R/out/results.tsv $R/truth.tsv -o $R/report --by donor frac --title "Contamination spike-in (#41)"
rm -rf $W
