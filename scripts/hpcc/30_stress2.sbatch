#!/bin/bash -l
#SBATCH -J gmg_stress2
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/30_stress2.%j.log
# Second stress test (#41, #45): chimeras and hybrids from 24 genome pairs, mitochondrion/plasmid/virus-only inputs, archaea and more
# non-fungal negatives. Per-input output stays here; scripts/stress2_score.py writes the aggregates ($R/tables/*.md, $R/report/).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
GMG_REPO=${STRESS2_REPO:-$GMG_REPO}   # a separate clone of the branch, so the shared checkout is not switched
module load ncbi_datasets/18.30.1
# newest release index (2026.10.3 = ambiguous-label genomes left out + distant-decoy fix), else 2026.10.2
INDEX=${INDEX:-$(ls -t $GR/release/2026.10.3/bfd-2026.10.3-s10000.gmgi $GR/release/2026.10.2/bfd-2026.10.2-s10000.gmgi 2>/dev/null | head -1)}
echo "index: $INDEX"; echo "gmg: $(readlink -f $GMG)"
R=$GR/stress2; W=$R/work; mkdir -p $W/genomes $W/small $W/neg $R/tables $R/report; cd $W
P=$GMG_REPO/data/validation/stress2_pairs.tsv

# 1. chimeras and hybrids
tail -n +2 $P | cut -f2,3 | tr '\t' '\n' | sort -u > pair_acc.txt
datasets download genome accession --inputfile pair_acc.txt --include genome --filename p.zip --no-progressbar || exit 1
unzip -oq p.zip -d px
for f in px/ncbi_dataset/data/GC*/*.fna; do mv "$f" genomes/$(basename $(dirname "$f")).fna; done
python3 $GMG_REPO/scripts/stress2_sim.py $P genomes sim || exit 1
$GMG --threads ${SLURM_CPUS_PER_TASK} run sim/genomes --index $INDEX -o $R/out_chim --no-report --top 10
cp sim/truth.tsv $R/truth_chim.tsv

# 2. mitochondrion / plasmid / virus only
python3 $GMG_REPO/scripts/fetch_small_inputs.py small
$GMG --threads ${SLURM_CPUS_PER_TASK} run small --index $INDEX -o $R/out_small --no-report
cp small/small.tsv $R/small.tsv

# 3. archaea and more non-fungal negatives
grep -v '^#' $GMG_REPO/data/validation/negatives_2.txt | cut -d' ' -f1 > neg_acc.txt
datasets download genome accession --inputfile neg_acc.txt --include genome --filename n.zip --no-progressbar
unzip -oq n.zip -d nx
for f in nx/ncbi_dataset/data/GC*/*.fna; do mv "$f" neg/$(basename $(dirname "$f")).fna; done
ls neg | wc -l
$GMG --threads ${SLURM_CPUS_PER_TASK} run neg --index $INDEX -o $R/out_neg --no-report

# 4. aggregates
python3 $GMG_REPO/scripts/stress2_score.py chim $R/out_chim/results.tsv $R/out_chim/hits.tsv $R/truth_chim.tsv $R/tables/chim.md
python3 $GMG_REPO/scripts/stress2_score.py small $R/out_small/results.tsv $R/small.tsv $R/tables/small.md
python3 $GMG_REPO/scripts/stress2_score.py neg $R/out_neg/results.tsv $GMG_REPO/data/validation/negatives_2.txt $R/tables/neg.md
$GMG validate $R/out_chim/results.tsv $R/truth_chim.tsv -o $R/report --by kind variant pct --title "Chimeras and hybrids (#41, #45)"
echo "$INDEX" > $R/tables/index_used.txt
rm -rf $W
