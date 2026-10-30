#!/bin/bash -l
#SBATCH -J gmg_indigo
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/19_indigo.%j.log
# External set (#40): ~200 endophyte/environmental ascomycete assemblies (IndigoAg_Public), independent of BFD.
# Labels (provider names from the submission sheets) are in $GR/external/indigo_labels.tsv (sample, name, n50, total_len).
# Read-only on the source directory; per-genome output stays on the HPCC, only aggregates go in the repo.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
SRC=/bigdata/stajichlab/shared/projects/IndigoAg_Public/genomes
R=$GR/external/indigo; mkdir -p $R/in; cd $R
for d in $SRC/STR*; do
  s=$(basename $d); f=$(ls $d/third_annotation/*.scaffolds.fa 2>/dev/null | head -1)
  [ -n "$f" ] && ln -sf "$f" in/$s.fa
done
ls in | wc -l
INDEX=${INDEX:-$GR/index/bfd_s10000.v1.gmgi}
$GMG --threads ${SLURM_CPUS_PER_TASK} run in --index $INDEX -o out --report-format html
cp out/results.tsv results.tsv
