#!/bin/bash -l
#SBATCH -J gmg_clskani
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/27_cluster_skani.%j.log
# #47: the 95% clusters of `gmg cluster` (release 2026.10.2 metadata = index/clusters_ani0.95.tsv) against skani.
# usage: sbatch 27_cluster_skani.sbatch [SCRIPT_DIR]   (SCRIPT_DIR holds cluster_vs_skani.py; default: the GmG checkout)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
SD=${1:-$GMG_REPO/scripts}
SK=${SKANI:-$HOME/projects/GmG/.pixi/envs/bio/bin/skani}
R=$GR/cluster_skani; mkdir -p $R/skani
$SK --version
python3 $SD/cluster_vs_skani.py sample $GR/index/clusters_ani0.95.tsv $BFD_GENOMES $R || exit 1
cd $R
ls lists/*.txt | xargs -P ${SLURM_CPUS_PER_TASK} -I{} bash -c 'n=$(basename {} .txt); [ -s skani/$n.tsv ] || '"$SK"' triangle -l {} -E -t 1 -o skani/$n.tsv 2> skani/$n.err'
python3 $SD/cluster_vs_skani.py score $GR/index/clusters_ani0.95.tsv $R $R/report.md > /dev/null
ls skani | wc -l
