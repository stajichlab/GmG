#!/bin/bash -l
#SBATCH -J rc_run
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --array=0-29%8
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/32_readscal.%A_%a.log
# #32: abundance cutoff for reads, fixed 1 / fixed 2 / auto (valley of the histogram), at several depths on three real read collections:
# ZyGoLife paired Illumina (titration sheet, genus truth), Coccidioides (species truth) and BFD RNA-seq (species truth).
# One array task = one (collection, depth, cutoff). The truth tables get a deterministic cal/eval split (scripts/reads_calibration_split.py).
# Uses a private binary (GMGBIN) so the shared $GR/bin/gmg symlink is not touched.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
GMGBIN=${GMGBIN:-$GR/bin/gmg-readscal}
REPO=${REPO:-$GR/reads_cal_repo}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
R=$GR/reads/calibration; mkdir -p $R

# collection | depths (reads per file) ; cutoffs
DS=(zygo cocci rna)
DEPTHS_zygo=(15000 75000 300000 1500000)
DEPTHS_cocci=(100000 500000 2000000)
DEPTHS_rna=(100000 500000 2000000)
CUT=(1 2 auto)
# enumerate the 30 combinations in a fixed order
i=0; PLAN=()
for ds in "${DS[@]}"; do
  eval "dl=(\${DEPTHS_${ds}[@]})"
  for n in "${dl[@]}"; do for c in "${CUT[@]}"; do PLAN+=("$ds $n $c"); done; done
done
read -r ds N ma <<< "${PLAN[$SLURM_ARRAY_TASK_ID]}"
echo "task $SLURM_ARRAY_TASK_ID: $ds reads=$N min-abund=$ma  binary $GMGBIN  index $INDEX"

case $ds in
  zygo)  INPUT="--sample-sheet $GR/reads/titration/sheet.tsv"; TRUTH=$GR/reads/titration/truth.tsv; BY="split genus_in_ref" ;;
  cocci) INPUT="$GR/reads/cocci_validation/farm"; TRUTH=$REPO/data/validation/cocci_reads.truth.tsv; BY="split" ;;
  rna)   INPUT="$GR/reads/rna_validation/farm"; TRUTH=$GR/reads/rna_validation/rna.truth.tsv; BY="split" ;;
esac
[ -s $R/$ds.truth.split.tsv ] || python3 $REPO/scripts/reads_calibration_split.py $TRUTH $R/$ds.truth.split.tsv
O=$R/$ds/n${N}_a${ma}; mkdir -p $R/$ds
if [ ! -s $O/results.tsv ]; then
  /usr/bin/time -f "wall=%e s maxrss=%M KB" $GMGBIN --threads ${SLURM_CPUS_PER_TASK} run $INPUT --index $INDEX -o $O --no-report -q --max-reads $N --min-abund $ma
fi
$GMGBIN validate $O/results.tsv $R/$ds.truth.split.tsv -o $O/v --by $BY --title "$ds, first $N reads, min-abund $ma"
