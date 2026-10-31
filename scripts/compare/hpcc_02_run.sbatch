#!/bin/bash -l
#SBATCH -J cmp_run
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --array=0-4
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/cmp_run.%A_%a.log
# Tool comparison (#46): run one tool on the query sets s1 s2 s3 (queries in $C/q/SET/*.fa.gz), wall time and peak memory per run
# in $C/out/TOOL/SET/run.time (GNU time -v). array: 0 skani search, 1 mash dist, 2 mash screen, 3 sourmash gather, 4 gmg run (reps index)
# usage: GMG_REPO=$GR/compare/repo sbatch hpcc_02_run.sbatch
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
C=$GR/compare; E=$GR/envs/skani/bin; T=${SLURM_ARRAY_TASK_ID}
NAMES=(skani mash mashscreen sourmash gmg); TOOL=${NAMES[$T]}
CPU=${SLURM_CPUS_PER_TASK:-8}
for S in s1 s2 s3; do
  O=$C/out/$TOOL/$S; mkdir -p $O; cd $O
  ls $C/q/$S/*.fa.gz > qlist.txt; N=$(wc -l < qlist.txt)
  echo "== $TOOL $S ($N queries)"
  case $TOOL in
  skani)
    /usr/bin/time -v -o run.time $E/skani search -d $C/ref/skani_db --ql qlist.txt -t $CPU -o hits.tsv 2> run.err ;;
  mash)
    module load mash/2.3
    /usr/bin/time -v -o run.time bash -c "mash sketch -p $CPU -s 10000 -k 21 -l qlist.txt -o qset >/dev/null 2>&1 && mash dist -p $CPU -d 0.25 $C/ref/mash_refs.msh qset.msh > hits.tsv" ;;
  mashscreen)
    module load mash/2.3
    scr() { f=$1; id=$(basename $f .fa.gz); mash screen -p 2 -w -i 0.75 $C/ref/mash_refs.msh $f > scr.$id.tsv 2>/dev/null; }
    export -f scr; export C
    /usr/bin/time -v -o run.time bash -c "cat qlist.txt | xargs -P $((CPU/2)) -I{} bash -c 'scr {}'"
    for f in scr.*.tsv; do id=${f#scr.}; id=${id%.tsv}; awk -v q=$id 'BEGIN{OFS="\t"}{print q,$0}' $f; done > hits.tsv; rm -f scr.*.tsv ;;
  sourmash)
    module load sourmash/4.9
    mkdir -p sig gat
    sk() { f=$1; id=$(basename $f .fa.gz); sourmash sketch dna -p k=21,scaled=10000 --name $id -o sig/$id.sig $f >/dev/null 2>&1; }
    ga() { f=$1; id=$(basename $f .fa.gz); sourmash gather sig/$id.sig $C/ref/sm_refs.zip --threshold-bp 100000 -k 21 --dna -o gat/$id.csv >/dev/null 2>&1; }
    export -f sk ga; export C
    /usr/bin/time -v -o run.time bash -c "cat qlist.txt | xargs -P $CPU -I{} bash -c 'sk {}' && cat qlist.txt | xargs -P $CPU -I{} bash -c 'ga {}'" ;;
  gmg)
    IDX=$GR/release/2026.10.3/bfd-reps-2026.10.3-s20000.gmgi; EXA=""
    if [ ! -s $IDX ]; then IDX=$GR/release/2026.10.2/bfd-reps-2026.10.2-s20000.gmgi; EXA="--exclude-ids $GMG_REPO/reports/label_audit/reference_exclusions.txt"; fi
    echo "index $IDX $EXA" > index.txt
    if [ $S = s2 ]; then
      for M in loo los log; do
        python3 $GMG_REPO/scripts/compare/make_sheet.py $C/sets/s2_truth.tsv $C/q/s2 $M > sheet.$M.tsv
        L=self; [ $M = los ] && L=species; [ $M = log ] && L=genus
        /usr/bin/time -v -o run.$M.time $GMG --threads $CPU run --sample-sheet sheet.$M.tsv --index $IDX $EXA --exclude-level $L -o gmg.$M --no-report -q 2> run.$M.err
      done
    else
      /usr/bin/time -v -o run.time $GMG --threads $CPU run $C/q/$S --index $IDX $EXA -o gmg --no-report -q 2> run.err
    fi ;;
  esac
  ls | head -30
done
