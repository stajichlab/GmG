#!/bin/bash -l
#SBATCH -p short -N1 -c4 --mem=8G -t 2:00:00 -J f1bridge -o %x_%A_%a.out
#SBATCH --array=0-39
# F1 step 3: blast the FusariumID loci against each fusarioid genome and call its FusariumID species (genome_markers.py).
# usage: REPO=<checkout> W=<build dir> sbatch run_bridge.sbatch     (after run_genomes.sbatch wrote $W/genomes/paths.txt)
#        LIMIT=2 sbatch --array=0 ... for a timing test
# Genomes are decompressed into $SCRATCH; only $W/bridge/calls_<task>.tsv is kept.
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}
module load miniconda3 ncbi-blast
F1="$REPO/bench/F1_fusariumid"
N=${SLURM_ARRAY_TASK_COUNT:-40}; I=${SLURM_ARRAY_TASK_ID:-0}; LIMIT=${LIMIT:-0}
mkdir -p "$W/bridge"
LOCI="tef1 rpb1 rpb2 CaM tub2"
Q=$SCRATCH/query.fa; : > "$Q"
for l in $LOCI; do cat "$W/fid/loci/$l.fa" >> "$Q"; done
OUT="$W/bridge/calls_$I.tsv"
python3 -I "$F1/genome_markers.py" --header x /dev/null "$W/fid/records.tsv" | head -1 > "$OUT"
n=0
mapfile -t FILES < <(awk -v n=$N -v i=$I '(NR-1) % n == i' "$W/genomes/paths.txt")
for f in "${FILES[@]}"; do
  id=$(basename "$f" .fa.gz)
  t0=$(date +%s)
  zcat "$f" > "$SCRATCH/g.fa"
  makeblastdb -in "$SCRATCH/g.fa" -dbtype nucl -out "$SCRATCH/g" > /dev/null
  blastn -task megablast -query "$Q" -db "$SCRATCH/g" -evalue 1e-20 -max_target_seqs 5 -max_hsps 1 -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send bitscore" > "$SCRATCH/hits.tsv"
  python3 -I "$F1/genome_markers.py" "$id" "$SCRATCH/hits.tsv" "$W/fid/records.tsv" >> "$OUT"
  echo "$id $(( $(date +%s) - t0 ))s"
  rm -f "$SCRATCH"/g.*
  n=$((n+1))
  if [ "$LIMIT" -gt 0 ] && [ $n -ge "$LIMIT" ]; then break; fi
done
