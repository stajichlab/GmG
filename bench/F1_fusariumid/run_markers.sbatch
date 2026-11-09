#!/bin/bash -l
#SBATCH -p short -N1 -c16 --mem=48G -t 2:00:00 -J f1markers -o %x_%j.out
# F1 step 1: marker-only resolution of fusarioid species in the FusariumID database.
# usage: REPO=<checkout> W=<build dir> FAS=<fusarioidID*.fas> sbatch run_markers.sbatch
# Writes $W/fid (parsed records), $W/blast (all-vs-all per locus) and $W/markers (tables); results are copied to the repo by hand.
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}; FAS=${FAS:?set FAS}
module load miniconda3 ncbi-blast
F1="$REPO/bench/F1_fusariumid"
python3 -I "$F1/parse_fusariumid.py" "$FAS" "$W/fid"
mkdir -p "$W/blast" "$W/db"
for fa in "$W"/fid/loci/*.fa; do
  l=$(basename "$fa" .fa)
  makeblastdb -in "$fa" -dbtype nucl -out "$W/db/$l" > /dev/null
  blastn -task blastn -query "$fa" -db "$W/db/$l" -dust no -evalue 1e-10 -max_target_seqs 500 -max_hsps 1 \
    -num_threads $SLURM_CPUS_PER_TASK -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen bitscore" \
    > "$W/blast/$l.tsv"
  echo "$l $(wc -l < "$W/blast/$l.tsv") hits"
done
python3 -I "$F1/marker_resolution.py" "$W/fid" "$W/blast" "$W/markers"
