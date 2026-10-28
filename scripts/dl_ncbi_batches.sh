#!/bin/bash
# usage: dl_ncbi_batches.sh ACCESSIONS.txt OUT_DIR   (needs `datasets`; run inside the popsketch conda env)
# Downloads genomic FASTA in batches of 20 accessions with retries; writes OUT_DIR/<accession>.fna
set -u
acc=$1; out=$2; mkdir -p $out/zips $out/tmp
split -l 20 -d -a 3 $acc $out/tmp/batch_
for b in $out/tmp/batch_*; do
  n=$(basename $b)
  [ -s $out/zips/$n.done ] && continue
  for try in 1 2 3 4; do
    if datasets download genome accession --inputfile $b --include genome --filename $out/zips/$n.zip --no-progressbar >/dev/null 2>$out/zips/$n.err; then
      rm -rf $out/zips/$n.x; unzip -oq $out/zips/$n.zip -d $out/zips/$n.x
      for f in $out/zips/$n.x/ncbi_dataset/data/GC*/*.fna; do cp $f $out/$(basename $(dirname $f)).fna; done
      rm -rf $out/zips/$n.x $out/zips/$n.zip; echo ok > $out/zips/$n.done; break
    fi
    sleep $((try*20))
  done
  [ -s $out/zips/$n.done ] || echo "FAILED batch $n" >> $out/failed.txt
done
echo finished > $out/dl.finished
