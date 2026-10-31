#!/usr/bin/env python3
"""Measure the error rate of long reads by aligning them to the genome they come from (mappy, the minimap2 binding; `pip install mappy`).

usage: long_reads_measure.py LABEL READS.fastq[.gz] GENOME.fa[.gz] [N=400] [--header]
Prints one TSV line: label, n_reads, n_aligned, median_len, mean_len, error_pct, mismatch_pct, indel_pct (of the aligned bases, primary hits with
mapq >= 20 that cover >= 50% of the read; a gz file cut in the middle is fine: the complete records read are used). error_pct = (NM / alignment length):
a lower bound of the true error of the read (soft-clipped ends and unmappable reads are not counted). Used on the first 4 MB of each candidate run (pilot),
on the first 200 Mb of each chosen run and on the simulated reads (to check the simulator hits its target).
"""
import sys, gzip, zlib, statistics

import mappy


def fastq_records(path, n):
    """first n complete FASTQ records of a (possibly truncated) plain or gz file"""
    if path.endswith(".gz"):
        d = zlib.decompressobj(31)
        buf = b""
        with open(path, "rb") as f:
            while True:
                chunk = f.read(1 << 20)
                if not chunk: break
                try:
                    buf += d.decompress(chunk)
                except zlib.error:
                    break
                if buf.count(b"\n") > 4 * (n + 2) * 1 and len(buf) > 40_000_000: break
        text = buf.decode("ascii", "replace")
    else:
        with open(path) as f:
            text = f.read(200_000_000)
    lines = text.split("\n")
    out = []
    for i in range(0, len(lines) - 3, 4):
        if lines[i].startswith("@") and lines[i + 2].startswith("+") and len(lines[i + 1]) == len(lines[i + 3]):
            out.append(lines[i + 1])
        if len(out) >= 4 * n: break
    return out


def main():
    a = [x for x in sys.argv[1:] if x != "--header"]
    if "--header" in sys.argv:
        print("label\tn_reads\tn_aligned\tmedian_len\tmean_len\terror_pct\tmismatch_pct\tindel_pct")
    label, reads, genome = a[:3]
    n = int(a[3]) if len(a) > 3 else 400
    recs = fastq_records(reads, n)
    step = max(1, len(recs) // n)
    recs = recs[::step][:n]            # spread over the sample, not only the first reads
    lens = [len(s) for s in recs]
    al = mappy.Aligner(genome, preset="map-ont")
    nm = blen = mism = nal = 0
    for s in recs:
        best = None
        for h in al.map(s):
            if h.is_primary and h.mapq >= 20 and (h.q_en - h.q_st) >= 0.5 * len(s):
                best = h; break
        if best is None: continue
        nal += 1; nm += best.NM; blen += best.blen
        # mismatches = NM - indels: from the CIGAR
        ind = sum(l for l, op in best.cigar if op in (1, 2))
        mism += best.NM - ind
    if not blen:
        print(f"{label}\t{len(recs)}\t0\t{int(statistics.median(lens)) if lens else 0}\t\tNA\tNA\tNA"); return
    print(f"{label}\t{len(recs)}\t{nal}\t{int(statistics.median(lens))}\t{int(statistics.mean(lens))}\t{100*nm/blen:.2f}\t{100*mism/blen:.2f}\t{100*(nm-mism)/blen:.2f}")


if __name__ == "__main__":
    main()
