#!/usr/bin/env python3
"""Stream the first N bases (default 200 Mb) of long-read ENA runs, never a whole run (#44).

usage: long_reads_harvest.py RUNS.tsv OUT_DIR [BASES=200e6] [WORKERS=3] [BUDGET_BYTES=10e9]
RUNS.tsv: data/validation/long_reads_runs.tsv (columns run, url = ENA path of the FASTQ, fastq_bytes, base_count). Each file is read over HTTPS
(ftp.sra.ebi.ac.uk) only until BASES bases of complete records have been decompressed, then the connection is closed; written as OUT_DIR/<run>.fastq.gz.
A dropped or early-ended connection resumes with an HTTP Range request at the next compressed byte and the same decompressor (no byte read twice); a run
that cannot reach BASES is `.failed` (files removed; the fallback is `fastq-dump -X`, scripts/hpcc/37_long_reads_real.sbatch). Compressed bytes read from
the network go to OUT_DIR/bytes.tsv (run, bytes, bases, reads, seconds, resumes). Before a run starts its expected cost (file size x BASES / run bases, +20%)
must fit in BUDGET_BYTES minus what has been read; otherwise the harvest stops. Read-only on ENA.
"""
import csv, gzip, http.client, os, sys, time, zlib, threading, urllib.request, concurrent.futures

runs_p, out = sys.argv[1:3]
BASES = int(float(sys.argv[3])) if len(sys.argv) > 3 else 200_000_000
WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else 3
BUDGET = float(sys.argv[5]) if len(sys.argv) > 5 else 10e9
os.makedirs(out, exist_ok=True)
lock = threading.Lock()
total = [0]
stop = [False]


def stream(path, dest, bases_wanted):
    url = path if "://" in path else "https://" + path
    got = bases = reads = resumes = stalled = 0
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    rest = b""                                   # an incomplete record carried to the next chunk
    with gzip.open(dest + ".part", "wb", compresslevel=1) as fo:
        while bases < bases_wanted:
            before = got
            hdr = {"User-Agent": "gmg-validation/1 (reads head only)"}
            if got: hdr["Range"] = f"bytes={got}-"
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=120) as resp:
                    if got and getattr(resp, "status", 206) == 200 and url.startswith("http"):
                        raise IOError("server ignored the Range request")
                    while bases < bases_wanted:
                        chunk = resp.read(1 << 18)
                        if not chunk: break
                        got += len(chunk)
                        data = dec.decompress(chunk)
                        while dec.eof and dec.unused_data:
                            r2 = dec.unused_data
                            dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
                            data += dec.decompress(r2)
                        if not data: continue
                        lines = (rest + data).split(b"\n")
                        n_complete = (len(lines) - 1) // 4 * 4        # the last element is the unfinished line
                        keep, tail = lines[:n_complete], lines[n_complete:]
                        rest = b"\n".join(tail)
                        buf = []
                        for i in range(0, n_complete, 4):
                            buf.append(b"\n".join(keep[i:i + 4]) + b"\n")
                            bases += len(keep[i + 1]); reads += 1
                            if bases >= bases_wanted: break
                        fo.write(b"".join(buf))
            except (OSError, zlib.error, http.client.HTTPException):
                if got == before:
                    stalled += 1
                    if stalled >= 5: raise
                    time.sleep(10 * stalled)
                else:
                    stalled = 0
                resumes += 1
                continue
            if bases < bases_wanted:
                if got == before:
                    stalled += 1
                    if stalled >= 3: raise IOError(f"stream ended after {bases} bases")
                else:
                    stalled = 0
                resumes += 1
    os.replace(dest + ".part", dest)
    return got, bases, reads, resumes


def one(r):
    run = r["run"]
    if os.path.exists(f"{out}/{run}.done") or os.path.exists(f"{out}/{run}.failed") or stop[0]: return
    expect = 1.2 * float(r["fastq_bytes"]) * min(1.0, BASES / float(r["base_count"]))
    with lock:
        if total[0] + expect > BUDGET:
            stop[0] = True; print("budget reached, stopping before", run, file=sys.stderr); return
        total[0] += expect                     # reserved; replaced by the real figure below
    t0 = time.time()
    try:
        got, bases, reads, resumes = stream(r["url"], f"{out}/{run}.fastq.gz", BASES)
    except Exception as e:
        print(run, "FAILED", repr(e)[:120], file=sys.stderr)
        for p in (f"{out}/{run}.fastq.gz", f"{out}/{run}.fastq.gz.part"):
            if os.path.exists(p): os.remove(p)
        open(f"{out}/{run}.failed", "w").close()
        with lock: total[0] -= expect
        return
    with lock:
        total[0] += got - expect
        with open(f"{out}/bytes.tsv", "a") as fh:
            fh.write(f"{run}\t{got}\t{bases}\t{reads}\t{time.time() - t0:.0f}\t{resumes}\n")
    open(f"{out}/{run}.done", "w").close()


rows = list(csv.DictReader(open(runs_p), delimiter="\t"))
if os.path.exists(f"{out}/bytes.tsv"):
    total[0] = sum(int(f[1]) for f in (l.split("\t") for l in open(f"{out}/bytes.tsv")) if f[1].isdigit())
else:
    open(f"{out}/bytes.tsv", "w").write("run\tbytes\tbases\treads\tseconds\tresumes\n")
with concurrent.futures.ThreadPoolExecutor(WORKERS) as ex:
    list(ex.map(one, rows))
done = sum(os.path.exists(os.path.join(out, r["run"] + ".done")) for r in rows)
failed = sum(os.path.exists(os.path.join(out, r["run"] + ".failed")) for r in rows)
print(f"transferred {total[0] / 1e9:.2f} GB (compressed bytes read) for {done} of {len(rows)} runs; {failed} failed", file=sys.stderr)
