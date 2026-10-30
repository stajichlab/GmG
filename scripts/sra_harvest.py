#!/usr/bin/env python3
"""Stream the first read pairs of ENA runs (never the whole run) for the real-reads validation (#42).

usage: sra_harvest.py RUNS.tsv OUT_DIR [PAIRS=500000] [WORKERS=3] [BUDGET_BYTES=24e9]
RUNS.tsv: the accession list of sra_plan.py (columns run, fastq_1, fastq_2: ENA FTP paths, fetched over HTTPS from ftp.sra.ebi.ac.uk). For each run and each
mate the connection is read only until PAIRS reads (4 x PAIRS lines) have been decompressed, then closed; the reads are written as OUT_DIR/<run>_1.fastq.gz and
_2.fastq.gz. If the server ends or drops a connection early, the download resumes with an HTTP Range request at the next compressed byte and the same decompressor
(no byte is read twice); a mate that still has fewer reads than asked after several resumes is an error, and a run is `.done` only when BOTH mates have PAIRS reads
(otherwise `.failed`, with the files removed). The compressed bytes actually read from the network, resumes included, go to OUT_DIR/bytes.tsv
(run, bytes_1, bytes_2, pairs, seconds, resumes); the total is checked against BUDGET_BYTES before every run and the harvest stops if it would be exceeded.
Read-only on ENA; nothing but these two streams per run is requested.
"""
import csv, gzip, http.client, os, sys, time, zlib, threading, urllib.request, concurrent.futures

runs_p, out = sys.argv[1:3]
PAIRS = int(sys.argv[3]) if len(sys.argv) > 3 else 500_000
WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else 3
BUDGET = float(sys.argv[5]) if len(sys.argv) > 5 else 24e9
os.makedirs(out, exist_ok=True)
lock = threading.Lock()
total = [0]
stop = [False]


def stream(path, dest, lines_wanted):
    """Decompress the gzip stream of one FASTQ until lines_wanted lines are written; returns (compressed bytes read, lines written, resumes)."""
    url = path if "://" in path else "https://" + path     # ENA lists FTP paths without a scheme; file:// URLs are for the tests
    got = lines = resumes = stalled = 0
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    with gzip.open(dest + ".part", "wb", compresslevel=1) as fo:
        while lines < lines_wanted:
            before = got
            hdr = {"User-Agent": "gmg-validation/1 (reads head only)"}
            if got: hdr["Range"] = f"bytes={got}-"
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=120) as resp:
                    if got and getattr(resp, "status", 206) == 200 and url.startswith("http"):
                        raise IOError("server ignored the Range request")
                    while lines < lines_wanted:
                        chunk = resp.read(1 << 18)
                        if not chunk: break
                        got += len(chunk)
                        data = dec.decompress(chunk)
                        while dec.eof and dec.unused_data:      # a multi-member gzip file
                            rest = dec.unused_data
                            dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
                            data += dec.decompress(rest)
                        if not data: continue
                        n = data.count(b"\n")
                        if lines + n > lines_wanted:
                            keep = lines_wanted - lines
                            pos = -1
                            for _ in range(keep): pos = data.index(b"\n", pos + 1)
                            data, n = data[:pos + 1], keep
                        fo.write(data); lines += n
            except (OSError, zlib.error, http.client.HTTPException) as e:
                if got == before:
                    stalled += 1
                    if stalled >= 5: raise
                    time.sleep(10 * stalled)
                else:
                    stalled = 0
                resumes += 1
                continue
            if lines < lines_wanted:       # the server ended the stream early: resume where it stopped
                if got == before:
                    stalled += 1
                    if stalled >= 3: raise IOError(f"stream ended after {lines // 4} reads")
                else:
                    stalled = 0
                resumes += 1
    os.replace(dest + ".part", dest)
    return got, lines, resumes


def one(r):
    run = r["run"]
    if os.path.exists(f"{out}/{run}.done") or os.path.exists(f"{out}/{run}.failed") or stop[0]: return
    with lock:
        if total[0] + 2 * 150e6 > BUDGET:      # a run costs about 70-150 MB; do not start one the budget cannot cover
            stop[0] = True; print("budget reached, stopping", file=sys.stderr); return
    t0 = time.time(); b = []
    for k in (1, 2):
        try:
            got, lines, resumes = stream(r[f"fastq_{k}"], f"{out}/{run}_{k}.fastq.gz", 4 * PAIRS)
        except Exception as e:
            print(run, k, "FAILED", repr(e)[:120], file=sys.stderr)
            for j in (1, 2):
                for p in (f"{out}/{run}_{j}.fastq.gz", f"{out}/{run}_{j}.fastq.gz.part"):
                    if os.path.exists(p): os.remove(p)
            open(f"{out}/{run}.failed", "w").close()
            return
        b.append((got, lines // 4, resumes))
        with lock: total[0] += got
    with lock:
        with open(f"{out}/bytes.tsv", "a") as fh:
            fh.write(f"{run}\t{b[0][0]}\t{b[1][0]}\t{min(b[0][1], b[1][1])}\t{time.time() - t0:.0f}\t{b[0][2] + b[1][2]}\n")
    open(f"{out}/{run}.done", "w").close()


rows = list(csv.DictReader(open(runs_p), delimiter="\t"))
if os.path.exists(f"{out}/bytes.tsv"):
    total[0] = sum(int(f[1]) + int(f[2]) for f in (l.split("\t") for l in open(f"{out}/bytes.tsv")) if f[1].isdigit())
else:
    open(f"{out}/bytes.tsv", "w").write("run\tbytes_1\tbytes_2\tpairs\tseconds\tresumes\n")
with concurrent.futures.ThreadPoolExecutor(WORKERS) as ex:
    list(ex.map(one, rows))
done = sum(os.path.exists(os.path.join(out, r["run"] + ".done")) for r in rows)
failed = sum(os.path.exists(os.path.join(out, r["run"] + ".failed")) for r in rows)
print(f"transferred {total[0] / 1e9:.2f} GB (compressed bytes read) for {done} of {len(rows)} runs; {failed} failed", file=sys.stderr)
