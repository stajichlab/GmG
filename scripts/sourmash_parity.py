#!/usr/bin/env python3
"""FracMinHash parity with Python sourmash (#49): the hash set of `gmg sketch` against `sourmash sketch dna` on real genomes.

usage: sourmash_parity.py GSK_DIR SIG_DIR OUT.tsv
  GSK_DIR: <id>.gsk written by `gmg sketch -k 21 -s 10000`; SIG_DIR: <id>.sig written by `sourmash sketch dna -p k=21,scaled=10000`
  (same k, scaled and the default seed 42). Compares the sorted hash sets; writes id, n_gmg, n_sourmash, shared, identical.
.gsk layout: "GMGSKT01" | k u32 | scaled u32 | seed u64 | name_len u32 | name | n u64 | hashes u64*n | has_abund u8 | ...
"""
import glob, gzip, json, os, struct, sys


def read_gsk(p):
    b = open(p, "rb").read()
    assert b[:8] == b"GMGSKT01", p
    k, scaled, seed, nl = struct.unpack_from("<IIQI", b, 8)
    off = 8 + 4 + 4 + 8 + 4 + nl
    (n,) = struct.unpack_from("<Q", b, off)
    off += 8
    return k, scaled, seed, list(struct.unpack_from(f"<{n}Q", b, off))


def read_sig(p):
    raw = open(p, "rb").read()
    try:
        raw = gzip.decompress(raw)
    except OSError:
        pass
    sigs = json.loads(raw)
    mh = sigs[0]["signatures"][0]
    return mh["ksize"], mh["max_hash"], mh["seed"], mh["mins"]


gsk, sig, out = sys.argv[1:4]
rows, bad = [], 0
for g in sorted(glob.glob(f"{gsk}/*.gsk")):
    i = os.path.basename(g)[:-4]
    sp = f"{sig}/{i}.sig"
    if not os.path.exists(sp):
        print("no sourmash signature for", i, file=sys.stderr)
        continue
    k, scaled, seed, h1 = read_gsk(g)
    k2, maxh, seed2, h2 = read_sig(sp)
    s1, s2 = set(h1), set(h2)
    same = s1 == s2 and (k, seed) == (k2, seed2)
    bad += not same
    rows.append((i, len(s1), len(s2), len(s1 & s2), "yes" if same else "NO", k, scaled, seed, seed2))
with open(out, "w") as f:
    f.write("id\tn_gmg\tn_sourmash\tshared\tidentical\tk\tscaled\tseed_gmg\tseed_sourmash\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print(f"{len(rows)} genomes compared, {len(rows) - bad} with identical hash sets, {bad} differ")
sys.exit(1 if bad else 0)
