#!/usr/bin/env python3
"""Pack a folder of sketches made by `popsketch sketch` into one .tar.gz that the popsketch web page can open
(popsketch.html?data=URL, or "Open a sketch bundle").

usage: popsketch_bundle.py SKETCH_DIR -o BUNDLE.tar.gz [--labels labels.tsv] [--note TEXT]

The bundle holds popsketch-bundle.json (format 1: k, scaled, the sample ids), the optional labels.tsv (first column the sample id), and for
each sample <id>.json and <id>.gsk, and <id>.half_a.gsk and <id>.half_b.gsk for the --floor-samples (default 20) samples that popsketch
would pick for the noise floor (the other halves would not be used, and are two thirds of the bytes). Every sketch in the folder must have been made
with the parameters in the folder's params.txt; a sample with a missing file is refused rather than skipped, so a bundle is complete.
"""
import argparse, io, json, os, re, sys, tarfile, time

FORMAT = 1


def fnv1a(s):
    h = 0xCBF29CE484222325
    for b in s.encode():
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def pick_for_halves(ids, n):
    """The same pseudo-random subset popsketch uses (floor.rs pick_for_halves): sorted by (fnv1a(id), id), first n."""
    return set(sorted(ids, key=lambda i: (fnv1a(i), i))[:n])


def parse_params(text):
    p = dict(kv.split("=", 1) for kv in text.strip().split(";"))
    return {"k": int(p["k"]), "scaled": int(p["scaled"]), "max_bases": int(p["max_bases"])}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sketch_dir")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--labels", help="label table to include (first column: sample id)")
    ap.add_argument("--note", default="")
    ap.add_argument("--floor-samples", type=int, default=20,
                    help="samples whose split halves are packed (default 20, popsketch's own default): the noise floor uses only these, "
                         "and the halves are about two thirds of a sketch's size")
    a = ap.parse_args(argv)
    d = a.sketch_dir
    try:
        params = parse_params(open(os.path.join(d, "params.txt")).read())
    except (OSError, KeyError, ValueError) as e:
        sys.exit(f"{d}: not a folder of sketches (params.txt missing or unreadable: {e})")
    ids = sorted(f[:-5] for f in os.listdir(d) if f.endswith(".json"))
    if not ids:
        sys.exit(f"{d}: no sketches (<id>.json)")
    bad = [i for i in ids if not re.fullmatch(r"[A-Za-z0-9_.-]+", i)]
    if bad:
        sys.exit(f"sample ids with characters a bundle cannot carry safely: {bad[:5]}")
    for i in ids:
        if not os.path.exists(os.path.join(d, i + ".gsk")):
            sys.exit(f"{i}: {i}.json without {i}.gsk (an unfinished sketch?)")
    picked = pick_for_halves(ids, a.floor_samples)
    halves = [i for i in ids if i in picked and os.path.exists(os.path.join(d, i + ".half_a.gsk")) and os.path.exists(os.path.join(d, i + ".half_b.gsk"))]
    if len(halves) < min(a.floor_samples, len(ids)):
        print(f"note: only {len(halves)} of the {min(a.floor_samples, len(ids))} samples that should carry halves have them; the noise floor falls back or uses fewer", file=sys.stderr)
    labels = None
    if a.labels:
        labels = open(a.labels, "rb").read()
        have = {l.split("\t")[0] for l in labels.decode().splitlines()[1:]}
        missing = [i for i in ids if i not in have]
        if missing:
            print(f"note: {len(missing)} of {len(ids)} samples have no row in the labels (first: {missing[:3]})", file=sys.stderr)
    manifest = {"format": FORMAT, "tool": "popsketch", **params, "samples": ids, "with_halves": len(halves),
                "labels": labels is not None, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "note": a.note}

    def add(tf, name, data):
        ti = tarfile.TarInfo(name)
        ti.size, ti.mtime, ti.mode = len(data), 0, 0o644  # fixed times: the same sketches pack to the same bytes
        tf.addfile(ti, io.BytesIO(data))

    with tarfile.open(a.out, "w:gz", compresslevel=6) as tf:
        add(tf, "popsketch-bundle.json", (json.dumps(manifest, indent=1) + "\n").encode())
        if labels is not None:
            add(tf, "labels.tsv", labels)
        for i in ids:
            for suffix in (".json", ".gsk", ".half_a.gsk", ".half_b.gsk"):
                if suffix.startswith(".half") and i not in halves:
                    continue
                p = os.path.join(d, i + suffix)
                if os.path.exists(p):
                    add(tf, i + suffix, open(p, "rb").read())
    print(f"{a.out}: {len(ids)} samples ({len(halves)} with split halves), k={params['k']} scaled={params['scaled']}, {os.path.getsize(a.out) / 1e6:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
