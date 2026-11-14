#!/usr/bin/env python3
"""Join the BioSample attributes that PopSeq has harvested (its biosample cache shards) to the candidate tables.

usage: popset_attrs.py WORKDIR SHARD_DIR OUT.csv

Reads every WORKDIR/cand_*.csv (output of `popseq runs`), collects their sample_accession values, and writes one row per
BioSample found in SHARD_DIR/batch_*.csv with the columns popset_select.py reads (isolation_source, env_*, geo_loc_name, ...).
Nothing is queried from an archive: this reads PopSeq's own cache. (Proposed to PopSeq: expose these columns in `popseq runs`.)
"""
import csv
import glob
import os
import sys

KEEP = ["isolation_source", "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "lat_lon", "host", "isolate", "strain",
        "ploidy", "mating_type", "culture_collection", "propagation", "sample_name", "common name", "collection_date"]


def main(work, shards, out):
    need = set()
    for f in glob.glob(os.path.join(work, "cand_*.csv")):
        with open(f, newline="") as fh:
            need.update(r["sample_accession"] for r in csv.DictReader(fh))
    found = {}
    for f in sorted(glob.glob(os.path.join(shards, "batch_*.csv"))):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                b = r.get("biosample")
                if b in need:
                    found[b] = {k: (r.get(k) or "").strip() for k in KEEP}
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["biosample"] + KEEP)
        for b, d in found.items():
            w.writerow([b] + [d[k] for k in KEEP])
    print(f"{len(need)} biosamples wanted, {len(found)} found in the cache", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        raise SystemExit(2)
    main(*sys.argv[1:])
