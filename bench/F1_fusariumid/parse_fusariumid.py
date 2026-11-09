#!/usr/bin/env python3
"""Parse the FusariumID (fusarium.org) sequence download into a record table and one FASTA per locus.

Input: the FASTA inside `sequencesindatabase<DATE>.zip` from https://www.fusarium.org/page/ID. Headers look like
    >1|CBS 131385|MN170445|Fusarium aberrans, <em>F. incarnatum-equiseti</em> species complex (FIESC), ..., ex-type of ...|tef1 (n1)
i.e. `seq_id | strain | accession | species, complex, [ex-type of X,] origin... | locus (nN)`.

Output (in OUTDIR):
    records.tsv       seq_id strain accession species genus complex ex_type locus length n_frac
    loci/<locus>.fa   sequences of that locus, ids = seq_id, upper case, gaps removed

usage: parse_fusariumid.py FASTA OUTDIR
"""
import html
import os
import re
import sys

TAG = re.compile(r"<[^>]+>")
AS_NAME = re.compile(r"\s*\(as '[^']*'\)")


def clean(s):
    return html.unescape(TAG.sub("", s)).strip()


def parse_header(h):
    """Return a dict of fields from one FusariumID header (without '>')."""
    parts = h.split("|")
    if len(parts) < 5:
        raise ValueError(f"unexpected header: {h[:80]}")
    seq_id, strain, acc = parts[0].strip(), parts[1].strip(), parts[2].strip()
    locus = parts[-1].split("(")[0].strip()
    desc = "|".join(parts[3:-1])
    fields = [clean(f) for f in desc.split(",")]
    species = AS_NAME.sub("", fields[0]).strip()
    species = " ".join(species.split()[:2]) if species.split()[:1] else species
    complex_ = fields[1] if len(fields) > 1 else ""
    ex_type = int("ex-type of" in desc)
    return {
        "seq_id": seq_id,
        "strain": strain,
        "accession": acc,
        "species": species,
        "genus": species.split()[0] if species else "",
        "complex": complex_,
        "ex_type": ex_type,
        "locus": locus,
    }


def read_fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main(fasta, outdir):
    os.makedirs(os.path.join(outdir, "loci"), exist_ok=True)
    handles = {}
    cols = ["seq_id", "strain", "accession", "species", "genus", "complex", "ex_type", "locus", "length", "n_frac"]
    n = 0
    with open(os.path.join(outdir, "records.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for h, s in read_fasta(fasta):
            r = parse_header(h)
            s = s.upper().replace("-", "").replace(".", "")
            if not s:
                continue
            r["length"] = len(s)
            r["n_frac"] = f"{sum(c not in 'ACGT' for c in s) / len(s):.4f}"
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
            loc = re.sub(r"[^A-Za-z0-9]+", "_", r["locus"])
            if loc not in handles:
                handles[loc] = open(os.path.join(outdir, "loci", f"{loc}.fa"), "w")
            handles[loc].write(f">{r['seq_id']}\n{s}\n")
            n += 1
    for fh in handles.values():
        fh.close()
    print(f"{n} records, {len(handles)} loci -> {outdir}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
