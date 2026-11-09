"""Fixture tests for bench/F1_fusariumid (FusariumID markers vs GmG sketches); no HPCC data needed."""
import importlib.util
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
F1 = ROOT / "bench" / "F1_fusariumid"


def load(name):
    spec = importlib.util.spec_from_file_location(name, F1 / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


H1 = ("1|CBS 131385|MN170445|Fusarium aberrans, <em>F. incarnatum-equiseti</em> species complex (FIESC), "
      "ex-type of <em>Fusarium aberrans</em>, Australia|tef1 (n1)")


def test_parse_header_and_species_names():
    p = load("parse_fusariumid")
    r = p.parse_header(H1)
    assert r["species"] == "Fusarium aberrans" and r["genus"] == "Fusarium"
    assert r["complex"] == "F. incarnatum-equiseti species complex (FIESC)"
    assert r["locus"] == "tef1" and r["ex_type"] == 1 and r["strain"] == "CBS 131385"
    n = p.normalize_species
    assert n("Fusarium callistephi (now synonym of Fusarium inflexum)") == "Fusarium inflexum"
    assert n("Fusarium clavus (as 'clavum')") == "Fusarium clavus"
    assert n("Fusarium caucasicum (status unclear)") == "Fusarium caucasicum"
    assert n("Fusarium oxysporum f.sp. dactyliferae (as f.sp. 'dactyliferum')") == "Fusarium oxysporum"
    # informal taxa stay distinct
    assert n("Fusarium sp. FTSC 18 (Laraba et al. 2022)") != n("Fusarium sp. FTSC 19 (Laraba et al. 2022)")


def test_parse_main_writes_loci(tmp_path):
    p = load("parse_fusariumid")
    fa = tmp_path / "in.fas"
    fa.write_text(f">{H1}\nacgt-ACGTN\n>2|X 1|AB1|Fusarium beta, FX|rpb2 (n3)\nGGGG\n")
    p.main(str(fa), str(tmp_path / "out"))
    rec = (tmp_path / "out" / "records.tsv").read_text().splitlines()
    assert len(rec) == 3 and rec[1].split("\t")[8] == "9"  # gap removed, length 9
    assert (tmp_path / "out" / "loci" / "tef1.fa").read_text() == ">1\nACGTACGTN\n"


def recs():
    # strain s1/s2: species A; s3: species B; s4: species C (singleton)
    rows = []
    for sid, strain, sp, locus in [
        ("1", "s1", "A", "tef1"), ("2", "s2", "A", "tef1"), ("3", "s3", "B", "tef1"), ("4", "s4", "C", "tef1"),
        ("5", "s1", "A", "rpb2"), ("6", "s2", "A", "rpb2"), ("7", "s3", "B", "rpb2"), ("8", "s4", "C", "rpb2"),
    ]:
        rows.append(dict(seq_id=sid, strain=strain, species=sp, complex="X", locus=locus, strain_key=(strain, sp)))
    return {r["seq_id"]: r for r in rows}


def test_marker_loo_and_multilocus():
    m = load("marker_resolution")
    rec = recs()
    # tef1: 1-2 differ by 1/500, 1-3 by 2/500, 3-4 by 1/500 ; rpb2: 5-6 by 0, 5-7 by 0 (tie A/B), 7-8 by 10
    def sym(pairs):
        h = {}
        for a, b, d in pairs:
            h.setdefault(a, {})[b] = (d, 500)
            h.setdefault(b, {})[a] = (d, 500)
        return h
    tef = sym([("1", "2", 1), ("1", "3", 2), ("2", "3", 3), ("3", "4", 1), ("1", "4", 5), ("2", "4", 6)])
    rpb = sym([("5", "6", 0), ("5", "7", 0), ("6", "7", 1), ("7", "8", 10), ("5", "8", 10), ("6", "8", 10)])
    rows = {r["query"]: r for r in m.per_locus(rec, tef, "tef1")}
    assert rows["1"]["outcome"] == "correct" and rows["1"]["answerable"]
    assert rows["3"]["answerable"] is False  # B has one strain
    assert rows["1"]["het_diff"] == 2
    rows = {r["query"]: r for r in m.per_locus(rec, rpb, "rpb2")}
    assert rows["5"]["outcome"] == "ambiguous"
    sp = m.strain_pairs(rec, {"tef1": tef, "rpb2": rpb})
    ml, near = m.multilocus(rec, sp, ["tef1", "rpb2"])
    by = {r["query"]: r for r in ml}
    assert by["s1"]["outcome"] == "correct"  # tef1 breaks the rpb2 tie: 1/1000 vs 2/1000
    assert near["A"] == 2
    s = m.summarise(ml, dict(set="x"))
    assert s["answerable"] == 2 and s["correct"] == 2


def test_kmer_hashes_canonical():
    m = load("marker_resolution")
    seq = "ACGTTGCAAGGCTTACCGATCGATCGGA"
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert np.array_equal(m.kmer_hashes(seq, 21), m.kmer_hashes(rc, 21))
    assert len(m.kmer_hashes("ACGTN" * 10, 21)) == 0


def test_genome_markers_call(tmp_path):
    g = load("genome_markers")
    rec = {k: v for k, v in recs().items()}
    for r in rec.values():
        r.pop("strain_key")
    blast = tmp_path / "hits.tsv"
    # qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send bitscore
    lines = [
        ("1", "c1", 99.8, 500, 1, 0, 1, 500, 500, 1, 500, 900),
        ("2", "c1", 99.6, 500, 2, 0, 1, 500, 500, 1, 500, 890),
        ("3", "c1", 99.0, 500, 5, 0, 1, 500, 500, 1, 500, 850),
        ("5", "c2", 100.0, 500, 0, 0, 1, 500, 500, 1, 500, 920),
        ("7", "c2", 100.0, 500, 0, 0, 1, 500, 500, 1, 500, 920),
        ("8", "c2", 99.0, 100, 1, 0, 1, 100, 500, 1, 100, 150),  # covers 20%: dropped
    ]
    blast.write_text("".join("\t".join(map(str, x)) + "\n" for x in lines))
    row = g.call("G1", str(blast), rec, ["tef1", "rpb2"])
    assert row["tef1_species"] == "A" and row["rpb2_species"] == "A;B"
    assert row["combined_species"] == "A" and row["combined_diffs"] == "1"
    assert row["combined_other_species"] == "B" and row["n_loci_hit"] == 2


def write_gsk(path, name, hashes, k=21, scaled=1000):
    with open(path, "wb") as f:
        f.write(b"GMGSKT01" + struct.pack("<II", k, scaled) + struct.pack("<Q", 42))
        f.write(struct.pack("<I", len(name)) + name.encode() + struct.pack("<Q", len(hashes)))
        f.write(np.array(sorted(hashes), dtype="<u8").tobytes() + b"\x00")


def test_genome_ani(tmp_path):
    a = load("genome_ani")
    write_gsk(tmp_path / "x.gsk", "x", [1, 2, 3, 4])
    name, k, sc, h = a.read_gsk(tmp_path / "x.gsk")
    assert (name, k, sc) == ("x", 21, 1000) and list(h) == [1, 2, 3, 4]
    assert a.count_shared(np.array([1, 3, 5], dtype=np.uint64), np.array([2, 3, 5, 9], dtype=np.uint64)) == 2
    big = np.uint64(2**64 - 1)
    assert list(a.downsample(np.array([5, big], dtype=np.uint64), 10)) == [5]


def test_genome_resolution_helpers():
    r = load("genome_resolution")
    assert r.agreement("Fusarium oxysporum", "Fusarium oxysporum") == "same"
    assert r.agreement("Fusarium oxysporum", "Fusarium nirenbergiae") == "differs"
    assert r.agreement("Fusarium cf. solani", "Neocosmospora solani") == "ncbi_unnamed"
    assert r.agreement("X y", "") == "no_fid_call"
    labels = {"a": "A", "b": "A", "c": "B", "d": "B", "e": "C"}
    pairs = [("a", "b", 0.99), ("a", "c", 0.98), ("c", "d", 0.995), ("b", "e", 0.97)]
    got = {g: (c, h) for g, _l, c, h in r.ani_gap(pairs, labels)}
    assert "e" not in got  # singleton label
    assert got["a"] == (0.99, 0.98) and got["b"] == (0.99, 0.97)
    calls = [dict(genome="g1", combined_species="A", combined_dist="0.002"),
             dict(genome="g2", combined_species="B", combined_dist="0.05")]
    assert r.fid_labels(calls, 0.01) == {"g1": "A"}
