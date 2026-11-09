"""scripts/fusariumid_markers.py (GmG marker mode inputs and scoring) on tiny fixtures."""
import argparse
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load():
    spec = importlib.util.spec_from_file_location("fm", ROOT / "scripts" / "fusariumid_markers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_fixture(d):
    rows = [  # seq_id strain species locus
        ("1", "S1", "Fusarium alpha", "tef1"), ("2", "S1", "Fusarium alpha", "rpb2"),
        ("3", "S2", "Fusarium alpha", "tef1"), ("4", "S3", "Fusarium beta", "tef1"),
        ("5", "S4", "Fusarium gamma", "tef1"),
    ]
    cols = ["seq_id", "strain", "accession", "species", "genus", "complex", "ex_type", "locus", "length", "n_frac"]
    with open(d / "records.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for sid, st, sp, l in rows:
            f.write("\t".join([sid, st, "X", sp, "Fusarium", "FX", "0", l, "4", "0"]) + "\n")
    (d / "loci").mkdir()
    (d / "loci" / "tef1.fa").write_text(">1\nACGT\n>3\nACGA\n>4\nTTTT\n>5\nGGGG\n")
    (d / "loci" / "rpb2.fa").write_text(">2\nCCCC\n")


def test_prepare_all_refs(tmp_path):
    m = load()
    write_fixture(tmp_path)
    a = argparse.Namespace(records=str(tmp_path / "records.tsv"), loci_dir=str(tmp_path / "loci"), out=str(tmp_path / "o"),
                           holdout=0.0, seed=1, loci="tef1,rpb2")
    m.prepare(a)
    meta = (tmp_path / "o" / "meta.tsv").read_text().splitlines()
    assert len(meta) == 5  # header + 4 strains
    assert meta[1].split("\t")[2:8] == ["Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae", "Fusarium",
                                        "Fusarium alpha"]
    s1 = (tmp_path / "o" / "refs" / "fid1_S1.fa").read_text()
    assert ">1 tef1" in s1 and ">2 rpb2" in s1  # one reference per strain, all its loci
    assert not (tmp_path / "o" / "queries").exists()


def test_prepare_holdout_splits_by_strain(tmp_path):
    m = load()
    write_fixture(tmp_path)
    a = argparse.Namespace(records=str(tmp_path / "records.tsv"), loci_dir=str(tmp_path / "loci"), out=str(tmp_path / "o"),
                           holdout=0.5, seed=3, loci="tef1,rpb2")
    m.prepare(a)
    refs = {l.split("\t")[0] for l in (tmp_path / "o" / "meta.tsv").read_text().splitlines()[1:]}
    truth = (tmp_path / "o" / "queries_truth.tsv").read_text().splitlines()[1:]
    held = {l.split("\t")[3] for l in truth}
    assert len(refs) == 2 and len(held) == 2
    assert not any(r.endswith("_" + s) for r in refs for s in held)  # a strain is either reference or query


def test_score(tmp_path):
    m = load()
    (tmp_path / "truth.tsv").write_text("query_id\tlocus\tspecies\tstrain\tcomplex\n"
                                        "q1\ttef1\tFusarium alpha\tS\tFX\nq2\ttef1\tFusarium alpha\tS\tFX\n"
                                        "q3\ttef1\tFusarium zeta\tS\tFX\nq4\ttef1\tFusarium oxysporum\tS\tFOSC\n")
    (tmp_path / "meta.tsv").write_text("id\tspecies\nr1\tFusarium alpha\nr2\tFusarium oxysporum\n")
    (tmp_path / "cx.tsv").write_text("species\tcomplex\tcollapse\nFusarium oxysporum\tFOSC complex\t1\n")
    (tmp_path / "res.tsv").write_text(
        "query_id\tcall_species\tani\ttop_hit_species\tspecies_complex\twithheld\n"
        "q1\tFusarium alpha\t0.999\tFusarium alpha\t\t\n"
        "q2\tFusarium beta\t0.991\tFusarium beta\t\t\n"
        "q3\tFusarium alpha\t0.990\tFusarium alpha\t\t\n"
        "q4\t\t0.999\tFusarium oxysporum\tFOSC complex\tspecies: species_complex (x)\n")
    a = argparse.Namespace(results=str(tmp_path / "res.tsv"), truth=str(tmp_path / "truth.tsv"),
                           complexes=str(tmp_path / "cx.tsv"), out=str(tmp_path / "cal.tsv"),
                           meta=str(tmp_path / "meta.tsv"), min_precision=0.95)
    m.score(a)
    rows = {(r.split("\t")[0], r.split("\t")[1]): r.split("\t") for r in (tmp_path / "cal.tsv").read_text().splitlines()[1:]}
    hdr = (tmp_path / "cal.tsv").read_text().splitlines()[0].split("\t")
    at = lambda c, k: rows[("tef1", c)][hdr.index(k)]
    # cutoff 0.99: q1 right, q2 wrong -> precision 0.5; q3 (species absent) is a false call; q4 a right complex call
    assert at("0.99", "species_precision") == "0.5000" and at("0.99", "absent_false_call") == "1.0000"
    assert at("0.99", "complex_calls") == "1" and at("0.99", "complex_precision") == "1.0000"
    # cutoff 0.995: only q1 is called
    assert at("0.995", "species_precision") == "1.0000" and at("0.995", "absent_false_call") == "0.0000"
