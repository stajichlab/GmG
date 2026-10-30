"""The Python bindings must give the records `gmg run` writes (same Rust code). Needs the `gmg` binary (env GMG_BIN, default
`gmg` on PATH) to build a small reference from random genomes."""
import gzip, json, os, random, subprocess

import pytest

import gmg

GMG = os.environ.get("GMG_BIN", "gmg")
LINEAGE = [
    "Ascomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tAspergillus alpha",
    "Ascomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tAspergillus beta",
    "Ascomycota\tSordariomycetes\tHypocreales\tNectriaceae\tFusarium\tFusarium gamma",
    "Ascomycota\tSordariomycetes\tHypocreales\tNectriaceae\tFusarium\tFusarium delta",
    "Basidiomycota\tAgaricomycetes\tAgaricales\tAgaricaceae\tAgaricus\tAgaricus eps",
    "Basidiomycota\tAgaricomycetes\tAgaricales\tAgaricaceae\tAgaricus\tAgaricus zeta",
]


def fasta(name, seq):
    return f">{name}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n"


def run(*args):
    r = subprocess.run([GMG, *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("gmg")
    rng = random.Random(1)
    genomes = ["".join(rng.choice("ACGT") for _ in range(200_000)) for _ in range(6)]
    meta = ["id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies"]
    refs = []
    (d / "refs").mkdir()
    for i, g in enumerate(genomes):
        p = d / "refs" / f"g{i}.fa"
        p.write_text(fasta(f"g{i}", g))
        refs.append(str(p))
        meta.append(f"g{i}\t{100 + i}\t{LINEAGE[i]}")
    (d / "meta.tsv").write_text("\n".join(meta) + "\n")
    run("sketch", "-k", "21", "-s", "10", "-o", str(d / "sk"), *refs)
    run("build", "--meta", str(d / "meta.tsv"), "--sketches", str(d / "sk"), "-o", str(d / "ix.gmgi"))
    (d / "in").mkdir()
    (d / "in" / "copy.fa").write_text(fasta("c", genomes[2]))
    (d / "in" / "other.fa.gz").write_bytes(gzip.compress(fasta("o", genomes[4]).encode()))
    (d / "in" / "random.fa").write_text(fasta("r", "".join(rng.choice("ACGT") for _ in range(150_000))))
    (d / "in" / "empty.fa").write_text("")
    run("run", str(d / "in"), "--index", str(d / "ix.gmgi"), "-o", str(d / "out"), "-q", "--no-report")
    cli = {}
    for line in (d / "out" / "results.jsonl").read_text().splitlines():
        rec = json.loads(line)
        cli[rec["id"]] = rec
    return d, genomes, cli


def test_index_properties_and_manifest(world):
    d, _, _ = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    assert (ix.ksize, ix.scaled, ix.n_genomes) == (21, 10, 6)
    assert isinstance(ix.manifest(), dict) and ix.manifest()["genomes"] == 6
    assert ix.thresholds == "known"
    ix.use_thresholds("novel")
    assert ix.thresholds == "novel"
    with pytest.raises(ValueError):
        ix.use_thresholds("nope")


@pytest.mark.parametrize("name", ["copy", "other", "random"])
def test_classify_equals_gmg_run(world, name):
    d, _, cli = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    path = next((d / "in").glob(f"{name}.fa*"))
    rec = ix.classify(str(path))
    want = cli[name]
    for k in ("id", "status", "best_rank", "kind", "guess", "qc"):
        assert rec.get(k) == want.get(k), k


def test_species_call_and_gzip(world):
    d, _, _ = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    assert ix.classify(str(d / "in" / "copy.fa"))["guess"]["calls"]["species"]["taxon"] == "Fusarium gamma"
    assert ix.classify(str(d / "in" / "other.fa.gz"))["guess"]["calls"]["species"]["taxon"] == "Agaricus eps"


def test_failures_are_records_not_exceptions(world):
    d, _, _ = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    rec = ix.classify(str(d / "in" / "empty.fa"))
    assert rec["status"] == "error" and rec["error_code"]
    with pytest.raises(OSError):
        ix.classify(str(d / "in" / "missing.fa"))
    with pytest.raises(ValueError):
        gmg.Index(str(d / "in" / "copy.fa"))  # not an index


def test_sketch_and_guess_hashes(world):
    d, genomes, _ = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    h = ix.sketch(fasta("q", genomes[3]).encode())
    assert len(h) > 1000 and h == sorted(set(h))
    g = ix.guess_hashes(h)
    assert g["calls"]["species"]["taxon"] == "Fusarium delta"


def test_tsv_text(world):
    d, _, cli = world
    ix = gmg.Index(str(d / "ix.gmgi"))
    recs = [ix.classify(str(d / "in" / "copy.fa"))]
    tsv = gmg.results_tsv_text(recs)
    lines = tsv.strip().split("\n")
    assert lines[0].startswith("query_id\t") and lines[1].startswith("copy\t")
    assert gmg.hits_tsv_text(recs).startswith("query_id\t")
