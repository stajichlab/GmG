"""Tests for the analysis scripts on tiny synthetic tables (no BFD data needed)."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_meta(path):
    rows = [  # id, taxid, phylum, class, order, family, genus, species
        ("a1", 1, "P", "C", "O", "F1", "G1", "G1 alpha"),
        ("a2", 1, "P", "C", "O", "F1", "G1", "G1 alpha"),
        ("a3", 1, "P", "C", "O", "F1", "G1", "G1 alpha"),
        ("b1", 2, "P", "C", "O", "F1", "G1", "G1 beta"),   # singleton species in genus G1
        ("c1", 3, "P", "C", "O", "F1", "G2", "G2 gamma"),
        ("c2", 3, "P", "C", "O", "F1", "G2", "G2 gamma"),
        ("d1", 4, "Q", "D", "N", "F9", "G9", "G9 delta"),  # alone in its phylum
    ]
    pd.DataFrame(rows, columns=["id", "taxid"] + RANKS).to_csv(path, sep="\t", index=False)


def bench_row(qid, truth, vote, ani=0.99, cons=1.0, top_id="x"):
    row = {"id": qid, "query_hashes": 1000, "top_id": top_id, "top_shared": 900, "containment_query": ani**21,
           "containment_ref": 0.9, "ani": ani, "n_ref": 7}
    for r, t, v in zip(RANKS, truth, vote):
        row.update({f"{r}_true": t, f"{r}_top": v, f"{r}_vote": v, f"{r}_cons": cons})
    return row


def test_answerable_rules(tmp_path):
    ab = load("analyze_bench")
    write_meta(tmp_path / "meta.tsv")
    ab.NC, ab.NSP, ab.NGE = ab.counts(str(tmp_path / "meta.tsv"))
    lin = {"a1": ("P", "C", "O", "F1", "G1", "G1 alpha"), "b1": ("P", "C", "O", "F1", "G1", "G1 beta"),
           "d1": ("Q", "D", "N", "F9", "G9", "G9 delta")}
    rows = []
    for mode in ("loo", "los", "log"):
        for qid, t in lin.items():
            r = bench_row(qid, t, t)
            r["mode"] = mode
            rows.append(r)
    d = ab.add_answerable(pd.DataFrame(rows))
    got = {(r["mode"], r["id"]): {k: bool(r[f"{k}_ans"]) for k in RANKS} for _, r in d.iterrows()}
    # leave-one-out: a1's species has other members, b1's species does not, d1 is alone at every rank
    assert got[("loo", "a1")]["species"] and not got[("loo", "b1")]["species"]
    assert got[("loo", "b1")]["genus"]
    assert not any(got[("loo", "d1")].values())
    # leave-species-out: species never answerable; genus answerable only if another species is in the genus
    assert not got[("los", "a1")]["species"] and got[("los", "a1")]["genus"]
    assert got[("los", "b1")]["genus"]
    # leave-genus-out: genus/species never answerable; family answerable (genus G2 remains), phylum d1 has no relative
    assert not got[("log", "a1")]["genus"] and got[("log", "a1")]["family"]
    assert not got[("log", "d1")]["phylum"]


def test_best_threshold_reaches_target():
    ab = load("analyze_bench")
    rng = np.random.default_rng(1)
    ani = rng.uniform(0.7, 1.0, 4000)
    ok = rng.uniform(size=4000) < np.where(ani > 0.9, 0.99, 0.5)   # reliable above ANI 0.9, a coin flip below
    x = pd.DataFrame({"ani": ani, "cons": 1.0, "voted": True, "vote_ok": ok})
    b = ab.best_threshold(x, 0.95)
    assert b is not None and b["precision"] >= 0.95 and b["coverage"] > 0.2
    assert 0.85 <= b["min_ani"] <= 0.95   # near the true change point
    x["vote_ok"] = rng.uniform(size=4000) < 0.5      # hopeless: no threshold reaches the target
    assert ab.best_threshold(x, 0.95) is None


def test_find_label_suspects_flags_conflict_and_hub(tmp_path):
    # 3 genomes (q1..q3) whose near-identical neighbour hub1 carries a different class; q4 is consistent
    rows = [bench_row(f"q{i}", ("P", "C", "O", "F", "G", f"S{i}"), ("P", "OTHER", "O", "F", "G", f"S{i}"),
                      ani=0.99, cons=1.0, top_id="hub1") for i in (1, 2, 3)]
    rows.append(bench_row("q4", ("P", "C", "O", "F", "G", "S4"), ("P", "C", "O", "F", "G", "S4"), top_id="ok1"))
    bench = tmp_path / "bench_loo.tsv"
    pd.DataFrame(rows).to_csv(bench, sep="\t", index=False)
    out = tmp_path / "labels"
    r = subprocess.run([sys.executable, "-I", str(ROOT / "scripts" / "find_label_suspects.py"), str(bench), str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    sus = pd.read_csv(f"{out}.suspects.tsv", sep="\t")
    assert set(sus[sus["rank"] == "class"]["id"]) == {"q1", "q2", "q3"}
    hubs = pd.read_csv(f"{out}.hubs.tsv", sep="\t")
    assert list(hubs["top_id"]) == ["hub1"] and int(hubs["n"].iloc[0]) == 3
    excl = set(Path(f"{out}.exclude.txt").read_text().split())
    assert {"q1", "q2", "q3", "hub1"} <= excl and "q4" not in excl
