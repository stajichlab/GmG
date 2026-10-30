"""mix_reads.py and mix_score.py (#45) on tiny synthetic inputs."""
import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def fastq(path, tag, n):
    with gzip.open(path, "wt") as f:
        for i in range(n):
            f.write(f"@{tag}{i}\nACGT\n+\nFFFF\n")


def run(*args):
    r = subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)
    return r


def test_mix_reads_takes_the_first_pairs_of_each_parent(tmp_path):
    for t in "AB":
        for k in ("R1", "R2"):
            fastq(tmp_path / f"{t}_{k}.fastq.gz", f"{t}{k}_", 1000)
    r = run(ROOT / "scripts/mix_reads.py", tmp_path, "m90", 100, "90,10",
            tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz", tmp_path / "B_R1.fastq.gz", tmp_path / "B_R2.fastq.gz")
    assert r.returncode == 0, r.stderr
    assert len(r.stdout.strip().split("\n")) == 2
    heads = [l for l in (tmp_path / "m90_R1.fastq").read_text().split("\n") if l.startswith("@")]
    assert len(heads) == 100
    assert sum(h.startswith("@AR1_") for h in heads) == 90 and sum(h.startswith("@BR1_") for h in heads) == 10
    assert heads[0] == "@AR1_0" and heads[90] == "@BR1_0"
    # R1 and R2 have the same number of pairs
    assert len([l for l in (tmp_path / "m90_R2.fastq").read_text().split("\n") if l.startswith("@")]) == 100
    # a pure control and a three-parent mixture
    r = run(ROOT / "scripts/mix_reads.py", tmp_path, "pure", 50, "100,0",
            tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz", tmp_path / "B_R1.fastq.gz", tmp_path / "B_R2.fastq.gz")
    assert r.returncode == 0 and "@BR1_" not in (tmp_path / "pure_R1.fastq").read_text()
    r = run(ROOT / "scripts/mix_reads.py", tmp_path, "t", 100, "60,30,10", *[tmp_path / f"{t}_{k}.fastq.gz" for t in "ABA" for k in ("R1", "R2")])
    assert r.returncode == 0
    assert len([l for l in (tmp_path / "t_R1.fastq").read_text().split("\n") if l.startswith("@")]) == 100


def test_mix_reads_refuses_a_parent_that_is_too_short(tmp_path):
    for t in "AB":
        for k in ("R1", "R2"):
            fastq(tmp_path / f"{t}_{k}.fastq.gz", t, 10)
    r = run(ROOT / "scripts/mix_reads.py", tmp_path, "x", 100, "50,50",
            tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz", tmp_path / "B_R1.fastq.gz", tmp_path / "B_R2.fastq.gz")
    assert r.returncode != 0 and "needed" in r.stderr


def table(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def test_mix_score_counts_majority_minority_and_the_minority_in_the_hits(tmp_path):
    A = dict(zip(RANKS, ["Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", ""]))
    B = dict(zip(RANKS, ["Chytridiomycota", "Chytridiomycetes", "Chytridiales", "Chytriomycetaceae", "Chytriomyces", ""]))
    ph = ["pair_id", "distance", "high", "a_id", "a_r1", "a_r2", "b_id", "b_r1", "b_r2"] + [f"a_{r}" for r in RANKS] + [f"b_{r}" for r in RANKS]
    table(tmp_path / "plan.tsv", ph, [["p00", "diff_phylum", 0, "a", "x", "x", "b", "x", "x"] + [A[r] for r in RANKS] + [B[r] for r in RANKS]])
    table(tmp_path / "trip.tsv", ["triple_id", "type"], [])
    rh = ["query_id", "status", "ani", "qc_flags"] + [f"call_{r}" for r in RANKS]

    def res(vid, ani, call):
        return [vid, "ok", ani, "weak_phylum" if ani < 0.85 else ""] + [call.get(r, "") for r in RANKS]

    rows = [res("p00__low__A100", 0.99, A), res("p00__low__A90", 0.99 * 0.9 ** (1 / 21), A), res("p00__low__A50", 0.80, B), res("p00__low__A0", 0.99, B)]
    table(tmp_path / "res.tsv", rh, rows)
    hh = ["query_id", "hit_rank", "hit_id"] + RANKS
    hits = []
    for vid in ("p00__low__A90", "p00__low__A50"):
        hits.append([vid, 1, "h1"] + [A[r] for r in RANKS])
        hits.append([vid, 2, "h2"] + [B[r] for r in RANKS])
    table(tmp_path / "hits.tsv", hh, hits)
    r = run(ROOT / "scripts/mix_score.py", tmp_path / "plan.tsv", tmp_path / "trip.tsv", "--run", f"t={tmp_path / 'res.tsv'}:{tmp_path / 'hits.tsv'}")
    assert r.returncode == 0, r.stderr
    rows = {tuple(l.split("|")[1:4]): l for l in r.stdout.split("\n") if l.startswith("| diff_phylum")}
    a90 = rows[(" diff_phylum ", " phylum ", " 90% ")]
    cells = [c.strip() for c in a90.split("|")[1:-1]]
    assert cells[4] == "100%" and cells[5] == "0%" and cells[8] == "100%"   # majority, minority, minority in the top-10 hits
    assert cells[9].split(" / ")[0] == cells[9].split(" / ")[1]               # observed ANI ratio equals p^(1/21) by construction
    a50 = [c.strip() for c in rows[(" diff_phylum ", " phylum ", " 50% ")].split("|")[1:-1]]
    assert a50[4] == "100%" and a50[5] == "0%"   # at 50:50 either parent is a correct answer (here B)
    assert a50[11] == "1"                          # the weak_phylum flag was counted
