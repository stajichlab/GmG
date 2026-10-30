"""Tests for the reads-validation scripts (#42, #45) on tiny synthetic inputs."""
import gzip
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def run(script, *args):
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)], capture_output=True, text=True, check=True).stdout


def test_reads_spike_adds_the_requested_share_of_foreign_pairs(tmp_path):
    """Foreign pairs are f / (100 - f) of the host pairs, host reads come first, names say donor and share."""
    rng = random.Random(1)

    def seq(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    for k in "12":
        with gzip.open(tmp_path / f"h_R{k}.fastq.gz", "wt") as fh:
            for i in range(200):
                fh.write(f"@h{i}\n{seq(150)}\n+\n{'F' * 150}\n")
    (tmp_path / "bac.fna").write_text(">b\n" + seq(4000) + "\n")
    (tmp_path / "hum.fna").write_text(">h\nNNNN" + seq(6000) + "\n")
    out = run("reads_spike.py", "hostA", tmp_path / "h_R1.fastq.gz", tmp_path / "h_R2.fastq.gz", tmp_path / "bac.fna", tmp_path / "hum.fna",
              tmp_path / "o", 200, "25,50")
    ids = {line.split("\t")[0] for line in out.splitlines()}
    assert ids == {"hostA__none_0", "hostA__bacteria_25", "hostA__bacteria_50", "hostA__human_25", "hostA__human_50"}
    heads = (tmp_path / "o" / "hostA__bacteria_50_R1.fastq").read_text().splitlines()[0::4]
    assert sum(h.startswith("@sim_") for h in heads) == 200 and sum(h.startswith("@h") for h in heads) == 200
    heads = (tmp_path / "o" / "hostA__human_25_R1.fastq").read_text().splitlines()[0::4]
    assert sum(h.startswith("@sim_") for h in heads) == 67   # 200 * 25 / 75
    seqs = (tmp_path / "o" / "hostA__human_50_R1.fastq").read_text().splitlines()[1::4]
    assert all("N" not in s for s in seqs)   # fragments over the N run are redrawn


def test_reads_summary_formats_a_row_per_run(tmp_path):
    rows = ["group\trank\tn\tcalled\tprecision\trecall", "all\tphylum\t10\t0.9000\t1.0000\t0.9000", "all\tgenus\t10\t0.5000\t0.8000\t0.4000"]
    (tmp_path / "a.tsv").write_text("\n".join(rows) + "\n")
    out = run("reads_summary.py", "all", f"x={tmp_path / 'a.tsv'}")
    assert "| x | 10 | 0.90 / 1.000 (9) |" in out and "0.50 / 0.800 (5)" in out


def test_name_truth_resolves_synonym_genera(tmp_path):
    """An old genus name resolves through a synonym to the current genus lineage; unknown genera are dropped."""
    td = tmp_path / "td"
    td.mkdir()
    (td / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tphylum\t|\n3\t|\t2\t|\tfamily\t|\n4\t|\t3\t|\tgenus\t|\n")
    (td / "names.dmp").write_text("2\t|\tPhyA\t|\t\t|\tscientific name\t|\n3\t|\tFamA\t|\t\t|\tscientific name\t|\n"
                                  "4\t|\tNewgenus\t|\t\t|\tscientific name\t|\n4\t|\tOldgenus\t|\t\t|\tsynonym\t|\n")
    (tmp_path / "labels.tsv").write_text("sample\tname\ns1\tOldgenus alpha\ns2\tNothere beta\n")
    (tmp_path / "meta.tsv").write_text("id\tgenus\nx\tNewgenus\n")
    run("name_truth.py", tmp_path / "labels.tsv", td, tmp_path / "t.tsv", "--index-meta", tmp_path / "meta.tsv")
    t = pd.read_csv(tmp_path / "t.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(t["id"]) == ["s1"] and t.loc[0, "genus"] == "Newgenus" and t.loc[0, "phylum"] == "PhyA" and t.loc[0, "genus_in_ref"] == "True"


def test_reads_vs_assembly_matches_strain_suffix(tmp_path):
    head = "query_id\tcall_phylum\tcall_genus\tstatus\n"
    (tmp_path / "r.tsv").write_text(head + "Genus_sp_NRRL_1\tP\tG1\tok\nGenus_sp_NRRL_2\tP\tG2\tok\n")
    (tmp_path / "a.tsv").write_text(head + "NRRL_1\tP\tG1\tok\n")
    (tmp_path / "t.tsv").write_text("id\tphylum\tclass\torder\tfamily\tgenus\tlabel_genus\nGenus_sp_NRRL_1\tP\t\t\t\tG1\tGenus\nGenus_sp_NRRL_2\tP\t\t\t\tG9\tGenus\n")
    out = run("reads_vs_assembly.py", tmp_path / "r.tsv", tmp_path / "a.tsv", tmp_path / "t.tsv")
    assert "1 of 2 read sets matched" in out and "genus: both called 1, same call 1" in out
