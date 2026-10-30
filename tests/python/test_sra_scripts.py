"""Tests for the real-reads (SRA/ENA) scripts of #42 on tiny synthetic inputs; no network."""
import csv
import gzip
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(script, *args):
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)], capture_output=True, text=True, check=True)


def write_tsv(path, rows):
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "\n".join("\t".join(str(r[c]) for c in cols) for r in rows) + "\n")


def test_sra_plan_keeps_only_runs_whose_species_three_records_agree_on(tmp_path):
    meta = [{"id": f"GCA_00000000{i}.1_asm", "taxid": 1, "phylum": "P", "class": "C", "order": "O", "family": "F", "genus": "Gen", "species": sp}
            for i, sp in [(1, "Gen alpha"), (2, "Gen beta"), (3, "Gen sp. X"), (4, "Gen delta"), (5, "Gen alpha")]]
    write_tsv(tmp_path / "meta.tsv", meta)
    asm = [("GCA_000000001.1", "SAMN1", "Gen alpha"), ("GCA_000000002.1", "SAMN2", "Gen gamma"),   # organism name disagrees
           ("GCA_000000003.1", "SAMN3", "Gen sp. X"), ("GCA_000000004.1", "SAMN4", "Gen delta"), ("GCA_000000005.1", "SAMN5", "Gen alpha")]
    (tmp_path / "asm.txt").write_text("# comment\n#assembly_accession\tbiosample\torganism_name\tgenome_size\n"
                                      + "\n".join(f"{a}\t{b}\t{o}\t30000000" for a, b, o in asm) + "\n")
    f = lambda r: f"ftp.sra.ebi.ac.uk/vol1/{r}_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/{r}_2.fastq.gz"
    ena = [dict(run_accession="R1", sample_accession="SAMN1", secondary_sample_accession="", scientific_name="Gen alpha", center_name="C1",
                instrument_model="HiSeq", read_count=900000, base_count=270000000, fastq_ftp=f("R1")),
           dict(run_accession="R2", sample_accession="SAMN2", secondary_sample_accession="", scientific_name="Gen beta", center_name="C2",
                instrument_model="HiSeq", read_count=900000, base_count=270000000, fastq_ftp=f("R2")),      # reference label Gen beta, NCBI Gen gamma
           dict(run_accession="R3", sample_accession="SAMN3", secondary_sample_accession="", scientific_name="Gen sp. X", center_name="C2",
                instrument_model="HiSeq", read_count=900000, base_count=270000000, fastq_ftp=f("R3")),      # "sp."
           dict(run_accession="R4", sample_accession="SAMN4", secondary_sample_accession="", scientific_name="Gen delta", center_name="C3",
                instrument_model="HiSeq", read_count=100000, base_count=30000000, fastq_ftp=f("R4")),       # too few pairs
           dict(run_accession="R5", sample_accession="SAMEA9", secondary_sample_accession="SAMN5", scientific_name="Gen alpha", center_name="C3",
                instrument_model="NovaSeq", read_count=900000, base_count=270000000, fastq_ftp=f("R5"))]    # matched through the secondary accession
    write_tsv(tmp_path / "ena.tsv", ena)
    (tmp_path / "excl.txt").write_text("")
    run("sra_plan.py", tmp_path / "ena.tsv", tmp_path / "asm.txt", tmp_path / "meta.tsv", tmp_path / "excl.txt", tmp_path / "out.tsv", 10)
    rows = list(csv.DictReader(open(tmp_path / "out.tsv"), delimiter="\t"))
    assert sorted(r["run"] for r in rows) == ["R1", "R5"]
    r1 = next(r for r in rows if r["run"] == "R1")
    assert r1["bfd_id"] == "GCA_000000001.1_asm" and r1["read_len"] == "150" and r1["species_others_in_ref"] == "1" and r1["genome_size"] == "30000000"


def test_sra_harvest_reads_only_the_head_and_counts_bytes(tmp_path):
    rng = random.Random(2)
    for k in "12":
        with gzip.open(tmp_path / f"RUN1_{k}.src.fastq.gz", "wt") as fh:
            for i in range(60000):
                s = "".join(rng.choice("ACGT") for _ in range(100))
                fh.write(f"@r{i}/{k}\n{s}\n+\n{'F' * 100}\n")
    runs = tmp_path / "runs.tsv"
    write_tsv(runs, [dict(run="RUN1", fastq_1=f"file://{tmp_path}/RUN1_1.src.fastq.gz", fastq_2=f"file://{tmp_path}/RUN1_2.src.fastq.gz")])
    run("sra_harvest.py", runs, tmp_path / "raw", 300, 2, 1e9)
    for k in "12":
        lines = gzip.open(tmp_path / "raw" / f"RUN1_{k}.fastq.gz", "rt").read().splitlines()
        assert len(lines) == 1200 and lines[0] == f"@r0/{k}" and lines[-4] == f"@r299/{k}"
    b = (tmp_path / "raw" / "bytes.tsv").read_text().splitlines()[1].split("\t")
    full = (tmp_path / "RUN1_1.src.fastq.gz").stat().st_size
    assert b[0] == "RUN1" and 0 < int(b[1]) < full / 2 and int(b[3]) == 300   # stopped long before the end of the file
    assert (tmp_path / "raw" / "RUN1.done").exists()
    before = (tmp_path / "raw" / "RUN1_1.fastq.gz").stat().st_mtime_ns
    run("sra_harvest.py", runs, tmp_path / "raw", 300, 2, 1e9)           # a finished run is not fetched again
    assert (tmp_path / "raw" / "RUN1_1.fastq.gz").stat().st_mtime_ns == before


def test_sra_depth_table_counts_recall_over_attainable_inputs(tmp_path):
    """Recall is over inputs whose truth can still be found in the reference after the exclusion; the minimum depth comes from it."""
    write_tsv(tmp_path / "runs.tsv", [dict(run="A", read_len=100, genome_size=10_000_000, species_others_in_ref=1, genus_others_in_ref=2),
                                      dict(run="B", read_len=100, genome_size=10_000_000, species_others_in_ref=0, genus_others_in_ref=0)])
    write_tsv(tmp_path / "truth.tsv", [dict(id="A", genus="G1", species="G1 a", phylum="P", center="X"), dict(id="B", genus="G2", species="G2 b", phylum="P", center="X")])
    hdr = ["query_id", "call_phylum", "call_genus", "call_species"]
    write_tsv(tmp_path / "lo.tsv", [dict(zip(hdr, ["A", "P", "", ""])), dict(zip(hdr, ["B", "", "", ""]))])
    write_tsv(tmp_path / "hi.tsv", [dict(zip(hdr, ["A", "P", "G1", "G1 a"])), dict(zip(hdr, ["B", "P", "G2", "G9 x"]))])
    out = run("sra_depth_table.py", tmp_path / "runs.tsv", tmp_path / "truth.tsv", tmp_path / "o.md",
              f"self:pair:5000={tmp_path / 'lo.tsv'}", f"self:pair:500000={tmp_path / 'hi.tsv'}").stdout
    # B's genus is not elsewhere in the reference (genus_others 0): only A counts for recall, so recall is 0 then 1 even though B is also called right
    assert "| 5,000 (median 0.10x) | 2 | 1 | 0.00 | - (0) | 0.00 |" in out
    assert "| 500,000 (median 10.00x) | 2 | 1 | 1.00 | 1.000 (2) | 1.00 |" in out
    assert "Minimum depth for genus (recall over attainable inputs >= 0.8 from there on): 500,000 pairs (median 10.00x)" in out
    assert "Minimum depth for species (recall over attainable inputs >= 0.8 from there on): 500,000 pairs" in out
    assert "**By phylum, 500,000 pairs**" in out and "| P | 2 | 1 |" in out and "centers with >= 3 runs" in out


def test_sra_harvest_resumes_a_stream_the_server_ends_early(tmp_path):
    """A server that closes every response after 60 kB (HTTP Range supported) must still deliver the requested pairs, without re-reading bytes."""
    import threading
    from http.server import BaseHTTPRequestHandler, HTTPServer
    rng = random.Random(5)
    src = tmp_path / "big.fastq.gz"
    with gzip.open(src, "wt") as fh:
        for i in range(20000):
            s = "".join(rng.choice("ACGT") for _ in range(100))
            fh.write(f"@r{i}\n{s}\n+\n{'F' * 100}\n")
    blob = src.read_bytes()
    hits = []

    class H(BaseHTTPRequestHandler):
        def do_GET(self):
            start = int(self.headers["Range"].split("=")[1].rstrip("-")) if self.headers.get("Range") else 0
            end = min(len(blob), start + 60_000)
            hits.append((self.path, start))
            self.send_response(206 if start else 200)
            self.send_header("Content-Length", str(end - start))
            self.end_headers()
            self.wfile.write(blob[start:end])

        def log_message(self, *a):
            pass

    srv = HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{srv.server_address[1]}/big.fastq.gz"
    runs = tmp_path / "runs.tsv"
    write_tsv(runs, [dict(run="RUN2", fastq_1=url, fastq_2=url.replace("big", "big2"))])
    run("sra_harvest.py", runs, tmp_path / "raw", 4000, 1, 1e9)
    srv.shutdown()
    lines = gzip.open(tmp_path / "raw" / "RUN2_1.fastq.gz", "rt").read().splitlines()
    assert len(lines) == 16000 and lines[0] == "@r0" and lines[-4] == "@r3999"
    h1 = [b for p, b in hits if p == "/big.fastq.gz"]
    assert len(h1) > 2 and h1 == sorted(h1) and len(set(h1)) == len(h1)       # several resumes, each at a later byte
    row = (tmp_path / "raw" / "bytes.tsv").read_text().splitlines()[1].split("\t")
    assert int(row[3]) == 4000 and int(row[5]) >= 2 and int(row[1]) <= max(h1) + 60_000
