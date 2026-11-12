"""scripts/popsra_retry_list.py: which runs of a manifest need another go (PopSeq spec 002's exit codes)."""
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load():
    spec = importlib.util.spec_from_file_location("rl", ROOT / "scripts" / "popsra_retry_list.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write(tmp, runs, lines, delim=","):
    m = tmp / "manifest.csv"
    m.write_text(delim.join(["run_accession", "study"]) + "\n" + "".join(f"{r}{delim}S{i}\n" for i, r in enumerate(runs)))
    res = tmp / "results.tsv"
    head = ["run", "tier", "route", "layout", "bases_kept", "reads_kept", "coverage_target", "coverage_reached", "seconds", "verified", "outcome"]
    res.write_text("\t".join(head) + "\n" + "".join("\t".join([r, "coverage", "aws_odp", "PAIRED", "", "", "", "", "1", "x", o]) + "\n" for r, o in lines))
    return m, res


def test_each_outcome_goes_to_the_documented_state(tmp_path):
    rl = load()
    runs = ["ok1", "tr1", "vf1", "hf1", "un1", "up1", "gone"]
    m, res = write(tmp_path, runs, [("ok1", "ok"), ("tr1", "transient"), ("vf1", "verify_failed"), ("hf1", "hook_failed"), ("un1", "unavailable_at_ncbi"), ("up1", "unpaired")])
    rows, _ = rl.read_table(str(res))
    st = rl.classify(runs, rows)
    assert st == {"ok1": "done", "tr1": "retry", "vf1": "retry", "hf1": "inspect", "un1": "final", "up1": "final", "gone": "retry"}


def test_the_last_line_decides_and_a_repeated_verify_failure_is_set_aside(tmp_path):
    rl = load()
    runs = ["a", "b", "c"]
    m, res = write(tmp_path, runs, [("a", "transient"), ("a", "ok"), ("b", "verify_failed"), ("b", "verify_failed"), ("c", "ok"), ("c", "transient")])
    rows, _ = rl.read_table(str(res))
    st = rl.classify(runs, rows)
    assert st == {"a": "done", "b": "quarantine", "c": "retry"}


def test_the_retry_manifest_keeps_rows_and_columns_in_manifest_order(tmp_path, capsys):
    rl = load()
    runs = ["r3", "r1", "r2", "r4"]
    m, res = write(tmp_path, runs, [("r1", "ok"), ("r2", "transient"), ("r4", "verify_failed")])
    out = tmp_path / "retry.csv"
    assert rl.main([str(m), str(res), "-o", str(out)]) == 0
    rows = list(csv.DictReader(open(out)))
    assert [r["run_accession"] for r in rows] == ["r3", "r2", "r4"]  # r3 never ran; r1 is done
    assert set(rows[0]) == {"run_accession", "study"}
    err = capsys.readouterr().err
    assert "done\t1" in err and "retry\t3" in err


def test_a_tab_separated_manifest_works_and_a_bad_manifest_is_refused(tmp_path, capsys):
    rl = load()
    m, res = write(tmp_path, ["x", "y"], [("x", "ok")], delim="\t")
    out = tmp_path / "r.tsv"
    assert rl.main([str(m), str(res), "-o", str(out)]) == 0
    assert out.read_text().splitlines()[0] == "run_accession\tstudy" and "y\tS1" in out.read_text()
    bad = tmp_path / "bad.csv"
    bad.write_text("acc\nx\n")
    assert rl.main([str(bad), str(res)]) == 2
