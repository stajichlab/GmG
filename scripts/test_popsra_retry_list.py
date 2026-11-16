#!/usr/bin/env python3
"""Tests for popsra_retry_list.py, popsra_merge.py and popsra_manifest.py: python3 scripts/test_popsra_retry_list.py"""
import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import popsra_manifest as man  # noqa: E402
import popsra_merge as merge  # noqa: E402
import popsra_retry_list as rl  # noqa: E402

HEADER = "run\ttier\troute\tlayout\tbases_kept\treads_kept\tcov_t\tcov_r\tseconds\tverified\toutcome\tsampling\tfallback_from"


def result_line(run, outcome):
    return f"{run}\tcoverage\taws_odp\tPAIRED\t1\t1\t25\t25\t10\ttrue\t{outcome}\tstride\t"


class Tmp(unittest.TestCase):
    def setUp(self):
        self._d = tempfile.TemporaryDirectory()
        self.t = self._d.name

    def tearDown(self):
        self._d.cleanup()

    def write(self, name, text):
        p = os.path.join(self.t, name)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def retry(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            rc = rl.main(list(argv))
        return rc, out.getvalue(), err.getvalue()


class RetryPolicy(Tmp):
    def manifest(self, runs):
        return self.write("m.csv", "run_accession,study\n" + "".join(f"{r},S\n" for r in runs))

    def results(self, lines):
        return self.write("r.tsv", HEADER + "\n" + "".join(l + "\n" for l in lines))

    def retried(self, rc_out):
        return [r["run_accession"] for r in csv.DictReader(io.StringIO(rc_out[1]))]

    def test_each_outcome_goes_to_its_state(self):
        m = self.manifest(["A", "B", "C", "D", "E", "F", "G"])
        r = self.results([
            result_line("A", "ok"), result_line("B", "transient"), result_line("C", "hook_failed"),
            result_line("D", "unavailable_at_ncbi"), result_line("E", "verify_failed"),
        ])
        # F has no line, G has a line this script does not know
        r2 = self.write("r2.tsv", result_line("G", "brand_new") + "\n")
        out = self.retry(m, r, r2)
        self.assertEqual(self.retried(out), ["B", "E", "F"])

    def test_an_unknown_outcome_is_set_aside_not_retried(self):
        m = self.manifest(["G"])
        r = self.results([result_line("G", "brand_new")])
        out = self.retry(m, r)
        self.assertEqual(self.retried(out), [])
        self.assertIn("inspect\t1", out[2])

    def test_a_transient_run_is_retried_only_up_to_the_attempt_limit(self):
        m = self.manifest(["T"])
        for n, want in ((1, ["T"]), (2, ["T"]), (3, [])):
            r = self.results([result_line("T", "transient")] * n)
            self.assertEqual(self.retried(self.retry(m, r)), want, f"{n} attempts")
        r = self.results([result_line("T", "transient")] * 3)
        self.assertEqual(self.retried(self.retry(m, r, "--max-attempts", "5")), ["T"])

    def test_two_verify_failures_are_quarantined(self):
        m = self.manifest(["V"])
        self.assertEqual(self.retried(self.retry(m, self.results([result_line("V", "verify_failed")]))), ["V"])
        two = self.results([result_line("V", "verify_failed")] * 2)
        self.assertEqual(self.retried(self.retry(m, two)), [])

    def test_the_last_line_decides(self):
        m = self.manifest(["R"])
        r = self.results([result_line("R", "transient"), result_line("R", "ok")])
        self.assertEqual(self.retried(self.retry(m, r)), [])

    def test_a_missing_results_file_means_everything_is_retried(self):
        m = self.manifest(["A", "B"])
        out = self.retry(m, os.path.join(self.t, "nope.tsv"))
        self.assertEqual(out[0], 0)
        self.assertEqual(self.retried(out), ["A", "B"])

    def test_a_results_file_with_a_repeated_or_missing_header_does_not_crash(self):
        m = self.manifest(["A", "B"])
        two_headers = self.write("th.tsv", f"{HEADER}\n{result_line('A', 'ok')}\n{HEADER}\n{result_line('B', 'ok')}\n")
        self.assertEqual(self.retried(self.retry(m, two_headers)), [])
        no_header = self.write("nh.tsv", result_line("A", "ok") + "\n" + result_line("B", "ok") + "\n")
        out = self.retry(m, no_header)  # read with the default columns: no line is lost
        self.assertEqual((out[0], self.retried(out)), (0, []))

    def test_a_manifest_with_no_rows_left_is_not_an_error(self):
        m = self.write("empty.csv", "run_accession,study\n")
        out = self.retry(m, self.results([]))
        self.assertEqual((out[0], self.retried(out)), (0, []))
        self.assertTrue(out[1].startswith("run_accession,study"))

    def test_a_manifest_without_the_column_is_exit_2(self):
        m = self.write("bad.csv", "acc\nA\n")
        self.assertEqual(self.retry(m, self.results([]))[0], 2)

    def test_exit_status_tells_a_driver_whether_runs_remain(self):
        m = self.manifest(["A", "B"])
        self.assertEqual(self.retry(m, self.results([result_line("A", "ok")]), "--exit-status")[0], 10)
        done = self.results([result_line("A", "ok"), result_line("B", "ok")])
        self.assertEqual(self.retry(m, done, "--exit-status")[0], 0)
        self.assertEqual(self.retry(m, self.results([result_line("A", "ok")]))[0], 0)  # without the flag: always 0

    def test_a_quoted_field_with_a_comma_survives_the_retry_manifest(self):
        m = self.write("q.csv", 'run_accession,geo\nSRR1,"USA: California, Riverside"\nSRR2,x\n')
        out = self.retry(m, self.results([result_line("SRR2", "ok")]))
        rows = list(csv.DictReader(io.StringIO(out[1])))
        self.assertEqual([(r["run_accession"], r["geo"]) for r in rows], [("SRR1", "USA: California, Riverside")])
        # and the chunk reader gets the same accession back
        p = self.write("retry.csv", out[1])
        self.assertEqual(man.runs(p), ["SRR1"])


class Manifest(Tmp):
    def test_a_quoted_comma_does_not_shift_the_columns(self):
        p = self.write("m.csv", 'sample,geo,run_accession\nS1,"USA: California, Riverside",SRR001\nS2,x,SRR002\n')
        self.assertEqual(man.runs(p), ["SRR001", "SRR002"])

    def test_tsv_blank_rows_and_a_header_only_file(self):
        p = self.write("m.tsv", "run_accession\tstudy\nA\tS\n\nB\tS\n")
        self.assertEqual(man.runs(p), ["A", "B"])
        self.assertEqual(man.runs(self.write("h.csv", "run_accession,study\n")), [])

    def test_no_column_is_an_error_and_the_program_exits_2(self):
        p = self.write("bad.csv", "acc\nA\n")
        with self.assertRaises(ValueError):
            man.runs(p)
        rc = subprocess.run([sys.executable, os.path.join(HERE, "popsra_manifest.py"), p], capture_output=True).returncode
        self.assertEqual(rc, 2)


class Merge(Tmp):
    def test_task_files_are_joined_in_task_order_with_one_header_and_nothing_modified(self):
        for task, runs in (("10", ["C"]), ("2", ["B"]), ("0", ["A"])):
            self.write(f"w/tasks/{task}/results.tsv", HEADER + "\n" + "".join(result_line(r, "ok") + "\n" for r in runs))
            self.write(f"w/tasks/{task}/exits.tsv", "".join(f"{r}\t0\tT\n" for r in runs))
        self.write("w/results.tsv", HEADER + "\n" + result_line("OLD", "ok") + "\n")  # the flat layout of an earlier run
        self.assertEqual(merge.main([os.path.join(self.t, "w")]), 0)
        lines = open(os.path.join(self.t, "w/merged/results.tsv")).read().splitlines()
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([l.split("\t")[0] for l in lines[1:]], ["OLD", "A", "B", "C"])  # flat first, tasks 0, 2, 10
        self.assertEqual(open(os.path.join(self.t, "w/merged/exits.tsv")).read().split(), "A 0 T B 0 T C 0 T".split())
        self.assertEqual(len(open(os.path.join(self.t, "w/results.tsv")).read().splitlines()), 2)  # input untouched
        # running it again gives the same file (the output is not an input)
        before = open(os.path.join(self.t, "w/merged/results.tsv")).read()
        merge.main([os.path.join(self.t, "w")])
        self.assertEqual(open(os.path.join(self.t, "w/merged/results.tsv")).read(), before)

    def test_a_task_run_twice_does_not_repeat_its_header_and_a_name_with_no_input_is_skipped(self):
        self.write("w/tasks/0/results.tsv", HEADER + "\n" + result_line("A", "transient") + "\n" + HEADER + "\n" + result_line("A", "ok") + "\n")
        merge.main([os.path.join(self.t, "w")])
        text = open(os.path.join(self.t, "w/merged/results.tsv")).read()
        self.assertEqual(text.count("run\ttier"), 1)
        self.assertFalse(os.path.exists(os.path.join(self.t, "w/merged/profile.tsv")))

    def test_merged_results_feed_the_retry_list(self):
        self.write("w/tasks/0/results.tsv", HEADER + "\n" + result_line("A", "ok") + "\n")
        self.write("w/tasks/1/results.tsv", HEADER + "\n" + result_line("B", "transient") + "\n")
        merge.main([os.path.join(self.t, "w")])
        m = self.write("m.csv", "run_accession\nA\nB\nC\n")
        out = io.StringIO()
        with redirect_stdout(out), redirect_stderr(io.StringIO()):
            rl.main([m, os.path.join(self.t, "w/merged/results.tsv")])
        self.assertEqual([r["run_accession"] for r in csv.DictReader(io.StringIO(out.getvalue()))], ["B", "C"])


if __name__ == "__main__":
    unittest.main(verbosity=1)
