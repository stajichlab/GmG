#!/usr/bin/env python3
"""Tests for popset_select.py: python3 scripts/test_popset_select.py"""
import csv
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import popset_select as ps  # noqa: E402


class Classify(unittest.TestCase):
    def cat(self, text, has=True):
        return ps.cultivation_of(text.lower(), has)[0]

    def test_the_named_cultivations(self):
        for text, want in [
            ("isolation_source: sake brewery", "sake_rice"), ("soy sauce moromi", "soy_other_ferment"), ("beer", "beer"),
            ("spontaneous wine fermentation", "wine"), ("vineyard", "wine"), ("fermenting dough", "bread_commercial"),
            ("active dried yeast packet", "bread_commercial"), ("liquid yeast", "bread_commercial"), ("fermented sorghum grain", "soy_other_ferment"), ("flor velum", "wine"), ("cachaca distillery", "spirits_ethanol"), ("sugarcane bioethanol plant", "spirits_ethanol"),
            ("clinical isolate from blood", "clinical"), ("oak bark", "wild_nature"), ("soil under oak", "wild_nature"),
            ("a thing nobody classified", "other_known"),
        ]:
            self.assertEqual(self.cat(text), want, text)

    def test_no_text_is_unannotated_not_other(self):
        self.assertEqual(ps.cultivation_of("", False)[0], "unannotated")

    def test_the_more_specific_cultivation_wins(self):
        self.assertEqual(self.cat("rice wine"), "sake_rice")  # not wine
        self.assertEqual(self.cat("beer brewery in a vineyard region"), "beer")

    def test_laboratory_and_experimental_samples_are_recognised(self):
        for t in ["laboratory strain", "evolution experiment", "chemostat", "haploid cells from tetrads", "by4741 derivative", "ypd plate",
                  "engineered strain", "pyes2-empty", "ko171 x ko124", "dz12-derived isolate", "sc-his 2mm",
                  "d141 deltaku70 no. 60-4", "nr-166 (cts2/ard1/cts3)", "derived from coccidioides posadasii strain silveira", "ku70\u0394 mutant", "aspergillus fumigatus akub deletion"]:
            self.assertTrue(ps.is_lab(t.lower()), t)
        # "SC-2010" in a collection sample name is S. cerevisiae, not a medium
        self.assertFalse(ps.is_lab("150a_sc_dbvpg1106-sc-2010-11-01-9578"))
        for t in ["wine fermentation in a vineyard", "oak bark", "sake brewery", "clinical isolate", "vpci 1882/14", "vpci 40/ei/13/3/b", "unknown patient"]:
            self.assertFalse(ps.is_lab(t.lower()), t)
        self.assertFalse(ps.is_lab("collaborative effort"))  # whole words only: not "lab" inside another word


class Geography(unittest.TestCase):
    def test_places_that_differ_only_in_spacing_are_one_place(self):
        self.assertEqual(ps.place_of({"geo_loc_name": "USA:Washington", "country_iso": "US"}, True),
                         ps.place_of({"geo_loc_name": "USA: Washington ", "country_iso": "US"}, True))
        self.assertEqual(ps.place_of({"geo_loc_name": "", "country_iso": "US"}, True), "US")
        self.assertEqual(ps.place_of({"geo_loc_name": "USA: Utah", "country_iso": "US"}, False), "US")

    def test_regions(self):
        self.assertEqual(ps.region_of("US"), "North America")
        self.assertEqual(ps.region_of("it"), "Europe")
        self.assertEqual(ps.region_of("CN"), "Asia")
        self.assertEqual(ps.region_of("ZA"), "Africa")
        self.assertEqual(ps.region_of("BR"), "South America")
        self.assertEqual(ps.region_of("AU"), "Oceania")
        self.assertEqual(ps.region_of(""), "unknown_region")
        self.assertEqual(ps.region_of("XX"), "unknown_region")


def make(tmp, rows, attrs):
    cand = os.path.join(tmp, "cand.csv")
    cols = ["run_accession", "sample_accession", "study_accession", "species_name", "strain", "country_iso", "country", "collection_date",
            "est_depth_x", "base_count", "library_layout", "fastq_ftp"]
    with open(cand, "w", newline="") as fh:
        fh.write("# popseq_schema: 1.1.0\n")
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    at = os.path.join(tmp, "attrs.csv")
    acols = ["biosample", "isolation_source", "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "host", "isolate",
             "strain", "culture_collection", "sample_name"]
    with open(at, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=acols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(attrs)
    return cand, at


def fixture(n_per=30):
    rows, attrs = [], []
    k = 0
    for cat, src in [("wine", "wine fermentation"), ("beer", "beer"), ("clinical", "clinical isolate"), ("wild", "oak bark")]:
        for i in range(n_per):
            k += 1
            run, bs = f"RUN{k:04d}", f"SAMN{k:04d}"
            iso = ["IT", "FR", "US", "CN", "ZA", "BR"][i % 6]
            rows.append({"run_accession": run, "sample_accession": bs, "study_accession": f"PRJ{(i // 5) % 4}{cat[0]}", "species_name": "X y",
                         "country_iso": iso, "est_depth_x": str(30 + i)})
            attrs.append({"biosample": bs, "isolation_source": src})
    # a laboratory sample and a sample without any text
    rows.append({"run_accession": "LAB1", "sample_accession": "SAMNLAB", "study_accession": "PRJL", "country_iso": "US", "est_depth_x": "99"})
    attrs.append({"biosample": "SAMNLAB", "isolation_source": "laboratory strain"})
    rows.append({"run_accession": "NOTEXT", "sample_accession": "SAMNNO", "study_accession": "PRJN", "country_iso": "", "est_depth_x": "40"})
    return rows, attrs


class Select(unittest.TestCase):
    def setUp(self):
        self._d = tempfile.TemporaryDirectory()
        self.t = self._d.name

    def tearDown(self):
        self._d.cleanup()

    def run_sel(self, rows, attrs, *argv):
        cand, at = make(self.t, rows, attrs)
        out = os.path.join(self.t, "sel.csv")
        with redirect_stderr(io.StringIO()):
            rc = ps.main([cand, "--biosamples", at, "--out", out, *argv])
        return rc, list(csv.DictReader(open(out)))

    def test_equal_shares_per_category_and_one_run_per_sample_and_lab_left_out(self):
        rows, attrs = fixture()
        rows.append(dict(rows[0], run_accession="RUN0001B", est_depth_x="5"))  # a second, shallower run of the same sample
        rc, sel = self.run_sel(rows, attrs, "-n", "40", "--by", "cultivation")
        self.assertEqual(rc, 0)
        self.assertEqual(len(sel), 40)
        self.assertNotIn("LAB1", {r["run_accession"] for r in sel})
        self.assertNotIn("RUN0001B", {r["run_accession"] for r in sel})
        from collections import Counter
        c = Counter(r["category"] for r in sel)
        # wine, beer, clinical, wild_nature and unannotated (1 run): the 39 seats the small category cannot use are shared
        self.assertEqual(c["unannotated"], 1)
        full = [c[k] for k in ("wine", "beer", "clinical", "wild_nature")]
        self.assertEqual(sum(full), 39)
        self.assertTrue(max(full) - min(full) <= 1, full)

    def test_the_same_seed_gives_the_same_set_and_another_seed_another(self):
        rows, attrs = fixture()
        a = self.run_sel(rows, attrs, "-n", "40", "--seed", "1")[1]
        b = self.run_sel(rows, attrs, "-n", "40", "--seed", "1")[1]
        c = self.run_sel(rows, attrs, "-n", "40", "--seed", "2")[1]
        self.assertEqual([r["run_accession"] for r in a], [r["run_accession"] for r in b])
        self.assertNotEqual([r["run_accession"] for r in a], [r["run_accession"] for r in c])

    def test_geography_is_spread_before_it_is_repeated(self):
        rows, attrs = fixture()
        rc, sel = self.run_sel(rows, attrs, "-n", "6", "--by", "cultivation", "--quota", "wine=6")
        self.assertEqual(rc, 0)
        self.assertEqual(len({r["country_iso"] for r in sel}), 6)  # six countries, one run each

    def test_quotas_and_the_cap_per_study(self):
        rows, attrs = fixture()
        rc, sel = self.run_sel(rows, attrs, "-n", "20", "--by", "cultivation", "--quota", "wine=15", "--quota", "beer=5", "--max-per-study", "4")
        from collections import Counter
        c = Counter(r["category"] for r in sel)
        self.assertEqual((c["wine"], c["beer"]), (15, 5))  # 4 wine studies x 4 = 16 >= 15
        self.assertTrue(max(Counter(r["study_accession"] for r in sel).values()) <= 4)

    def test_a_quota_a_small_category_cannot_fill_is_given_to_the_others(self):
        rows, attrs = fixture(n_per=30)
        rows = rows[:62]  # wine 30, beer 30, clinical 2
        rc, sel = self.run_sel(rows, attrs, "-n", "60", "--by", "cultivation", "--quota", "wine=20", "--quota", "beer=20", "--quota", "clinical=20")
        from collections import Counter
        c = Counter(r["category"] for r in sel)
        self.assertEqual((rc, len(sel)), (0, 60))
        self.assertEqual(c["clinical"], 2)
        self.assertEqual(c["wine"] + c["beer"], 58)

    def test_a_quota_above_what_exists_is_capped_and_reported_by_the_exit_code(self):
        rows, attrs = fixture(n_per=3)
        rc, sel = self.run_sel(rows, attrs, "-n", "50", "--by", "cultivation")
        self.assertEqual(rc, 3)  # fewer than asked for
        self.assertEqual(len(sel), len({r["run_accession"] for r in sel}))

    def test_keeping_lab_samples_is_an_option(self):
        rows, attrs = fixture(n_per=1)
        _, sel = self.run_sel(rows, attrs, "-n", "50", "--lab", "keep")
        self.assertIn("LAB1", {r["run_accession"] for r in sel})

    def test_a_sample_with_no_text_is_unannotated_and_gets_a_region_of_unknown(self):
        rows, attrs = fixture(n_per=1)
        _, sel = self.run_sel(rows, attrs, "-n", "50", "--by", "cultivation")
        r = next(x for x in sel if x["run_accession"] == "NOTEXT")
        self.assertEqual((r["category"], r["region"]), ("unannotated", "unknown_region"))

    def test_place_spreads_a_one_country_species_over_its_states(self):
        rows, attrs = [], []
        for i in range(24):
            rows.append({"run_accession": f"R{i}", "sample_accession": f"S{i}", "study_accession": "P1", "country_iso": "US", "est_depth_x": "40"})
            attrs.append({"biosample": f"S{i}", "isolation_source": "clinical", "geo_loc_name": ["USA: Arizona", "USA: California", "USA: Texas"][i % 3]})
        # without --place all runs are "US" and the pick is by hash; with it each state gets a third
        _, sel = self.run_sel(rows, attrs, "-n", "9", "--by", "species", "--place")
        from collections import Counter
        self.assertEqual(sorted(Counter(r["geo_loc_name"] for r in sel).values()), [3, 3, 3])

    def test_among_equal_places_the_shallower_run_is_taken(self):
        rows = [{"run_accession": f"R{i}", "sample_accession": f"S{i}", "study_accession": f"P{i}", "country_iso": "US", "est_depth_x": str(d)}
                for i, d in enumerate([300, 28, 160, 90])]
        attrs = [{"biosample": f"S{i}", "isolation_source": "clinical"} for i in range(4)]
        _, sel = self.run_sel(rows, attrs, "-n", "1", "--by", "species")
        self.assertEqual(sel[0]["run_accession"], "R1")  # 28x, the shallowest: the least to download

    def test_by_region_uses_the_continent(self):
        rows, attrs = fixture()
        rc, sel = self.run_sel(rows, attrs, "-n", "30", "--by", "region")
        from collections import Counter
        c = Counter(r["category"] for r in sel)
        self.assertTrue(max(c.values()) - min(c.values()) <= 6, c)  # about equal shares of the continents present


if __name__ == "__main__":
    unittest.main(verbosity=1)
