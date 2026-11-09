"""scripts/make_species_complexes.py on a tiny records table."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load():
    spec = importlib.util.spec_from_file_location("msc", ROOT / "scripts" / "make_species_complexes.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_clean_complex():
    m = load()
    assert m.clean_complex("F. oxysporum species complex (FOSC)") == ("Fusarium oxysporum species complex (FOSC)", "FOSC")
    assert m.clean_complex("F. fujikuroi species complex (FFSC")[0] == "Fusarium fujikuroi species complex (FFSC)"
    assert m.clean_complex("Neocosmospora (previously the F. solani species complex")[0] == (
        "Neocosmospora (formerly Fusarium solani species complex)")
    assert m.clean_complex("Fusarium not assigned") == (None, None)


def test_build():
    m = load()
    rows = [
        dict(species="Fusarium oxysporum", complex="F. oxysporum species complex (FOSC)"),
        dict(species="Fusarium nirenbergiae", complex="F. oxysporum species complex (FOSC)"),
        dict(species="Fusarium graminearum", complex="F. sambucinum species complex (FSAMSC) including ..."),
        dict(species="Neocosmospora solani", complex="Neocosmospora (previously the F. solani species complex"),
        dict(species="Atractium stilbaster", complex="Atractium"),
        dict(species="Fusarium sp. FTSC 18 (Laraba et al. 2022)", complex="F. tricinctum species complex (FTSC)"),
    ]
    t = m.build(rows, {"FOSC"})
    assert t["Fusarium oxysporum"][1] == 1 and t["Fusarium nirenbergiae"][1] == 1
    assert t["Fusarium graminearum"] == ("Fusarium sambucinum species complex (FSAMSC)", 0)
    assert t["Fusarium solani"] == t["Neocosmospora solani"]  # NCBI-style alias
    assert "Atractium stilbaster" not in t  # group is just the genus
    assert not any(" sp." in k for k in t)
