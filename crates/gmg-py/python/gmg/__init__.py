"""Guess My Genome: fungal taxonomy from FracMinHash sketches.

    import gmg
    ix = gmg.Index("bfd-reps-2026.10.5-s20000.gmgi")
    rec = ix.classify("genome.fa.gz")          # the record `gmg run` writes to results.jsonl
    rec["status"], rec["guess"]["calls"]["genus"]["taxon"]
"""
from ._gmg import Index, SCHEMA_VERSION, __version__, hits_tsv_text, results_tsv_text

__all__ = ["Index", "SCHEMA_VERSION", "__version__", "results_tsv_text", "hits_tsv_text"]
