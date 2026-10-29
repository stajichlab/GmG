# Reference indexes: provenance and distribution

An index (`.gmgi`) carries a one-line JSON **manifest** written by `gmg build` (`gmg info INDEX` prints it as `manifest.*`; `gmg run` copies it into `run.json` under `index.manifest`):

| key | meaning |
|---|---|
| `gmg`, `commit` | tool version and git commit that built it |
| `built_utc` | build time |
| `format` | `GMGIDX01` (raw) or `GMGIDX02` (compact, after `gmg convert`); `converted_utc` is added by convert |
| `k`, `scaled`, `seed` | sketch parameters; queries with a different k or seed are refused, a coarser query `scaled` is refused, a finer one is downsampled |
| `genomes`, `skipped` | genomes indexed, and metadata rows without a usable sketch |
| `meta_sha256`, `meta_name` | checksum and name of the metadata table the index was built from |
| `taxonomy_version`, `note` | set with `gmg build --taxonomy-version "NCBI 2026-10" --note "..."` |

Indexes built before manifests existed print `manifest (none ...)`; rebuild (`gmg build`, minutes on the HPCC) to get one.

## Proposed release layout (needs a decision on hosting)

Indexes are data, released separately from code (`v*` tags are code). Suggested: a GitHub release per index version with tag `index-YYYY.MM`, assets:

* `bfd-s10000.v2.gmgi`: full reference, ~200 MB (23.7k genomes, ANI clusters included once `gmg cluster` output is merged into the metadata);
* `bfd-reps-s20000.v2.gmgi`: representatives only, for browsers and laptops (target <= 100 MB);
* `SHA256SUMS` and `MANIFEST.json` (the `gmg info` output of each file).

`scripts/make_index_release.sh INDEX...` writes `SHA256SUMS` and `MANIFEST.txt` next to the files. Release assets are limited to 2 GB each, which these indexes are far below.
