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

## Releases (GitHub release assets)

Indexes are data, released separately from code (`v*` tags are code): one GitHub release per index version, tag `index-YYYY.MM`, not marked as the latest release. First release: [index-2026.10](https://github.com/stajichlab/GmG/releases/tag/index-2026.10). Assets:

* `bfd-YYYY.MM-s10000.gmgi`: full reference, ~200 MB (23.7k genomes, ANI clusters at 0.95 in the `cluster` column);
* `bfd-reps-YYYY.MM-s20000.gmgi`: representatives only (best genome per species), for browsers and laptops, ~60 MB; fewer species calls (see the release notes);
* `SHA256SUMS` and `MANIFEST.json` (the `gmg info` output of each file).

`scripts/make_index_release.sh INDEX...` writes `SHA256SUMS` and `MANIFEST.txt` next to the files. Release assets are limited to 2 GB each, which these indexes are far below.

Building a release: `scripts/hpcc/15_release_index.sbatch` (HPCC) writes the two indexes, `SHA256SUMS` and `MANIFEST.txt`; `gh release create index-YYYY.MM <files> --latest=false` publishes them.
