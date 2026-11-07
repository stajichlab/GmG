# Plan: compact sketches and hosted sketch sets for the popsketch web page

Status: draft v2, 2026-11-06, revised after an Opus review checked against the code (nothing built or measured). Numbers marked (est.) are estimates. v1 had the size about 3 times too low and missed the index screen; both are fixed below.

## Goal
Let a person open `popsketch.html`, pick a species or a set of runs (SRA accessions), and see the popsketch analysis (QC flags, PCA/PCoA, clusters, tree, label overlays) without sketching anything in the browser and without downloading reads. The sketches are built ahead of time on HPCC and served as static files.

## What exists
- 2,002 sketches of 1 Gb of bases each (7 fungal species), about 8 MB per sample on disk: whole `.gsk` about 3 MB, two split halves about 2 MB each, plus `<id>.json`. 24 GB in all.
- `.gsk` (gmg-core `skfile.rs`): magic, k, scaled, seed, name, `n` hashes as raw u64, optional abundances as raw u32.
- Web page opens a bundle (`.tar.gz` of the above, `scripts/popsketch_bundle.py`, `web/bundle.js`) from a file or `?data=URL`, then `Session.add_sketch` + `analyze()` in a worker (wasm). The 692 MB *C. neoformans* bundle works but is too large to serve.
- `analyze()` uses per sample: the hashes with their abundances (abundance cutoff is chosen from the histogram in `coverage.rs`, so the kept hashes depend on the cutoff), the spectrum (histogram) for coverage/size/shoulder flags, the read counts and base counts in the `.json`, and, for 20 samples picked by `pick_for_halves`, the split halves for the noise floor.
- Estimated: about 175,000 of 320,000 hashes per 1 Gb sample survive the abundance cutoff; the rest are error k-mers seen once or twice.

## Design

### 1. Compact sample file (`.gsc`, format 1)
Per sample, one file holding what analysis needs:
- header: k, scaled, seed, id, read and base counts, the sketch parameters, the code version, the abundance cutoff used, the split-half distance `d` (section 2), and a checksum of the payload (a corrupt varint stream decodes to plausible hashes);
- the abundance histogram, stored to its true maximum (pooled at 4096): `assess` uses the abundances as a multiset and `hashes.len()`, and `valley_peak` depends on the histogram's length, so do not pad a fixed array;
- the kept hashes (abundance >= the cutoff), sorted, delta-coded. No per-hash abundances: kept_hashes, hash PCA and the decontaminate vote use presence only (qc.rs:302, pipeline.rs:812, :659).

Size, corrected: at scaled 200 the hashes lie below 2^64/200 = 9.2e16; with 175,000 kept the mean gap is 5.3e11 (2^38.9). Varint needs 5 to 6 bytes per gap; the entropy floor is about 5.0 bytes (Golomb-Rice or Elias-Fano get close). So about 0.9 to 1.0 MB per sample (0.88 MB at the entropy floor), not the 0.35 to 0.6 MB of v1. 2,002 samples are about 2 GB, 300 samples about 300 MB (est.), still 3 to 4 times smaller than the whole `.gsk`.

**Index screen.** `screen()` classifies against the reference index using every spectrum hash at or below the index's max hash, with a fixed abundance cutoff of 2 for reads (pipeline.rs:246-282), because the auto cutoff drops the low-abundance k-mers of a light contaminant. Compact input without those hashes silently loses `mixture_suspected` and `different_organism`. Choose before building:
- (a) a second layer: hashes with abundance >= 2 at or below the index's max hash (about scaled_sample / scaled_index, 200/10000 = 2% of the sample's hashes), stored with their abundance;
- (b) no screen on compact input, and the page says so.
Prefer (a) if it costs under 0.1 MB per sample (measure).

### 2. Noise floor: store the split-half distance per sample (option E)
`floor::estimate` needs one number per sample, the split-half distance `d` (floor.rs:46-64), and the floor is the 90th percentile of `d` over the picked samples that are not flagged (pipeline.rs:700-704). So each `.gsc` carries `d`, computed at packing time from the halves (cutoff applied then, so halves need no histogram of their own). `Session` still picks its 20 with `pick_for_halves` for the chosen set and takes their stored `d`. This is exact under any subset, fetches nothing extra and keeps the meaning of the floor. Rejected: a per-species floor (changes what the floor is), a fixed 20 samples (changes which samples count, and the count still varies when some are flagged), halves for everyone (adds about 2 MB per sample). Cost: `d` is computed only for samples whose halves exist; every one of the 2,002 has them.

### 3. Reader and parity
- `compact_io` in popsketch (writer and reader), `Session.add_compact` in wasm beside `add_sketch`, and a packer (`popsketch pack` or an extension of `scripts/popsketch_bundle.py`).
- Change in `SampleSketch`: a histogram field (`Spectrum` becomes an enum, or `hist: Option<Vec<u64>>`) and branches at `assess` (qc.rs:146-171), `kept_hashes` (qc.rs:302), `screen` (pipeline.rs:247) and `floor::estimate` (take stored `d`). Do not fake a spectrum with dummy hashes: `screen` would pick them up.
- The cutoff is stored with the code version. On read, recompute it from the histogram and fail on a mismatch, so a later change to `auto_min_abund` (#159) cannot leave stale files unnoticed.
- Parity tests: native `run --from-sketches` on the original sketches against the compact files, covering random subsets of different sizes (QC is cohort-relative: size_z, shoulder, `unrelated_to_cohort`, and the picked halves all depend on the subset), not only the seven whole species; plus wasm against native, as `parity.mjs` does. Expect exact identity of ani, flags and clusters; any difference blocks the plan until explained. Not yet checked: whether `dist::compute` is bit-identical between native and wasm.

### 4. Catalog
`catalog.json` per species: accession, biosample, study, species, country, strain, host, bases, file size, URL, plus sketch parameters. Built from `runs.tsv` and ENA metadata. Samples excluded for failure (the 21 *P. oryzae* runs that are not at NCBI or ENA; dropped by decision, 2026-11-07) are listed with the reason, not silently absent.

### 5. Page
- Species picker or accession paste; shows the catalog rows, total download size, a cap (for example 400 samples), a progress bar, and cancels cleanly.
- One whole file per sample (no Range), 6 in parallel, retry with backoff, SHA-256 checked for every file (cheap at about 1 MB), cached in the browser (Cache API) so a re-run does not refetch. Range is only needed if per-species packed blobs are added later; GCS ignores Range on gzip-transcoded objects.
- Runs `analyze()` in the existing worker; labels from the catalog columns.
- Falls back to the existing bundle and upload paths; nothing is removed.
- A visible provenance line: sketch parameters, build date, binary version, so a result can be tied to its inputs.

### 6. Hosting
- **Decided (2026-11-07): NRP Ceph S3**, existing bucket `s3://stajichlab`, prefix `sketch/fungi/`, uploaded with `s3cmd -P` (public read) using the owner's `~/.s3cfg` (laptop or HPCC). GCP or R2 stay the fallback. Public read bucket, CORS allowing GET/HEAD from `https://gmg.fungalgenomes.org` (and localhost for tests), exposing `Content-Length` and `ETag`. Serve the files uncompressed: delta-coded hashes barely gzip.
- Content-hashed file names with `Cache-Control: immutable`; short TTL on the catalog.
- One config value (`SKETCH_BASE_URL`) on the site; catalogs hold relative paths.
- GCP: egress is the cost (about $0.12/GB, so about $0.05 for a 400-sample run); a budget alert does not cap spend, so consider R2 (no egress fees) or a CDN with rate limits. NRP: unknown limits and end of the free period; the sketches are rebuildable but URLs in docs and bundles would break, so keep the base URL in one place.
- Versioned prefixes (`v1/…`) so the format can change without breaking old links.
- Integrity: a manifest with SHA-256 for every file; the page verifies what it fetches for the catalog, not necessarily every sample.

### 7. HPCC work
- Packing the existing 2,002: CPU only, no downloads, a few minutes to an hour (est.). Run as an array with a low priority and `ArrayTaskThrottle`, outside the queue crowding of other jobs.
- Credentials: the owner's `~/.s3cfg` is on the laptop and on HPCC (owner's home only); run `s3cmd` from either. Do not copy it into shared folders, job scripts or CI.
- Re-sketching is not needed. Measurement and parity are the only compute.

## Steps and gates
1. Measure on one species (*C. neoformans*): write the kept-hash payload, report the real bytes per sample for hashes, histogram, `d`, and the index-screen layer. Gate: whole-sample file under about 1.5 MB. (Reads existing data only.)
2. Decide the index screen (a or b) from that measurement.
3. Implement `compact_io`, `add_compact` and the packer; parity on random subsets (native and wasm). Gate: exact, or each difference explained and accepted.
4. Pack all 2,002 on HPCC (CPU only, array, throttled); catalog; page against a local server with CORS headers; headless-browser e2e including a memory check at the intended cap.
5. Bucket, CORS, upload, real-page test, cost check.
6. Docs and release.
HPCC load is needed for steps 1 and 4 only, and is small.

## Not covered / risks
- Only the 7 species sketched so far. More species or all of SRA fungi (tens of thousands of runs) is a separate project (fetch time, storage).
- Sketches are 1 Gb subsamples at k=21, scaled 200; the depth findings in `docs/popsketch.md` apply. #159 (estimator validation) is open; the abundance cutoff in compact files inherits any remaining estimator errors (two known runs).
- Public data only; the catalog carries SRA metadata, nothing private.
- Browser memory: the hash arrays are about 560 MB for 400 samples, but `analyze` copies the kept sets three times (pipeline.rs:595, :615, :769) plus the gram matrix map, so the peak may reach about 1.5 GB plus fetched buffers (reviewer estimate, not measured); wasm32 is limited to 4 GB. Measure, and decide the cap (the CLI's `max_samples` is 1500, so say what the page does for larger sets).
- A format change means the old bundle path must keep working or be retired explicitly.
- Anyone can pull from a public bucket: costs and abuse limits (budget alerts, rate limits).
- Ownership: who maintains the bucket after the free period; where the packing scripts live (this repo).

## Alternatives considered (from review)
- Precomputed per-species results for the 'whole species' case: no per-sample download, but no subsets or user labels.
- An ANI matrix plus per-sample scalars (cutoff, size, coverage, shoulder pair, `d`) reproduces everything except hash PCA (the default), decontaminate and the screen.
- Scaled 1000: about 0.2 MB per sample, but changes cutoffs and distance resolution and breaks parity with results so far.
- Serving the existing gzipped `.gsk`: about 2 to 3 MB per sample, no new code.
