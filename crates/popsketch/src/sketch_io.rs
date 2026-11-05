//! Sketch one sample (the files of one sample merged): hashes with abundances for reads, plain hashes for assemblies, an
//! optional pair of half-sketches (alternate reads) for the split-half noise floor, and a small on-disk cache.

#[cfg(feature = "native")]
use anyhow::Context;
use anyhow::{bail, Result};
#[cfg(feature = "native")]
use gmg_core::coverage::profile;
#[cfg(feature = "native")]
use gmg_core::SketchFile;
use gmg_core::{FracSketch, DEFAULT_SEED};
#[cfg(feature = "native")]
use std::fs::{self, File};
#[cfg(feature = "native")]
use std::io::{BufReader, BufWriter};
#[cfg(feature = "native")]
use std::path::Path;
use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq)]
pub struct SketchParams {
    pub ksize: u32,
    pub scaled: u32,
    /// bases read per sample (the files of a sample share it); 0 = no limit. Reads only: assemblies are read whole.
    pub max_bases: u64,
    /// read only the first file of each sample (R1 or the first lane)
    pub r1_only: bool,
    /// reads: stop once the k-mer coverage estimated from the reads so far reaches this (checked every `chunk_reads` reads);
    /// `max_bases` still applies as a ceiling
    pub target_cov: Option<f64>,
    /// reads between coverage checks when `target_cov` is set
    pub chunk_reads: usize,
}

impl Default for SketchParams {
    fn default() -> Self {
        SketchParams {
            ksize: 21,
            scaled: 200,
            max_bases: 1_000_000_000,
            r1_only: false,
            target_cov: None,
            chunk_reads: 200_000,
        }
    }
}

impl SketchParams {
    #[cfg(feature = "native")]
    fn tag(&self) -> String {
        format!(
            "k={};scaled={};max_bases={};r1_only={};target_cov={}",
            self.ksize,
            self.scaled,
            self.max_bases,
            self.r1_only,
            self.target_cov
                .map_or("none".to_string(), |c| format!("{c}"))
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Counts {
    pub reads: u64,
    pub bases: u64,
    pub non_acgt: u64,
}

/// Hashes (sorted, unique) and how often each was seen.
#[derive(Clone, Debug, Default)]
pub struct Spectrum {
    pub hashes: Vec<u64>,
    pub abunds: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct SampleSketch {
    pub id: String,
    /// `reads` or `assembly`
    pub kind: String,
    pub files: Vec<PathBuf>,
    pub spectrum: Spectrum,
    pub counts: Counts,
    /// two sketches of alternate reads of the same sample (reads only, a subset of samples)
    pub halves: Option<(Spectrum, Spectrum)>,
}

fn finish(sk: FracSketch, reads: bool) -> Spectrum {
    if reads {
        let (hashes, abunds) = sk.finish_with_abundance();
        Spectrum { hashes, abunds }
    } else {
        let hashes = sk.finish();
        let abunds = vec![1; hashes.len()];
        Spectrum { hashes, abunds }
    }
}

/// Sketch a sample held in memory: `parts` are the (already decompressed) FASTQ or FASTA contents of its files. The same
/// counting as `sketch_sample` without a coverage target: reads up to `max_bases` shared by the parts, optional split halves.
pub fn sketch_bytes(
    id: &str,
    kind: &str,
    parts: &[&[u8]],
    p: &SketchParams,
    want_halves: bool,
) -> Result<SampleSketch> {
    let reads = kind == "reads";
    let used: &[&[u8]] = if reads && p.r1_only {
        &parts[..1.min(parts.len())]
    } else {
        parts
    };
    if used.is_empty() {
        bail!("{id}: no files");
    }
    let per_file = if reads && p.max_bases > 0 {
        p.max_bases.div_ceil(used.len() as u64)
    } else {
        u64::MAX
    };
    let k = p.ksize as usize;
    let mut main = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let halves_wanted = reads && want_halves;
    let mut ha = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let mut hb = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let mut counts = Counts::default();
    for bytes in used {
        let (mut bases, mut idx) = (0u64, 0u64);
        gmg_core::for_each_record(bytes, &mut |seq: &[u8]| {
            main.add_sequence(seq);
            if halves_wanted {
                if idx & 1 == 0 {
                    ha.add_sequence(seq);
                } else {
                    hb.add_sequence(seq);
                }
            }
            idx += 1;
            counts.reads += 1;
            counts.bases += seq.len() as u64;
            counts.non_acgt += seq
                .iter()
                .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                .count() as u64;
            bases += seq.len() as u64;
            bases < per_file
        });
    }
    if counts.reads == 0 {
        bail!("{id}: no sequence records");
    }
    Ok(SampleSketch {
        id: id.to_string(),
        kind: kind.to_string(),
        files: Vec::new(),
        spectrum: finish(main, reads),
        counts,
        halves: halves_wanted.then(|| (finish(ha, true), finish(hb, true))),
    })
}

#[cfg(feature = "native")]
/// Sketch the files of one sample. For reads the base limit is shared by the files (each gets an equal part), and with
/// `want_halves` the reads of each file are also dealt alternately into two further sketches.
pub fn sketch_sample(
    id: &str,
    kind: &str,
    files: &[PathBuf],
    p: &SketchParams,
    want_halves: bool,
) -> Result<SampleSketch> {
    let reads = kind == "reads";
    let used: &[PathBuf] = if reads && p.r1_only {
        &files[..1.min(files.len())]
    } else {
        files
    };
    if used.is_empty() {
        bail!("{id}: no files");
    }
    let per_file = if reads && p.max_bases > 0 {
        p.max_bases.div_ceil(used.len() as u64)
    } else {
        u64::MAX
    };
    if reads && p.target_cov.is_some() {
        return sketch_reads_chunked(id, kind, used, p, per_file, want_halves);
    }
    let k = p.ksize as usize;
    let mut main = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let halves_wanted = reads && want_halves;
    let mut ha = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let mut hb = FracSketch::new(k, p.scaled, DEFAULT_SEED);
    let mut counts = Counts::default();
    for f in used {
        let mut rd = needletail::parse_fastx_file(f)
            .with_context(|| format!("{id}: opening {}", f.display()))?;
        let (mut bases, mut idx) = (0u64, 0u64);
        while let Some(rec) = rd.next() {
            let rec = rec.with_context(|| format!("{id}: reading {}", f.display()))?;
            let seq = rec.seq();
            main.add_sequence(&seq);
            if halves_wanted {
                if idx & 1 == 0 {
                    ha.add_sequence(&seq);
                } else {
                    hb.add_sequence(&seq);
                }
            }
            idx += 1;
            counts.reads += 1;
            counts.bases += seq.len() as u64;
            counts.non_acgt += seq
                .iter()
                .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                .count() as u64;
            bases += seq.len() as u64;
            if bases >= per_file {
                break;
            }
        }
    }
    if counts.reads == 0 {
        bail!("{id}: no sequence records");
    }
    Ok(SampleSketch {
        id: id.to_string(),
        kind: kind.to_string(),
        files: used.to_vec(),
        spectrum: finish(main, reads),
        counts,
        halves: halves_wanted.then(|| (finish(ha, true), finish(hb, true))),
    })
}

#[cfg(feature = "native")]
/// Merge a chunk's sorted (hash, count) lists into the running sorted lists.
fn merge_counts(h: &mut Vec<u64>, a: &mut Vec<u32>, nh: Vec<u64>, na: Vec<u32>) {
    let cap = h.len() + nh.len();
    let (mut oh, mut oa) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
    let (mut i, mut j) = (0, 0);
    while i < h.len() && j < nh.len() {
        match h[i].cmp(&nh[j]) {
            std::cmp::Ordering::Less => {
                oh.push(h[i]);
                oa.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                oh.push(nh[j]);
                oa.push(na[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                oh.push(h[i]);
                oa.push(a[i] + na[j]);
                i += 1;
                j += 1;
            }
        }
    }
    oh.extend_from_slice(&h[i..]);
    oa.extend_from_slice(&a[i..]);
    oh.extend_from_slice(&nh[j..]);
    oa.extend_from_slice(&na[j..]);
    *h = oh;
    *a = oa;
}

#[cfg(feature = "native")]
/// Reads sketched in chunks, with the files of the sample read in lockstep (R1 and R2 alike), stopping once the k-mer coverage
/// estimated from the reads so far reaches `target_cov`, or at the per-file base ceiling.
fn sketch_reads_chunked(
    id: &str,
    kind: &str,
    used: &[PathBuf],
    p: &SketchParams,
    per_file: u64,
    want_halves: bool,
) -> Result<SampleSketch> {
    let target = p.target_cov.unwrap_or(f64::INFINITY);
    let k = p.ksize as usize;
    let mut readers = Vec::new();
    for f in used {
        readers.push(Some(
            needletail::parse_fastx_file(f)
                .with_context(|| format!("{id}: opening {}", f.display()))?,
        ));
    }
    let new_sketches = || {
        (
            FracSketch::new(k, p.scaled, DEFAULT_SEED),
            FracSketch::new(k, p.scaled, DEFAULT_SEED),
            FracSketch::new(k, p.scaled, DEFAULT_SEED),
        )
    };
    let (mut cm, mut ca, mut cb) = new_sketches();
    let (mut mh, mut ma) = (Vec::new(), Vec::new());
    let (mut ah, mut aa) = (Vec::new(), Vec::new());
    let (mut bh, mut ba) = (Vec::new(), Vec::new());
    let mut counts = Counts::default();
    let mut file_bases = vec![0u64; used.len()];
    let mut idx = vec![0u64; used.len()];
    let mut in_chunk = 0usize;
    let mut flush =
        |cm: FracSketch, ca: FracSketch, cb: FracSketch, mh: &mut Vec<u64>, ma: &mut Vec<u32>| {
            let (h, a) = cm.finish_with_abundance();
            merge_counts(mh, ma, h, a);
            if want_halves {
                let (h, a) = ca.finish_with_abundance();
                merge_counts(&mut ah, &mut aa, h, a);
                let (h, a) = cb.finish_with_abundance();
                merge_counts(&mut bh, &mut ba, h, a);
            }
        };
    loop {
        let mut progressed = false;
        for fi in 0..readers.len() {
            let Some(rd) = readers[fi].as_mut() else {
                continue;
            };
            if file_bases[fi] >= per_file {
                readers[fi] = None;
                continue;
            }
            let Some(rec) = rd.next() else {
                readers[fi] = None;
                continue;
            };
            let rec = rec.with_context(|| format!("{id}: reading {}", used[fi].display()))?;
            let seq = rec.seq();
            cm.add_sequence(&seq);
            if want_halves {
                if idx[fi] & 1 == 0 {
                    ca.add_sequence(&seq);
                } else {
                    cb.add_sequence(&seq);
                }
            }
            idx[fi] += 1;
            counts.reads += 1;
            counts.bases += seq.len() as u64;
            counts.non_acgt += seq
                .iter()
                .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                .count() as u64;
            file_bases[fi] += seq.len() as u64;
            in_chunk += 1;
            progressed = true;
        }
        if in_chunk >= p.chunk_reads {
            let fresh = new_sketches();
            let done = (
                std::mem::replace(&mut cm, fresh.0),
                std::mem::replace(&mut ca, fresh.1),
                std::mem::replace(&mut cb, fresh.2),
            );
            flush(done.0, done.1, done.2, &mut mh, &mut ma);
            in_chunk = 0;
            if profile(&ma).is_some_and(|c| c.kmer >= target) {
                break;
            }
        }
        if !progressed {
            break;
        }
    }
    flush(cm, ca, cb, &mut mh, &mut ma);
    if counts.reads == 0 {
        bail!("{id}: no sequence records");
    }
    Ok(SampleSketch {
        id: id.to_string(),
        kind: kind.to_string(),
        files: used.to_vec(),
        spectrum: Spectrum {
            hashes: mh,
            abunds: ma,
        },
        counts,
        halves: want_halves.then_some((
            Spectrum {
                hashes: ah,
                abunds: aa,
            },
            Spectrum {
                hashes: bh,
                abunds: ba,
            },
        )),
    })
}

// ---- cache: out/sketches/{id}.gsk, {id}.half_a.gsk, {id}.half_b.gsk, {id}.json, params.txt ---------------------------------

#[cfg(feature = "native")]
fn gsk_path(dir: &Path, id: &str, part: &str) -> PathBuf {
    dir.join(format!("{id}{part}.gsk"))
}

#[cfg(feature = "native")]
fn write_spectrum(path: &Path, p: &SketchParams, name: &str, s: &Spectrum) -> Result<()> {
    let sf = SketchFile {
        ksize: p.ksize,
        scaled: p.scaled,
        seed: DEFAULT_SEED,
        name: name.to_string(),
        hashes: s.hashes.clone(),
        abunds: Some(s.abunds.clone()),
    };
    let mut w = BufWriter::new(File::create(path).with_context(|| path.display().to_string())?);
    sf.write(&mut w)?;
    Ok(())
}

#[cfg(feature = "native")]
fn read_spectrum(path: &Path) -> Result<Spectrum> {
    let mut r = BufReader::new(File::open(path)?);
    let sf = SketchFile::read(&mut r)?;
    let abunds = sf.abunds.unwrap_or_else(|| vec![1; sf.hashes.len()]);
    Ok(Spectrum {
        hashes: sf.hashes,
        abunds,
    })
}

#[cfg(feature = "native")]
pub fn cache_ready(dir: &Path, p: &SketchParams) -> bool {
    fs::read_to_string(dir.join("params.txt")).is_ok_and(|t| t.trim() == p.tag())
}

#[cfg(feature = "native")]
pub fn cache_init(dir: &Path, p: &SketchParams) -> Result<()> {
    fs::create_dir_all(dir)?;
    fs::write(dir.join("params.txt"), p.tag())?;
    Ok(())
}

#[cfg(feature = "native")]
pub fn cache_write(dir: &Path, s: &SampleSketch, p: &SketchParams) -> Result<()> {
    write_spectrum(&gsk_path(dir, &s.id, ""), p, &s.id, &s.spectrum)?;
    if let Some((a, b)) = &s.halves {
        write_spectrum(&gsk_path(dir, &s.id, ".half_a"), p, &s.id, a)?;
        write_spectrum(&gsk_path(dir, &s.id, ".half_b"), p, &s.id, b)?;
    }
    let files: Vec<String> = s.files.iter().map(|f| f.display().to_string()).collect();
    let meta = serde_json::json!({"id": s.id, "kind": s.kind, "reads": s.counts.reads, "bases": s.counts.bases,
                                  "non_acgt": s.counts.non_acgt, "files": files});
    fs::write(
        dir.join(format!("{}.json", s.id)),
        serde_json::to_string(&meta)?,
    )?;
    Ok(())
}

#[cfg(feature = "native")]
/// A cached sketch, if it is there and has the halves when they are wanted.
pub fn cache_read(dir: &Path, id: &str, want_halves: bool) -> Option<SampleSketch> {
    let meta: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join(format!("{id}.json"))).ok()?).ok()?;
    let spectrum = read_spectrum(&gsk_path(dir, id, "")).ok()?;
    let kind = meta["kind"].as_str()?.to_string();
    let halves = if want_halves && kind == "reads" {
        Some((
            read_spectrum(&gsk_path(dir, id, ".half_a")).ok()?,
            read_spectrum(&gsk_path(dir, id, ".half_b")).ok()?,
        ))
    } else {
        None
    };
    Some(SampleSketch {
        id: id.to_string(),
        kind,
        files: meta["files"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_str().map(PathBuf::from))
                    .collect()
            })
            .unwrap_or_default(),
        spectrum,
        counts: Counts {
            reads: meta["reads"].as_u64()?,
            bases: meta["bases"].as_u64()?,
            non_acgt: meta["non_acgt"].as_u64().unwrap_or(0),
        },
        halves,
    })
}

#[cfg(feature = "native")]
#[cfg(test)]
mod tests {
    use super::*;

    fn write_fastq(path: &Path, n: usize, len: usize) {
        let mut s = String::new();
        let mut x: u64 = 12345;
        for i in 0..n {
            let mut seq = String::new();
            for _ in 0..len {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                seq.push(b"ACGT"[(x % 4) as usize] as char);
            }
            s.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(len)));
        }
        fs::write(path, s).unwrap();
    }

    #[test]
    fn base_limit_is_shared_by_the_files_of_a_sample() {
        let dir = std::env::temp_dir().join(format!("popsketch_io_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (r1, r2) = (dir.join("a_R1.fastq"), dir.join("a_R2.fastq"));
        write_fastq(&r1, 2000, 100);
        write_fastq(&r2, 2000, 100);
        let p = SketchParams {
            ksize: 21,
            scaled: 10,
            max_bases: 100_000,
            r1_only: false,
            ..SketchParams::default()
        };
        let s = sketch_sample("a", "reads", &[r1.clone(), r2.clone()], &p, true).unwrap();
        // 100,000 bases split over two files = 50,000 bases = 500 reads each
        assert_eq!(s.counts.reads, 1000);
        assert_eq!(s.counts.bases, 100_000);
        let (a, b) = s.halves.as_ref().unwrap();
        assert!(!a.hashes.is_empty() && !b.hashes.is_empty());
        let solo = SketchParams {
            r1_only: true,
            ..p.clone()
        };
        let s1 = sketch_sample("a", "reads", &[r1, r2], &solo, false).unwrap();
        assert_eq!(s1.counts.reads, 1000); // one file takes the whole limit
        assert_eq!(s1.files.len(), 1);
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn cache_round_trips_and_checks_parameters() {
        let dir = std::env::temp_dir().join(format!("popsketch_cache_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let f = dir.join("s.fastq");
        write_fastq(&f, 800, 100);
        let p = SketchParams {
            ksize: 21,
            scaled: 5,
            max_bases: 0,
            r1_only: false,
            ..SketchParams::default()
        };
        let s = sketch_sample("s", "reads", &[f], &p, true).unwrap();
        let cd = dir.join("sketches");
        cache_init(&cd, &p).unwrap();
        cache_write(&cd, &s, &p).unwrap();
        assert!(cache_ready(&cd, &p));
        assert!(!cache_ready(
            &cd,
            &SketchParams {
                scaled: 6,
                ..p.clone()
            }
        ));
        let back = cache_read(&cd, "s", true).unwrap();
        assert_eq!(back.spectrum.hashes, s.spectrum.hashes);
        assert_eq!(back.spectrum.abunds, s.spectrum.abunds);
        assert_eq!(back.counts, s.counts);
        assert!(back.halves.is_some());
        assert!(cache_read(&cd, "missing", false).is_none());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn assemblies_have_unit_abundance_and_no_halves() {
        let dir = std::env::temp_dir().join(format!("popsketch_asm_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let f = dir.join("g.fa");
        let seq = "ACGTACGGTCAGTCAGGATCGATCGGATCGATTTAGCGCGATATCGGCTAGCTAGGATCGA";
        fs::write(&f, format!(">c1\n{seq}\n")).unwrap();
        let p = SketchParams {
            ksize: 11,
            scaled: 1,
            max_bases: 5,
            r1_only: true,
            ..SketchParams::default()
        };
        let s = sketch_sample("g", "assembly", &[f], &p, true).unwrap();
        assert!(s.halves.is_none());
        assert!(s.spectrum.abunds.iter().all(|&a| a == 1));
        assert_eq!(s.counts.bases, seq.len() as u64); // the limit does not apply to assemblies
        fs::remove_dir_all(&dir).ok();
    }

    /// FASTQ of `n` reads of 150 bases from random positions of a random genome of `len` bases (both files get different reads)
    fn write_reads(path: &Path, genome: &[u8], n: usize, seed: u64) {
        let mut x = seed;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut s = String::new();
        for i in 0..n {
            let start = (next() % (genome.len() - 150) as u64) as usize;
            s.push_str(&format!(
                "@r{i}\n{}\n+\n{}\n",
                String::from_utf8_lossy(&genome[start..start + 150]),
                "I".repeat(150)
            ));
        }
        fs::write(path, s).unwrap();
    }

    fn genome(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn target_coverage_stops_reading_early_and_is_otherwise_the_same_sketch() {
        let dir = std::env::temp_dir().join(format!("popsketch_target_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let g = genome(60_000, 5);
        let (r1, r2) = (dir.join("t_R1.fastq"), dir.join("t_R2.fastq"));
        write_reads(&r1, &g, 6_000, 11); // two files of 6,000 reads: 30x in total
        write_reads(&r2, &g, 6_000, 12);
        let files = [r1, r2];
        let plain = SketchParams {
            scaled: 1,
            max_bases: 0,
            ..SketchParams::default()
        };
        let all = sketch_sample("t", "reads", &files, &plain, false).unwrap();
        assert_eq!(all.counts.reads, 12_000);
        // a target that is never reached: the chunked path gives exactly the plain sketch
        let never = SketchParams {
            target_cov: Some(1e9),
            chunk_reads: 500,
            ..plain.clone()
        };
        let chunked = sketch_sample("t", "reads", &files, &never, true).unwrap();
        assert_eq!(chunked.counts, all.counts);
        assert_eq!(chunked.spectrum.hashes, all.spectrum.hashes);
        assert_eq!(chunked.spectrum.abunds, all.spectrum.abunds);
        let (ha, hb) = chunked.halves.as_ref().unwrap();
        assert!(ha.hashes.len() > 1000 && hb.hashes.len() > 1000);
        // a reachable target stops early, in whole rounds over the files, at about that coverage
        let early = SketchParams {
            target_cov: Some(8.0),
            chunk_reads: 500,
            ..plain.clone()
        };
        let s = sketch_sample("t", "reads", &files, &early, false).unwrap();
        assert!(
            s.counts.reads < 8_000 && s.counts.reads > 2_000,
            "reads {}",
            s.counts.reads
        );
        assert_eq!(s.counts.reads % 2, 0, "R1 and R2 are read in lockstep");
        let cov = gmg_core::coverage::profile(&s.spectrum.abunds)
            .unwrap()
            .kmer;
        assert!((8.0..=12.0).contains(&cov), "coverage {cov}");
        // the base ceiling still applies
        let capped = SketchParams {
            max_bases: 150 * 2_000,
            ..early.clone()
        };
        let c = sketch_sample("t", "reads", &files, &capped, false).unwrap();
        assert!(c.counts.reads <= 2_000, "reads {}", c.counts.reads);
        // the cache knows the setting
        let cd = dir.join("sk");
        cache_init(&cd, &early).unwrap();
        assert!(cache_ready(&cd, &early) && !cache_ready(&cd, &plain));
        fs::remove_dir_all(&dir).ok();
    }
}
