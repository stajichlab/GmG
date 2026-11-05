//! Sketch one sample (the files of one sample merged): hashes with abundances for reads, plain hashes for assemblies, an
//! optional pair of half-sketches (alternate reads) for the split-half noise floor, and a small on-disk cache.

use anyhow::{bail, Context, Result};
use gmg_core::{FracSketch, SketchFile, DEFAULT_SEED};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq)]
pub struct SketchParams {
    pub ksize: u32,
    pub scaled: u32,
    /// bases read per sample (the files of a sample share it); 0 = no limit. Reads only: assemblies are read whole.
    pub max_bases: u64,
    /// read only the first file of each sample (R1 or the first lane)
    pub r1_only: bool,
}

impl SketchParams {
    fn tag(&self) -> String {
        format!(
            "k={};scaled={};max_bases={};r1_only={}",
            self.ksize, self.scaled, self.max_bases, self.r1_only
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
                if idx % 2 == 0 {
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

// ---- cache: out/sketches/{id}.gsk, {id}.half_a.gsk, {id}.half_b.gsk, {id}.json, params.txt ---------------------------------

fn gsk_path(dir: &Path, id: &str, part: &str) -> PathBuf {
    dir.join(format!("{id}{part}.gsk"))
}

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

fn read_spectrum(path: &Path) -> Result<Spectrum> {
    let mut r = BufReader::new(File::open(path)?);
    let sf = SketchFile::read(&mut r)?;
    let abunds = sf.abunds.unwrap_or_else(|| vec![1; sf.hashes.len()]);
    Ok(Spectrum {
        hashes: sf.hashes,
        abunds,
    })
}

pub fn cache_ready(dir: &Path, p: &SketchParams) -> bool {
    fs::read_to_string(dir.join("params.txt")).is_ok_and(|t| t.trim() == p.tag())
}

pub fn cache_init(dir: &Path, p: &SketchParams) -> Result<()> {
    fs::create_dir_all(dir)?;
    fs::write(dir.join("params.txt"), p.tag())?;
    Ok(())
}

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
        };
        let s = sketch_sample("g", "assembly", &[f], &p, true).unwrap();
        assert!(s.halves.is_none());
        assert!(s.spectrum.abunds.iter().all(|&a| a == 1));
        assert_eq!(s.counts.bases, seq.len() as u64); // the limit does not apply to assemblies
        fs::remove_dir_all(&dir).ok();
    }
}
