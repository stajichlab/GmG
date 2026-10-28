//! Reference index: genomes with a 6-rank lineage and an inverted hash -> genomes map.
//!
//! In memory and on disk the map is three parallel arrays: sorted unique `hashes`,
//! `offsets` (len = hashes + 1) and `postings` (genome ids, ascending within a hash).
//! A search of one query sketch is one binary search per query hash.

use crate::skfile::{rd_u32, rd_u32_vec, rd_u64, rd_u64_vec, wr_u32_slice, wr_u64_slice};
use std::collections::HashMap;
use std::io::{self, Read, Write};

pub const RANKS: [&str; 6] = ["phylum", "class", "order", "family", "genus", "species"];
const MAGIC: &[u8; 8] = b"GMGIDX01";
const MAGIC_V1: &[u8; 8] = b"GMGIDX02"; // delta + LEB128 varint coded postings

#[derive(Clone, Debug, PartialEq)]
pub struct GenomeMeta {
    pub id: String,
    pub taxid: u32,
    /// phylum, class, order, family, genus, species; empty string = unknown
    pub lineage: [String; 6],
    /// number of retained hashes of this genome (filled by the builder)
    pub n_hashes: u32,
}

pub struct RefIndex {
    pub ksize: u32,
    pub scaled: u32,
    pub seed: u64,
    pub genomes: Vec<GenomeMeta>,
    /// per genome, per rank: interned taxon id, 0 = unknown
    pub lineage_ids: Vec<[u32; 6]>,
    /// interned taxon names per rank (index = id, 0 = "")
    pub taxa: [Vec<String>; 6],
    hashes: Vec<u64>,
    offsets: Vec<u64>,
    postings: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub genome: u32,
    pub shared: u32,
}

pub struct IndexBuilder {
    ksize: u32,
    scaled: u32,
    seed: u64,
    genomes: Vec<GenomeMeta>,
    sketches: Vec<Vec<u64>>,
}

impl IndexBuilder {
    pub fn new(ksize: u32, scaled: u32, seed: u64) -> Self {
        Self { ksize, scaled, seed, genomes: vec![], sketches: vec![] }
    }

    /// `hashes` must be sorted and unique (as returned by `FracSketch::finish`).
    pub fn add(&mut self, mut meta: GenomeMeta, hashes: Vec<u64>) {
        meta.n_hashes = hashes.len() as u32;
        self.genomes.push(meta);
        self.sketches.push(hashes);
    }

    /// Build the inverted index in `buckets` hash ranges to bound peak memory.
    pub fn build(self, buckets: usize) -> RefIndex {
        let max_hash = crate::max_hash_for_scaled(self.scaled);
        let buckets = buckets.max(1) as u64;
        let width = max_hash / buckets + 1;
        let mut hashes = Vec::new();
        let mut offsets = vec![0u64];
        let mut postings: Vec<u32> = Vec::new();
        for b in 0..buckets {
            let lo = b * width;
            let hi = if b + 1 == buckets { u64::MAX } else { (b + 1) * width - 1 };
            let mut pairs: Vec<(u64, u32)> = Vec::new();
            for (g, s) in self.sketches.iter().enumerate() {
                let a = s.partition_point(|&h| h < lo);
                let z = s.partition_point(|&h| h <= hi);
                pairs.extend(s[a..z].iter().map(|&h| (h, g as u32)));
            }
            pairs.sort_unstable();
            let mut prev: Option<u64> = None;
            for (h, g) in pairs {
                if prev != Some(h) {
                    if prev.is_some() {
                        offsets.push(postings.len() as u64);
                    }
                    hashes.push(h);
                    prev = Some(h);
                }
                postings.push(g);
            }
            if prev.is_some() {
                offsets.push(postings.len() as u64);
            }
        }
        RefIndex::assemble(self.ksize, self.scaled, self.seed, self.genomes, hashes, offsets, postings)
    }
}

impl RefIndex {
    fn assemble(
        ksize: u32,
        scaled: u32,
        seed: u64,
        genomes: Vec<GenomeMeta>,
        hashes: Vec<u64>,
        offsets: Vec<u64>,
        postings: Vec<u32>,
    ) -> Self {
        let mut taxa: [Vec<String>; 6] = Default::default();
        let mut maps: Vec<HashMap<String, u32>> = vec![HashMap::new(); 6];
        for r in 0..6 {
            taxa[r].push(String::new());
        }
        let mut lineage_ids = Vec::with_capacity(genomes.len());
        for g in &genomes {
            let mut ids = [0u32; 6];
            for r in 0..6 {
                let name = &g.lineage[r];
                if name.is_empty() {
                    continue;
                }
                let next = taxa[r].len() as u32;
                let id = *maps[r].entry(name.clone()).or_insert(next);
                if id == next {
                    taxa[r].push(name.clone());
                }
                ids[r] = id;
            }
            lineage_ids.push(ids);
        }
        Self { ksize, scaled, seed, genomes, lineage_ids, taxa, hashes, offsets, postings }
    }

    pub fn n_unique_hashes(&self) -> usize {
        self.hashes.len()
    }
    pub fn n_postings(&self) -> usize {
        self.postings.len()
    }

    /// Count shared hashes between `query` (sorted, unique) and every genome. Genomes for
    /// which `skip(genome)` is true are left out. Returned hits are sorted by `shared` descending
    /// (ties by genome id).
    pub fn search(&self, query: &[u64], skip: &dyn Fn(u32) -> bool) -> Vec<Hit> {
        let mut counts = vec![0u32; self.genomes.len()];
        let mut touched: Vec<u32> = Vec::new();
        let mut lo = 0usize;
        for &q in query {
            // query is sorted, so the search window only moves forward
            let pos = match self.hashes[lo..].binary_search(&q) {
                Ok(p) => lo + p,
                Err(p) => {
                    lo += p;
                    continue;
                }
            };
            lo = pos + 1;
            let (a, z) = (self.offsets[pos] as usize, self.offsets[pos + 1] as usize);
            for &g in &self.postings[a..z] {
                if skip(g) {
                    continue;
                }
                if counts[g as usize] == 0 {
                    touched.push(g);
                }
                counts[g as usize] += 1;
            }
        }
        let mut hits: Vec<Hit> =
            touched.into_iter().map(|g| Hit { genome: g, shared: counts[g as usize] }).collect();
        hits.sort_unstable_by(|a, b| b.shared.cmp(&a.shared).then(a.genome.cmp(&b.genome)));
        hits
    }

    pub fn write<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let mut meta = String::new();
        for g in &self.genomes {
            meta.push_str(&format!("{}\t{}\t{}\t{}\n", g.id, g.taxid, g.n_hashes, g.lineage.join("\t")));
        }
        w.write_all(MAGIC)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.genomes.len() as u32).to_le_bytes())?;
        w.write_all(&(meta.len() as u64).to_le_bytes())?;
        w.write_all(meta.as_bytes())?;
        w.write_all(&(self.hashes.len() as u64).to_le_bytes())?;
        wr_u64_slice(w, &self.hashes)?;
        wr_u64_slice(w, &self.offsets)?;
        w.write_all(&(self.postings.len() as u64).to_le_bytes())?;
        wr_u32_slice(w, &self.postings)
    }

    /// Write the compact v1 format: hashes as deltas, posting lists as (length, id deltas), all LEB128.
    /// `read` accepts both formats.
    pub fn write_v1<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let mut meta = String::new();
        for g in &self.genomes {
            meta.push_str(&format!("{}\t{}\t{}\t{}\n", g.id, g.taxid, g.n_hashes, g.lineage.join("\t")));
        }
        w.write_all(MAGIC_V1)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.genomes.len() as u32).to_le_bytes())?;
        w.write_all(&(meta.len() as u64).to_le_bytes())?;
        w.write_all(meta.as_bytes())?;
        w.write_all(&(self.hashes.len() as u64).to_le_bytes())?;
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
        let mut prev = 0u64;
        for (i, &h) in self.hashes.iter().enumerate() {
            put_varint(&mut buf, h - prev);
            prev = h;
            let (a, z) = (self.offsets[i] as usize, self.offsets[i + 1] as usize);
            put_varint(&mut buf, (z - a) as u64);
            let mut pg = 0u32;
            for (j, &g) in self.postings[a..z].iter().enumerate() {
                put_varint(&mut buf, if j == 0 { g as u64 } else { (g - pg) as u64 });
                pg = g;
            }
            if buf.len() >= (1 << 20) - 64 {
                w.write_all(&buf)?;
                buf.clear();
            }
        }
        w.write_all(&buf)
    }

    pub fn read<R: Read>(r: &mut R) -> io::Result<Self> {
        let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
        let mut m = [0u8; 8];
        r.read_exact(&mut m)?;
        let v1 = &m == MAGIC_V1;
        if &m != MAGIC && !v1 {
            return Err(bad("not a GmG index"));
        }
        let ksize = rd_u32(r)?;
        let scaled = rd_u32(r)?;
        let seed = rd_u64(r)?;
        let n_genomes = rd_u32(r)? as usize;
        let meta_len = rd_u64(r)? as usize;
        let mut mb = vec![0u8; meta_len];
        r.read_exact(&mut mb)?;
        let meta = String::from_utf8(mb).map_err(|_| bad("meta not utf-8"))?;
        let mut genomes = Vec::with_capacity(n_genomes);
        for line in meta.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 9 {
                return Err(bad("bad genome record"));
            }
            let mut lineage: [String; 6] = Default::default();
            for r in 0..6 {
                lineage[r] = f[3 + r].to_string();
            }
            genomes.push(GenomeMeta {
                id: f[0].to_string(),
                taxid: f[1].parse().map_err(|_| bad("taxid"))?,
                n_hashes: f[2].parse().map_err(|_| bad("n_hashes"))?,
                lineage,
            });
        }
        if genomes.len() != n_genomes {
            return Err(bad("genome count mismatch"));
        }
        let nh = rd_u64(r)? as usize;
        if v1 {
            let mut rest = Vec::new();
            r.read_to_end(&mut rest)?;
            let mut pos = 0usize;
            let mut hashes = Vec::with_capacity(nh);
            let mut offsets = Vec::with_capacity(nh + 1);
            let mut postings: Vec<u32> = Vec::new();
            offsets.push(0u64);
            let mut prev = 0u64;
            for _ in 0..nh {
                prev += get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))?;
                hashes.push(prev);
                let len = get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))? as usize;
                let mut g = 0u32;
                for j in 0..len {
                    let v = get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))? as u32;
                    g = if j == 0 { v } else { g + v };
                    postings.push(g);
                }
                offsets.push(postings.len() as u64);
            }
            return Ok(Self::assemble(ksize, scaled, seed, genomes, hashes, offsets, postings));
        }
        let hashes = rd_u64_vec(r, nh)?;
        let offsets = rd_u64_vec(r, nh + 1)?;
        let np = rd_u64(r)? as usize;
        let postings = rd_u32_vec(r, np)?;
        Ok(Self::assemble(ksize, scaled, seed, genomes, hashes, offsets, postings))
    }
}

fn put_varint(buf: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        buf.push((v as u8) | 0x80);
        v >>= 7;
    }
    buf.push(v as u8);
}

fn get_varint(buf: &[u8], pos: &mut usize) -> Option<u64> {
    let (mut v, mut shift) = (0u64, 0u32);
    loop {
        let b = *buf.get(*pos)?;
        *pos += 1;
        v |= ((b & 0x7f) as u64) << shift;
        if b < 0x80 {
            return Some(v);
        }
        shift += 7;
        if shift > 63 {
            return None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn v1_roundtrip_matches_v0() {
        let ix = small(3);
        let mut a = Vec::new();
        ix.write(&mut a).unwrap();
        let mut b = Vec::new();
        ix.write_v1(&mut b).unwrap();
        let back = RefIndex::read(&mut &b[..]).unwrap();
        assert_eq!(back.genomes, ix.genomes);
        assert_eq!(back.n_unique_hashes(), ix.n_unique_hashes());
        assert_eq!(back.n_postings(), ix.n_postings());
        assert_eq!(back.search(&[10, 20, 30], &|_| false), ix.search(&[10, 20, 30], &|_| false));
        let mut c = Vec::new();
        back.write(&mut c).unwrap();
        assert_eq!(a, c);
    }

    #[test]
    fn varint_edges() {
        for v in [0u64, 1, 127, 128, 300, u32::MAX as u64, u64::MAX] {
            let mut b = Vec::new();
            put_varint(&mut b, v);
            let mut p = 0;
            assert_eq!(get_varint(&b, &mut p), Some(v));
            assert_eq!(p, b.len());
        }
    }

    fn meta(id: &str, species: &str) -> GenomeMeta {
        let mut lineage: [String; 6] = Default::default();
        lineage[4] = "G".into();
        lineage[5] = species.into();
        GenomeMeta { id: id.into(), taxid: 1, lineage, n_hashes: 0 }
    }

    fn small(buckets: usize) -> RefIndex {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(meta("a", "G a"), vec![10, 20, 30, 40]);
        b.add(meta("b", "G b"), vec![20, 30, 50]);
        b.add(meta("c", "G a"), vec![10, 60]);
        b.build(buckets)
    }

    #[test]
    fn search_counts_shared() {
        for buckets in [1, 3, 8] {
            let ix = small(buckets);
            let hits = ix.search(&[10, 20, 30, 99], &|_| false);
            assert_eq!(hits[0], Hit { genome: 0, shared: 3 });
            assert_eq!(hits[1], Hit { genome: 1, shared: 2 });
            assert_eq!(hits[2], Hit { genome: 2, shared: 1 });
            let hits = ix.search(&[10, 20, 30], &|g| g == 0);
            assert_eq!(hits[0], Hit { genome: 1, shared: 2 });
        }
    }

    #[test]
    fn roundtrip_and_interning() {
        let ix = small(2);
        let mut buf = Vec::new();
        ix.write(&mut buf).unwrap();
        let back = RefIndex::read(&mut &buf[..]).unwrap();
        assert_eq!(back.genomes, ix.genomes);
        assert_eq!(back.search(&[10, 20], &|_| false), ix.search(&[10, 20], &|_| false));
        assert_eq!(back.lineage_ids[0][5], back.lineage_ids[2][5]);
        assert_ne!(back.lineage_ids[0][5], back.lineage_ids[1][5]);
        assert_eq!(back.lineage_ids[0][0], 0);
    }
}
