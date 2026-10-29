//! Reference index: genomes with a 6-rank lineage and an inverted hash -> genomes map.
//!
//! In memory and on disk the map is three parallel arrays: sorted unique `hashes`,
//! `offsets` (len = hashes + 1) and `postings` (genome ids, ascending within a hash).
//! A search of one query sketch is one binary search per query hash.
//!
//! Format v3 (`GMGIDX03`) keeps the map packed (see `packed.rs`): loading copies the file and reads a small block
//! directory, and a query decodes only the blocks it touches, so memory stays near the file size.

use crate::packed::Packed;
use crate::skfile::{rd_u32, rd_u32_vec, rd_u64, rd_u64_vec, wr_u32_slice, wr_u64_slice};
use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{self, Read, Write};

pub const RANKS: [&str; 6] = ["phylum", "class", "order", "family", "genus", "species"];
const MAGIC: &[u8; 8] = b"GMGIDX01";
const MAGIC_V1: &[u8; 8] = b"GMGIDX02"; // delta + LEB128 varint coded postings
const MAGIC_V3: &[u8; 8] = b"GMGIDX03"; // the same coding, packed in blocks with a directory (loaded as is)

#[derive(Clone, Debug, PartialEq, Default)]
pub struct GenomeMeta {
    pub id: String,
    pub taxid: u32,
    /// phylum, class, order, family, genus, species; empty string = unknown
    pub lineage: [String; 6],
    /// number of retained hashes of this genome (filled by the builder)
    pub n_hashes: u32,
    /// ANI-cluster id (the centroid genome's id from `gmg cluster`); empty if not clustered
    pub cluster: String,
    /// NCBI taxid of the taxon at each rank (0 = unknown); `lineage` holds the names
    pub rank_taxids: [u32; 6],
    /// the original label of the genome in the source collection (e.g. BFD species string), kept next to the NCBI lineage
    pub label: String,
}

/// hashes, offsets and postings, borrowed or expanded
type Arrays<'a> = (Cow<'a, [u64]>, Cow<'a, [u64]>, Cow<'a, [u32]>);

pub struct RefIndex {
    pub ksize: u32,
    pub scaled: u32,
    pub seed: u64,
    pub genomes: Vec<GenomeMeta>,
    /// per genome, per rank: interned taxon id, 0 = unknown
    pub lineage_ids: Vec<[u32; 6]>,
    /// interned taxon names per rank (index = id, 0 = "")
    pub taxa: [Vec<String>; 6],
    /// NCBI taxid of each interned taxon (parallel to `taxa`; 0 = unknown)
    pub taxon_taxids: [Vec<u32>; 6],
    /// optional provenance record (one line of JSON) stored in the index: how and from what it was built
    pub manifest: Option<String>,
    hashes: Vec<u64>,
    offsets: Vec<u64>,
    postings: Vec<u32>,
    /// set instead of the three arrays when loaded from a v3 file
    packed: Option<Packed>,
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
        Self {
            ksize,
            scaled,
            seed,
            genomes: vec![],
            sketches: vec![],
        }
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
        let width = (max_hash / buckets).saturating_add(1);
        let mut hashes = Vec::new();
        let mut offsets = vec![0u64];
        let mut postings: Vec<u32> = Vec::new();
        for b in 0..buckets {
            let lo = b * width;
            let hi = if b + 1 == buckets {
                u64::MAX
            } else {
                (b + 1) * width - 1
            };
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
        RefIndex::assemble(
            self.ksize,
            self.scaled,
            self.seed,
            self.genomes,
            hashes,
            offsets,
            postings,
        )
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
        let mut taxon_taxids: [Vec<u32>; 6] = Default::default();
        for t in taxa.iter_mut() {
            t.push(String::new());
        }
        for t in taxon_taxids.iter_mut() {
            t.push(0);
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
                    taxon_taxids[r].push(0);
                }
                // homonyms share a taxon here; the first non-zero taxid wins
                if taxon_taxids[r][id as usize] == 0 {
                    taxon_taxids[r][id as usize] = g.rank_taxids[r];
                }
                ids[r] = id;
            }
            lineage_ids.push(ids);
        }
        Self {
            ksize,
            scaled,
            seed,
            genomes,
            lineage_ids,
            taxa,
            taxon_taxids,
            manifest: None,
            hashes,
            offsets,
            postings,
            packed: None,
        }
    }

    fn manifest_line(&self) -> String {
        match &self.manifest {
            Some(m) => format!("#manifest\t{}\n", m.replace(['\n', '\r'], " ")),
            None => String::new(),
        }
    }

    /// Attach a provenance record (JSON on one line) that is written with the index.
    pub fn set_manifest(&mut self, json: String) {
        self.manifest = Some(json);
    }

    /// Expected genome size (Mb, estimated from the sketch: hashes x scaled) of reference genomes, as the 1st-99th percentile
    /// of the called phylum when it has at least `MIN_FOR_RANGE` genomes, else the 0.5th-99.5th percentile of all genomes.
    /// None when the reference is too small to say. Returns (low, high, basis).
    pub fn size_range_mb(&self, phylum: Option<&str>) -> Option<(f64, f64, String)> {
        const MIN_FOR_RANGE: usize = 50;
        let mb = |g: &GenomeMeta| g.n_hashes as f64 * self.scaled as f64 / 1e6;
        let pct = |v: &mut Vec<f64>, p: f64| {
            v.sort_by(|a, b| a.total_cmp(b));
            v[((p * v.len() as f64) as usize).min(v.len() - 1)]
        };
        if let Some(ph) = phylum.filter(|p| !p.is_empty()) {
            let mut v: Vec<f64> = self
                .genomes
                .iter()
                .filter(|g| g.lineage[0] == ph)
                .map(mb)
                .collect();
            if v.len() >= MIN_FOR_RANGE {
                let n = v.len();
                return Some((
                    pct(&mut v, 0.01),
                    pct(&mut v, 0.99),
                    format!("{ph} (n = {n})"),
                ));
            }
        }
        let mut v: Vec<f64> = self.genomes.iter().map(mb).collect();
        if v.len() < MIN_FOR_RANGE {
            return None;
        }
        let n = v.len();
        Some((
            pct(&mut v, 0.005),
            pct(&mut v, 0.995),
            format!("all reference genomes (n = {n})"),
        ))
    }

    pub fn n_unique_hashes(&self) -> usize {
        match &self.packed {
            Some(p) => p.n_hashes as usize,
            None => self.hashes.len(),
        }
    }
    pub fn n_postings(&self) -> usize {
        match &self.packed {
            Some(p) => p.n_postings as usize,
            None => self.postings.len(),
        }
    }

    /// True when the postings are kept packed (loaded from a v3 file).
    pub fn is_packed(&self) -> bool {
        self.packed.is_some()
    }

    /// Expand a packed index into the in-memory arrays (about 3x the memory, faster for very many queries: bench, cluster).
    /// An index that is already expanded is returned unchanged; a corrupt packed stream is an error.
    pub fn into_memory(mut self) -> io::Result<RefIndex> {
        if let Some(p) = self.packed.take() {
            let (h, o, g) = p.expand().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "corrupt index (postings)")
            })?;
            self.hashes = h;
            self.offsets = o;
            self.postings = g;
        }
        Ok(self)
    }

    /// Call `f` with the genome ids of every hash of `query` (sorted, unique) that is in the index.
    fn visit(&self, query: &[u64], f: &mut dyn FnMut(&[u32])) {
        if let Some(p) = &self.packed {
            return p.visit(query, f);
        }
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
            f(&self.postings[a..z]);
        }
    }

    /// For each phylum, the number of `query` hashes found only in genomes of that phylum (sorted, unique query).
    /// Hashes shared across phyla or only with genomes of unknown phylum are not counted. Returns (phylum name, count)
    /// sorted by count descending; a single-organism input has one dominant phylum, a mixture has two.
    pub fn phylum_exclusive(&self, query: &[u64]) -> Vec<(String, u32)> {
        let mut counts = vec![0u32; self.taxa[0].len().max(1)];
        let n = self.lineage_ids.len();
        self.visit(query, &mut |gs| {
            let mut ph = None;
            let mut single = true;
            for &g in gs {
                let id = if (g as usize) < n {
                    self.lineage_ids[g as usize][0]
                } else {
                    0
                };
                if id == 0 || ph.is_some_and(|p| p != id) {
                    single = false;
                    break;
                }
                ph = Some(id);
            }
            if let (true, Some(id)) = (single, ph) {
                counts[id as usize] += 1;
            }
        });
        let mut v: Vec<(String, u32)> = counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(i, &c)| (self.taxa[0][i].clone(), c))
            .collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        v
    }

    /// Count shared hashes between `query` (sorted, unique) and every genome. Genomes for
    /// which `skip(genome)` is true are left out. Returned hits are sorted by `shared` descending
    /// (ties by genome id).
    pub fn search(&self, query: &[u64], skip: &dyn Fn(u32) -> bool) -> Vec<Hit> {
        let mut counts = vec![0u32; self.genomes.len()];
        let mut touched: Vec<u32> = Vec::new();
        self.visit(query, &mut |gs| {
            for &g in gs {
                if (g as usize) >= counts.len() || skip(g) {
                    continue;
                }
                if counts[g as usize] == 0 {
                    touched.push(g);
                }
                counts[g as usize] += 1;
            }
        });
        let mut hits: Vec<Hit> = touched
            .into_iter()
            .map(|g| Hit {
                genome: g,
                shared: counts[g as usize],
            })
            .collect();
        hits.sort_unstable_by(|a, b| b.shared.cmp(&a.shared).then(a.genome.cmp(&b.genome)));
        hits
    }

    /// The metadata section of the file: manifest line, then one tab-separated line per genome.
    fn meta_text(&self) -> String {
        let mut meta = self.manifest_line();
        for g in &self.genomes {
            meta.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                g.id,
                g.taxid,
                g.n_hashes,
                g.lineage.join("\t"),
                g.cluster,
                g.rank_taxids.map(|t| t.to_string()).join("\t"),
                g.label.replace(['\t', '\n', '\r'], " ")
            ));
        }
        meta
    }

    /// The three arrays, expanded on the fly when the index is packed.
    fn arrays(&self) -> io::Result<Arrays<'_>> {
        match &self.packed {
            None => Ok((
                Cow::Borrowed(&self.hashes),
                Cow::Borrowed(&self.offsets),
                Cow::Borrowed(&self.postings),
            )),
            Some(p) => {
                let (h, o, g) = p.expand().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "corrupt index (postings)")
                })?;
                Ok((Cow::Owned(h), Cow::Owned(o), Cow::Owned(g)))
            }
        }
    }

    /// Write format v3 (packed). Loads fast and stays small in memory; `gmg convert` writes it.
    pub fn write_v3<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let meta = self.meta_text();
        let encoded;
        let packed = match &self.packed {
            Some(p) => p,
            None => {
                encoded = Packed::encode(&self.hashes, &self.offsets, &self.postings);
                &encoded
            }
        };
        w.write_all(MAGIC_V3)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.genomes.len() as u32).to_le_bytes())?;
        w.write_all(&(meta.len() as u64).to_le_bytes())?;
        w.write_all(meta.as_bytes())?;
        w.write_all(&packed.n_hashes.to_le_bytes())?;
        w.write_all(&packed.n_postings.to_le_bytes())?;
        w.write_all(&(packed.first.len() as u64).to_le_bytes())?;
        wr_u64_slice(w, &packed.first)?;
        wr_u64_slice(w, &packed.off)?;
        w.write_all(&(packed.data.len() as u64).to_le_bytes())?;
        w.write_all(&packed.data)
    }

    pub fn write<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let (hashes, offsets, postings) = self.arrays()?;
        let meta = self.meta_text();
        w.write_all(MAGIC)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.genomes.len() as u32).to_le_bytes())?;
        w.write_all(&(meta.len() as u64).to_le_bytes())?;
        w.write_all(meta.as_bytes())?;
        w.write_all(&(hashes.len() as u64).to_le_bytes())?;
        wr_u64_slice(w, &hashes)?;
        wr_u64_slice(w, &offsets)?;
        w.write_all(&(postings.len() as u64).to_le_bytes())?;
        wr_u32_slice(w, &postings)
    }

    /// Write the compact v1 format: hashes as deltas, posting lists as (length, id deltas), all LEB128.
    /// `read` accepts both formats.
    pub fn write_v1<W: Write>(&self, w: &mut W) -> io::Result<()> {
        let meta = self.meta_text();
        let (hashes, offsets, postings) = self.arrays()?;
        w.write_all(MAGIC_V1)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.genomes.len() as u32).to_le_bytes())?;
        w.write_all(&(meta.len() as u64).to_le_bytes())?;
        w.write_all(meta.as_bytes())?;
        w.write_all(&(hashes.len() as u64).to_le_bytes())?;
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
        let mut prev = 0u64;
        for (i, &h) in hashes.iter().enumerate() {
            put_varint(&mut buf, h - prev);
            prev = h;
            let (a, z) = (offsets[i] as usize, offsets[i + 1] as usize);
            put_varint(&mut buf, (z - a) as u64);
            let mut pg = 0u32;
            for (j, &g) in postings[a..z].iter().enumerate() {
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
        Ok(Self::read_impl(r, false)?.0)
    }

    /// Load from an owned buffer. A v3 file keeps the buffer as its packed data (no second copy: the header is
    /// dropped in place); older formats are decoded as in `read`.
    pub fn read_owned(mut bytes: Vec<u8>) -> io::Result<Self> {
        if bytes.len() < 8 || &bytes[..8] != MAGIC_V3 {
            return Self::read(&mut &bytes[..]);
        }
        let mut cur = io::Cursor::new(&bytes[..]);
        let (mut ix, dl) = Self::read_impl(&mut cur, true)?;
        let pos = cur.position() as usize;
        if bytes.len() < pos.saturating_add(dl) {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "file is truncated",
            ));
        }
        bytes.truncate(pos + dl);
        bytes.drain(..pos);
        if let Some(p) = ix.packed.as_mut() {
            p.data = bytes;
            p.check()?;
        }
        Ok(ix)
    }

    /// `skip_data`: leave the v3 data section unread and return its length (for `read_owned`).
    fn read_impl<R: Read>(r: &mut R, skip_data: bool) -> io::Result<(Self, usize)> {
        let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
        let mut m = [0u8; 8];
        r.read_exact(&mut m)?;
        let v1 = &m == MAGIC_V1;
        let v3 = &m == MAGIC_V3;
        if &m != MAGIC && !v1 && !v3 {
            return Err(bad("not a GmG index"));
        }
        let ksize = rd_u32(r)?;
        let scaled = rd_u32(r)?;
        let seed = rd_u64(r)?;
        crate::skfile::check_params(ksize, scaled)?;
        let n_genomes = rd_u32(r)? as usize;
        let meta_len = usize::try_from(rd_u64(r)?).map_err(|_| bad("implausible length"))?;
        let mb = crate::skfile::rd_bytes(r, meta_len)?;
        let meta = String::from_utf8(mb).map_err(|_| bad("meta not utf-8"))?;
        let mut genomes = Vec::with_capacity(n_genomes.min(1 << 20));
        let mut manifest = None;
        for line in meta.lines() {
            if let Some(m) = line.strip_prefix("#manifest\t") {
                manifest = Some(m.to_string());
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            // 9 fields = indexes written before clusters existed; 10 = with the cluster column;
            // 17 = also the NCBI taxid of every rank and the original label
            if f.len() != 9 && f.len() != 10 && f.len() != 17 {
                return Err(bad("bad genome record"));
            }
            let mut lineage: [String; 6] = Default::default();
            for r in 0..6 {
                lineage[r] = f[3 + r].to_string();
            }
            let mut rank_taxids = [0u32; 6];
            if f.len() == 17 {
                for r in 0..6 {
                    rank_taxids[r] = f[10 + r].parse().map_err(|_| bad("rank taxid"))?;
                }
            }
            genomes.push(GenomeMeta {
                id: f[0].to_string(),
                taxid: f[1].parse().map_err(|_| bad("taxid"))?,
                n_hashes: f[2].parse().map_err(|_| bad("n_hashes"))?,
                lineage,
                cluster: f.get(9).map(|c| c.to_string()).unwrap_or_default(),
                rank_taxids,
                label: f.get(16).map(|c| c.to_string()).unwrap_or_default(),
            });
        }
        if genomes.len() != n_genomes {
            return Err(bad("genome count mismatch"));
        }
        let nh = usize::try_from(rd_u64(r)?).map_err(|_| bad("implausible length"))?;
        if v3 {
            let np = rd_u64(r)?;
            let nb = usize::try_from(rd_u64(r)?).map_err(|_| bad("implausible length"))?;
            let first = rd_u64_vec(r, nb)?;
            let off = rd_u64_vec(r, nb)?;
            let dl = usize::try_from(rd_u64(r)?).map_err(|_| bad("implausible length"))?;
            let data = if skip_data {
                Vec::new()
            } else {
                crate::skfile::rd_bytes(r, dl)?
            };
            let packed = Packed {
                n_hashes: nh as u64,
                n_postings: np,
                first,
                off,
                data,
            };
            if !skip_data {
                packed.check()?;
            }
            let mut ix = Self::assemble(
                ksize,
                scaled,
                seed,
                genomes,
                Vec::new(),
                vec![0],
                Vec::new(),
            );
            ix.packed = Some(packed);
            ix.manifest = manifest;
            return Ok((ix, dl));
        }
        if v1 {
            let mut rest = Vec::new();
            r.read_to_end(&mut rest)?;
            let mut pos = 0usize;
            if nh > rest.len() {
                return Err(bad("truncated index"));
            }
            let mut hashes = Vec::with_capacity(nh);
            let mut offsets = Vec::with_capacity(nh + 1);
            let mut postings: Vec<u32> = Vec::new();
            offsets.push(0u64);
            let mut prev = 0u64;
            for _ in 0..nh {
                prev = prev
                    .checked_add(get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))?)
                    .ok_or_else(|| bad("corrupt index"))?;
                hashes.push(prev);
                let len =
                    get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))? as usize;
                let mut g = 0u32;
                for j in 0..len {
                    let v = u32::try_from(
                        get_varint(&rest, &mut pos).ok_or_else(|| bad("truncated index"))?,
                    )
                    .map_err(|_| bad("corrupt index"))?;
                    g = if j == 0 {
                        v
                    } else {
                        g.checked_add(v).ok_or_else(|| bad("corrupt index"))?
                    };
                    postings.push(g);
                }
                offsets.push(postings.len() as u64);
            }
            validate(genomes.len(), &hashes, &offsets, &postings)?;
            let mut ix = Self::assemble(ksize, scaled, seed, genomes, hashes, offsets, postings);
            ix.manifest = manifest;
            return Ok((ix, 0));
        }
        let hashes = rd_u64_vec(r, nh)?;
        let offsets = rd_u64_vec(
            r,
            nh.checked_add(1).ok_or_else(|| bad("implausible length"))?,
        )?;
        let np = usize::try_from(rd_u64(r)?).map_err(|_| bad("implausible length"))?;
        let postings = rd_u32_vec(r, np)?;
        validate(genomes.len(), &hashes, &offsets, &postings)?;
        let mut ix = Self::assemble(ksize, scaled, seed, genomes, hashes, offsets, postings);
        ix.manifest = manifest;
        Ok((ix, 0))
    }
}

/// Structural checks so a corrupt or hostile index fails to load instead of panicking during a search.
fn validate(n_genomes: usize, hashes: &[u64], offsets: &[u64], postings: &[u32]) -> io::Result<()> {
    let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
    if offsets.len() != hashes.len() + 1 || offsets.first() != Some(&0) {
        return Err(bad("corrupt index (offsets)"));
    }
    if offsets.last().copied() != Some(postings.len() as u64)
        || offsets.windows(2).any(|w| w[0] > w[1])
    {
        return Err(bad("corrupt index (offsets)"));
    }
    if hashes.windows(2).any(|w| w[0] >= w[1]) {
        return Err(bad("corrupt index (hashes not sorted)"));
    }
    if postings.iter().any(|&g| g as usize >= n_genomes) {
        return Err(bad("corrupt index (genome id out of range)"));
    }
    Ok(())
}

pub(crate) fn put_varint(buf: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        buf.push((v as u8) | 0x80);
        v >>= 7;
    }
    buf.push(v as u8);
}

pub(crate) fn get_varint(buf: &[u8], pos: &mut usize) -> Option<u64> {
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
        assert_eq!(
            back.search(&[10, 20, 30], &|_| false),
            ix.search(&[10, 20, 30], &|_| false)
        );
        let mut c = Vec::new();
        back.write(&mut c).unwrap();
        assert_eq!(a, c);
    }

    #[test]
    fn phylum_exclusive_separates_a_mixture() {
        let meta = |id: &str, ph: &str| {
            let mut l: [String; 6] = Default::default();
            l[0] = ph.into();
            GenomeMeta {
                id: id.into(),
                taxid: 1,
                lineage: l,
                n_hashes: 0,
                cluster: String::new(),
                ..Default::default()
            }
        };
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(meta("a", "Ascomycota"), (0..100).collect());
        b.add(meta("b", "Basidiomycota"), (1000..1100).collect());
        b.add(meta("c", "Basidiomycota"), (50..150).collect()); // 50..100 shared across phyla
        let ix = b.build(4);
        let mut q: Vec<u64> = (0..30).collect();
        q.extend(1000..1020);
        q.extend(60..70);
        q.sort_unstable();
        let v = ix.phylum_exclusive(&q);
        assert_eq!(
            v,
            vec![("Ascomycota".into(), 30), ("Basidiomycota".into(), 20)]
        );
    }

    #[test]
    fn size_range_uses_the_phylum_when_it_is_large_enough_else_everything() {
        let meta = |i: usize, ph: &str| {
            let mut l: [String; 6] = Default::default();
            l[0] = ph.into();
            GenomeMeta {
                id: format!("g{i}"),
                taxid: 1,
                lineage: l,
                n_hashes: 0,
                cluster: String::new(),
                ..Default::default()
            }
        };
        let mut b = IndexBuilder::new(21, 1000, 42);
        // 60 Ascomycota of 1000..1059 hashes (1.0-1.06 Mb at scaled 1000), 5 Basidiomycota of 9000 hashes
        for i in 0..60 {
            b.add(meta(i, "Ascomycota"), (0..1000 + i as u64).collect());
        }
        for i in 60..65 {
            b.add(meta(i, "Basidiomycota"), (0..9000).collect());
        }
        let ix = b.build(4);
        let (lo, hi, basis) = ix.size_range_mb(Some("Ascomycota")).unwrap();
        assert!(
            (1.0..=1.01).contains(&lo) && (1.05..=1.06).contains(&hi),
            "{lo} {hi}"
        );
        assert!(basis.starts_with("Ascomycota (n = 60"));
        // too few Basidiomycota: falls back to all genomes
        let (_, hi_all, basis) = ix.size_range_mb(Some("Basidiomycota")).unwrap();
        assert!(
            basis.starts_with("all reference genomes") && hi_all >= 9.0,
            "{hi_all}"
        );
        assert!(ix.size_range_mb(None).is_some());
    }

    #[test]
    fn manifest_survives_both_formats_and_does_not_count_as_a_genome() {
        let mut ix = small(3);
        assert!(ix.manifest.is_none());
        ix.set_manifest("{\"k\":21,\"note\":\"two\\nlines\"}".to_string());
        for v1 in [false, true] {
            let mut b = Vec::new();
            if v1 {
                ix.write_v1(&mut b).unwrap();
            } else {
                ix.write(&mut b).unwrap();
            }
            let back = RefIndex::read(&mut &b[..]).unwrap();
            assert_eq!(back.genomes, ix.genomes);
            assert!(back.manifest.as_deref().unwrap().starts_with("{\"k\":21"));
        }
    }

    #[test]
    fn rank_taxids_and_labels_survive_both_formats_and_follow_the_taxon() {
        let mut b = IndexBuilder::new(21, 1, 42);
        let mk = |id: &str, genus: &str, gt: u32, species: &str, st: u32, label: &str| {
            let mut l: [String; 6] = Default::default();
            l[0] = "Ascomycota".into();
            l[4] = genus.into();
            l[5] = species.into();
            GenomeMeta {
                id: id.into(),
                taxid: st,
                lineage: l,
                rank_taxids: [4890, 0, 0, 0, gt, st],
                label: label.into(),
                ..Default::default()
            }
        };
        b.add(
            mk(
                "a",
                "Fusarium",
                5506,
                "Fusarium oxysporum",
                5507,
                "F. oxysporum f. sp. x",
            ),
            (0..50).collect(),
        );
        b.add(
            mk(
                "b",
                "Fusarium",
                5506,
                "Fusarium solani",
                169388,
                "Fusarium solani",
            ),
            (40..90).collect(),
        );
        let ix = b.build(2);
        let g = ix.taxa[4].iter().position(|t| t == "Fusarium").unwrap();
        assert_eq!(ix.taxon_taxids[4][g], 5506);
        for v1 in [false, true] {
            let mut buf = Vec::new();
            if v1 {
                ix.write_v1(&mut buf).unwrap();
            } else {
                ix.write(&mut buf).unwrap();
            }
            let back = RefIndex::read(&mut &buf[..]).unwrap();
            assert_eq!(back.genomes, ix.genomes);
            assert_eq!(back.genomes[0].rank_taxids[5], 5507);
            assert_eq!(back.genomes[0].label, "F. oxysporum f. sp. x");
            let g = back.taxa[4].iter().position(|t| t == "Fusarium").unwrap();
            assert_eq!(back.taxon_taxids[4][g], 5506);
        }
    }

    #[test]
    fn packed_v3_answers_exactly_like_the_expanded_index() {
        // 40 genomes x ~900 hashes drawn from a shared pool so lists have 1..many genomes; > 3 blocks of 256 hashes
        let mut rng = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let mut b = IndexBuilder::new(21, 1, 42);
        let pool: Vec<u64> = (0..2500).map(|_| next() >> 3).collect();
        for i in 0..40 {
            let mut h: Vec<u64> = (0..900).map(|_| pool[(next() % 2500) as usize]).collect();
            h.sort_unstable();
            h.dedup();
            b.add(meta(&format!("g{i}"), "s"), h);
        }
        let ix = b.build(4);
        assert!(ix.n_unique_hashes() > 3 * 256);
        let mut buf = Vec::new();
        ix.write_v3(&mut buf).unwrap();
        let packed = RefIndex::read(&mut &buf[..]).unwrap();
        assert!(packed.is_packed() && !ix.is_packed());
        assert_eq!(packed.n_unique_hashes(), ix.n_unique_hashes());
        assert_eq!(packed.n_postings(), ix.n_postings());
        // queries: present hashes, absent hashes, below the smallest and above the largest, several sizes
        for size in [1usize, 5, 300, 2000] {
            let mut q: Vec<u64> = (0..size)
                .map(|i| match i % 4 {
                    0 => pool[(next() % 2500) as usize],
                    1 => next() >> 3,
                    2 => i as u64,
                    _ => u64::MAX - i as u64,
                })
                .collect();
            q.sort_unstable();
            q.dedup();
            assert_eq!(
                packed.search(&q, &|_| false),
                ix.search(&q, &|_| false),
                "size {size}"
            );
            assert_eq!(
                packed.search(&q, &|g| g % 3 == 0),
                ix.search(&q, &|g| g % 3 == 0)
            );
            assert_eq!(packed.phylum_exclusive(&q), ix.phylum_exclusive(&q));
        }
        // expanding gives back the arrays: writing v1 from either is identical, and v3 survives a rewrite
        let (mut a, mut c) = (Vec::new(), Vec::new());
        ix.write_v1(&mut a).unwrap();
        packed.write_v1(&mut c).unwrap();
        assert_eq!(a, c);
        let mut again = Vec::new();
        packed.write_v3(&mut again).unwrap();
        assert_eq!(again, buf);
        // the owned loader gives the same index without a second copy
        let owned = RefIndex::read_owned(buf.clone()).unwrap();
        assert!(owned.is_packed());
        let q = [pool[1], pool[2], pool[3]];
        assert_eq!(owned.search(&q, &|_| false), ix.search(&q, &|_| false));
        let mem = RefIndex::read(&mut &buf[..])
            .unwrap()
            .into_memory()
            .unwrap();
        assert!(!mem.is_packed());
        assert_eq!(
            mem.search(&[pool[0]], &|_| false),
            ix.search(&[pool[0]], &|_| false)
        );
        // and the packed file is no bigger than v1 (plus the small directory)
        assert!(
            buf.len() < a.len() + a.len() / 20 + 4096,
            "{} vs {}",
            buf.len(),
            a.len()
        );
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
        GenomeMeta {
            id: id.into(),
            taxid: 1,
            lineage,
            n_hashes: 0,
            cluster: String::new(),
            ..Default::default()
        }
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
            assert_eq!(
                hits[0],
                Hit {
                    genome: 0,
                    shared: 3
                }
            );
            assert_eq!(
                hits[1],
                Hit {
                    genome: 1,
                    shared: 2
                }
            );
            assert_eq!(
                hits[2],
                Hit {
                    genome: 2,
                    shared: 1
                }
            );
            let hits = ix.search(&[10, 20, 30], &|g| g == 0);
            assert_eq!(
                hits[0],
                Hit {
                    genome: 1,
                    shared: 2
                }
            );
        }
    }

    #[test]
    fn roundtrip_and_interning() {
        let ix = small(2);
        let mut buf = Vec::new();
        ix.write(&mut buf).unwrap();
        let back = RefIndex::read(&mut &buf[..]).unwrap();
        assert_eq!(back.genomes, ix.genomes);
        assert_eq!(
            back.search(&[10, 20], &|_| false),
            ix.search(&[10, 20], &|_| false)
        );
        assert_eq!(back.lineage_ids[0][5], back.lineage_ids[2][5]);
        assert_ne!(back.lineage_ids[0][5], back.lineage_ids[1][5]);
        assert_eq!(back.lineage_ids[0][0], 0);
    }
}
