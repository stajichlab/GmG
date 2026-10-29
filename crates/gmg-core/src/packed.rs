//! Packed postings (index format v3): the inverted index stays delta/varint coded in memory (or in a mapped file) and a
//! block directory lets a query decode only the blocks that hold its hashes. Loading is a copy of the file plus a
//! small directory; resident memory is about the file size instead of three expanded arrays (3x larger).
//!
//! Layout of the data stream: entries in hash order, `BLOCK` entries per block. Entry = [hash delta from the previous
//! entry of the same block, omitted for a block's first entry (its hash is in the directory)] [number of genomes]
//! [genome ids, delta coded]. The directory holds, per block, the first hash and the byte offset of its first entry.

use crate::index::{get_varint, put_varint};
use std::io;

/// hashes per block: small enough that a lookup decodes ~256 entries, large enough that the directory is tiny
pub(crate) const BLOCK: usize = 256;

#[derive(Clone, Debug)]
pub(crate) struct Packed {
    pub(crate) n_hashes: u64,
    pub(crate) n_postings: u64,
    pub(crate) first: Vec<u64>,
    pub(crate) off: Vec<u64>,
    pub(crate) data: Vec<u8>,
}

fn bad(m: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m.to_string())
}

impl Packed {
    /// Build from the expanded arrays of a freshly built index.
    pub(crate) fn encode(hashes: &[u64], offsets: &[u64], postings: &[u32]) -> Packed {
        let mut data = Vec::with_capacity(postings.len() * 2 + hashes.len() * 2);
        let (mut first, mut off) = (Vec::new(), Vec::new());
        let mut prev = 0u64;
        for (i, &h) in hashes.iter().enumerate() {
            if i % BLOCK == 0 {
                first.push(h);
                off.push(data.len() as u64);
            } else {
                put_varint(&mut data, h - prev);
            }
            prev = h;
            let (a, z) = (offsets[i] as usize, offsets[i + 1] as usize);
            put_varint(&mut data, (z - a) as u64);
            let mut g = 0u32;
            for (j, &v) in postings[a..z].iter().enumerate() {
                put_varint(&mut data, if j == 0 { v } else { v - g } as u64);
                g = v;
            }
        }
        Packed {
            n_hashes: hashes.len() as u64,
            n_postings: postings.len() as u64,
            first,
            off,
            data,
        }
    }

    /// Cheap structural checks done at load time (the stream itself is only decoded on demand, defensively).
    pub(crate) fn check(&self) -> io::Result<()> {
        // every entry takes at least one byte and every posting at least one: counts beyond the data are impossible
        if self.n_hashes > self.data.len() as u64 || self.n_postings > self.data.len() as u64 {
            return Err(bad("corrupt index (counts exceed the data)"));
        }
        let nb = self.n_hashes.div_ceil(BLOCK as u64) as usize;
        if self.first.len() != nb || self.off.len() != nb {
            return Err(bad("corrupt index (block directory)"));
        }
        if self.first.windows(2).any(|w| w[0] >= w[1])
            || self.off.windows(2).any(|w| w[0] > w[1])
            || self
                .off
                .last()
                .is_some_and(|&o| o as usize > self.data.len())
        {
            return Err(bad("corrupt index (block directory order)"));
        }
        Ok(())
    }

    fn entries_in(&self, b: usize) -> usize {
        if b + 1 < self.first.len() {
            BLOCK
        } else {
            (self.n_hashes as usize)
                .saturating_sub(b * BLOCK)
                .min(BLOCK)
        }
    }

    /// Call `f` with the genome ids of every query hash present (query sorted, unique). A corrupt stream ends the walk
    /// early instead of panicking.
    pub(crate) fn visit(&self, query: &[u64], f: &mut dyn FnMut(&[u32])) {
        let mut buf: Vec<u32> = Vec::new();
        let mut qi = 0usize;
        while qi < query.len() {
            let q = query[qi];
            let b = self.first.partition_point(|&x| x <= q);
            if b == 0 {
                // below the smallest hash in the index
                match self.first.first() {
                    Some(&f0) => {
                        while qi < query.len() && query[qi] < f0 {
                            qi += 1;
                        }
                        continue;
                    }
                    None => return,
                }
            }
            let b = b - 1;
            let next_first = self.first.get(b + 1).copied().unwrap_or(u64::MAX);
            let mut pos = self.off[b] as usize;
            let mut h = self.first[b];
            for i in 0..self.entries_in(b) {
                if i > 0 {
                    match get_varint(&self.data, &mut pos) {
                        Some(d) => h = h.saturating_add(d),
                        None => return,
                    }
                }
                let Some(len) = get_varint(&self.data, &mut pos) else {
                    return;
                };
                while qi < query.len() && query[qi] < h {
                    qi += 1;
                }
                let want = qi < query.len() && query[qi] == h;
                if want {
                    buf.clear();
                }
                let mut g = 0u32;
                for j in 0..len {
                    let Some(v) = get_varint(&self.data, &mut pos) else {
                        return;
                    };
                    if want {
                        g = if j == 0 {
                            v as u32
                        } else {
                            g.wrapping_add(v as u32)
                        };
                        buf.push(g);
                    }
                }
                if want {
                    f(&buf);
                    qi += 1;
                }
                if qi >= query.len() || query[qi] >= next_first {
                    break;
                }
            }
            // queries in this block's range that were not matched are not in the index
            while qi < query.len() && query[qi] < next_first {
                qi += 1;
            }
        }
    }

    /// Expand into the three arrays (for write, bench, cluster). None if the stream is corrupt.
    pub(crate) fn expand(&self) -> Option<(Vec<u64>, Vec<u64>, Vec<u32>)> {
        let mut hashes = Vec::with_capacity(self.n_hashes as usize);
        let mut offsets = Vec::with_capacity(self.n_hashes as usize + 1);
        let mut postings: Vec<u32> = Vec::with_capacity(self.n_postings as usize);
        offsets.push(0u64);
        for b in 0..self.first.len() {
            let mut pos = self.off[b] as usize;
            let mut h = self.first[b];
            for i in 0..self.entries_in(b) {
                if i > 0 {
                    h = h.checked_add(get_varint(&self.data, &mut pos)?)?;
                }
                hashes.push(h);
                let len = get_varint(&self.data, &mut pos)?;
                let mut g = 0u32;
                for j in 0..len {
                    let v = get_varint(&self.data, &mut pos)? as u32;
                    g = if j == 0 { v } else { g.checked_add(v)? };
                    postings.push(g);
                }
                offsets.push(postings.len() as u64);
            }
        }
        Some((hashes, offsets, postings))
    }
}
