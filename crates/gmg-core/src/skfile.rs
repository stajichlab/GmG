//! Single-sketch file (`.gsk`): little-endian, no compression.
//!
//! magic "GMGSKT01" | k u32 | scaled u32 | seed u64 | name_len u32 | name | n u64 | hashes u64*n |
//! has_abund u8 | abund u32*n (only if has_abund)

use std::io::{self, Read, Write};

const MAGIC: &[u8; 8] = b"GMGSKT01";

#[derive(Clone, Debug, PartialEq)]
pub struct SketchFile {
    pub ksize: u32,
    pub scaled: u32,
    pub seed: u64,
    pub name: String,
    pub hashes: Vec<u64>,
    pub abunds: Option<Vec<u32>>,
}

pub(crate) fn rd_u32<R: Read>(r: &mut R) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
pub(crate) fn rd_u64<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
/// Read exactly `len` bytes without trusting `len` for the allocation: a hostile header cannot make us allocate more
/// than the file actually holds.
/// k must fit a 64-bit k-mer and scaled must be at least 1.
pub(crate) fn check_params(ksize: u32, scaled: u32) -> io::Result<()> {
    if !(1..=32).contains(&ksize) || scaled == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid k-mer size or scaled value",
        ));
    }
    Ok(())
}

pub(crate) fn rd_bytes<R: Read>(r: &mut R, len: usize) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    r.take(len as u64).read_to_end(&mut buf)?;
    if buf.len() != len {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "file is truncated",
        ));
    }
    Ok(buf)
}
fn too_big() -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, "implausible length in header")
}
pub(crate) fn rd_u64_vec<R: Read>(r: &mut R, n: usize) -> io::Result<Vec<u64>> {
    let buf = rd_bytes(r, n.checked_mul(8).ok_or_else(too_big)?)?;
    Ok(buf
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect())
}
pub(crate) fn rd_u32_vec<R: Read>(r: &mut R, n: usize) -> io::Result<Vec<u32>> {
    let buf = rd_bytes(r, n.checked_mul(4).ok_or_else(too_big)?)?;
    Ok(buf
        .as_chunks::<4>()
        .0
        .iter()
        .map(|c| u32::from_le_bytes(*c))
        .collect())
}
pub(crate) fn wr_u64_slice<W: Write>(w: &mut W, v: &[u64]) -> io::Result<()> {
    let mut buf = Vec::with_capacity(v.len() * 8);
    for x in v {
        buf.extend_from_slice(&x.to_le_bytes());
    }
    w.write_all(&buf)
}
pub(crate) fn wr_u32_slice<W: Write>(w: &mut W, v: &[u32]) -> io::Result<()> {
    let mut buf = Vec::with_capacity(v.len() * 4);
    for x in v {
        buf.extend_from_slice(&x.to_le_bytes());
    }
    w.write_all(&buf)
}

impl SketchFile {
    pub fn write<W: Write>(&self, w: &mut W) -> io::Result<()> {
        w.write_all(MAGIC)?;
        w.write_all(&self.ksize.to_le_bytes())?;
        w.write_all(&self.scaled.to_le_bytes())?;
        w.write_all(&self.seed.to_le_bytes())?;
        w.write_all(&(self.name.len() as u32).to_le_bytes())?;
        w.write_all(self.name.as_bytes())?;
        w.write_all(&(self.hashes.len() as u64).to_le_bytes())?;
        wr_u64_slice(w, &self.hashes)?;
        match &self.abunds {
            Some(a) => {
                assert_eq!(a.len(), self.hashes.len());
                w.write_all(&[1])?;
                wr_u32_slice(w, a)
            }
            None => w.write_all(&[0]),
        }
    }

    pub fn read<R: Read>(r: &mut R) -> io::Result<Self> {
        let mut m = [0u8; 8];
        r.read_exact(&mut m)?;
        if &m != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a GmG sketch file",
            ));
        }
        let ksize = rd_u32(r)?;
        let scaled = rd_u32(r)?;
        let seed = rd_u64(r)?;
        check_params(ksize, scaled)?;
        let nl = rd_u32(r)? as usize;
        let nb = rd_bytes(r, nl)?;
        let name =
            String::from_utf8(nb).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let n = rd_u64(r)? as usize;
        let hashes = rd_u64_vec(r, n)?;
        let mut f = [0u8; 1];
        r.read_exact(&mut f)?;
        let abunds = if f[0] == 1 {
            Some(rd_u32_vec(r, n)?)
        } else {
            None
        };
        Ok(Self {
            ksize,
            scaled,
            seed,
            name,
            hashes,
            abunds,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let s = SketchFile {
            ksize: 21,
            scaled: 1000,
            seed: 42,
            name: "GCA_x".into(),
            hashes: vec![1, 5, 9],
            abunds: Some(vec![3, 1, 2]),
        };
        let mut buf = Vec::new();
        s.write(&mut buf).unwrap();
        assert_eq!(SketchFile::read(&mut &buf[..]).unwrap(), s);
    }
}
