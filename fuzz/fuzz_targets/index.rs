#![no_main]
//! Any bytes as an index file: loading (streamed and owned/packed) and a search must end in an error or a result, never a panic.
use gmg_core::{guess, RefIndex, Thresholds};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let q: Vec<u64> = (0..50).map(|x| x * 7).collect();
    if let Ok(ix) = RefIndex::read(&mut &data[..]) {
        let _ = guess(&ix, &q, &|_| false, &Thresholds::default());
    }
    if let Ok(ix) = RefIndex::read_owned(data.to_vec()) {
        let _ = guess(&ix, &q, &|_| false, &Thresholds::default());
    }
});
