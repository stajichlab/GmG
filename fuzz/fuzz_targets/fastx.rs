#![no_main]
//! Any bytes as FASTA/FASTQ text (the same entry points the web app calls).
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = gmg_core::sketch_fastx(data, 21, 10, 42, 1, 0);
    let _ = gmg_core::sketch_fastx_stats(data, 21, 10, false, 0);
    let _ = gmg_core::sketch_fastx_stats(data, 21, 10, true, 5);
});
