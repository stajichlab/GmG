#![no_main]
//! Any bytes as a .gmgs sketch file.
use gmg_core::SketchFile;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = SketchFile::read(&mut &data[..]);
});
