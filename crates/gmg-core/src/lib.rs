//! Core of Guess My Genome (GmG): a sourmash-compatible FracMinHash for DNA.
//!
//! Deliberately has no I/O and no heavy dependencies so the same code builds natively, as a
//! Python extension and for wasm. Sketches are interchangeable with sourmash 4.9.x
//! `k=<k>,scaled=<scaled>,seed=42` DNA sketches (verified against the `sourmash` crate in
//! `tests/oracle.rs`).

mod sketch;

pub use sketch::{containment, jaccard, max_hash_for_scaled, FracSketch, DEFAULT_SEED};
