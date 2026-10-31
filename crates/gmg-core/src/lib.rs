//! Core of Guess My Genome (GmG): a sourmash-compatible FracMinHash for DNA, plus the
//! reference index, search and taxon-calling logic.
//!
//! Deliberately has no filesystem access and no heavy dependencies so the same code builds
//! natively, as a Python extension and for wasm. Sketches are interchangeable with sourmash
//! 4.9.x `k=<k>,scaled=<scaled>,seed=42` DNA sketches (verified against the `sourmash` crate
//! in `tests/oracle.rs`).

mod call;
#[cfg(feature = "json")]
pub mod classify;
pub mod coverage;
mod fastx;
mod index;
#[cfg(feature = "json")]
mod json;
mod packed;
mod sketch;
mod skfile;

pub use call::{guess, guess_from_hits, Guess, RankVote, Tentative, Thresholds, TopHit, Withheld};
pub use fastx::{sketch_fastx, sketch_fastx_stats, SeqStats};
pub use index::{Explained, GenomeMeta, Hit, IndexBuilder, RefIndex, RANKS};
#[cfg(feature = "json")]
pub use json::guess_json;
pub use sketch::{containment, downsample, jaccard, max_hash_for_scaled, FracSketch, DEFAULT_SEED};
pub use skfile::SketchFile;
