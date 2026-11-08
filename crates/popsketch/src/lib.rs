//! popsketch: audit a set of read files and draw coarse population structure from FracMinHash sketches.
//! See `docs/popstruct-plan.md` in the repository for the design and the evidence behind it.

pub mod advisor;
pub mod cluster;
pub mod compare;
pub mod dapc;
#[cfg(feature = "native")]
pub mod demo;
pub mod dist;
pub mod floor;
pub mod labels;
pub mod ordination;
pub mod par;
pub mod permanova;
pub mod pipeline;
pub mod qc;
pub mod report;
pub mod sketch_io;
pub mod stats;
pub mod tree;
