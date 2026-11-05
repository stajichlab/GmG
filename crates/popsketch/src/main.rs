use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use popsketch::pipeline::{run, FlagMode, Options, OrdKind};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "popsketch",
    version,
    about = "Audit a set of read files and draw coarse population structure from FracMinHash sketches"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Clone, Copy, ValueEnum)]
enum FlaggedArg {
    /// leave flagged samples out of distances, PCA, clusters and tree (default)
    Exclude,
    /// keep them and draw them with their own marker
    Mark,
    /// remove contamination with the cohort vote and put them back (other flagged samples stay out)
    Decontaminate,
}

#[derive(Clone, Copy, ValueEnum)]
enum OrdArg {
    HashPca,
    Pcoa,
}

#[derive(Subcommand)]
enum Cmd {
    /// Sketch the samples, audit them, and write distances, PCA, clusters, a tree and a report
    #[command(
        after_long_help = "EXAMPLES\n  popsketch run reads/ -o out/ --labels labels.tsv\n  popsketch run --sample-sheet samples.tsv -o out/ --flagged mark --index bfd.gmgi"
    )]
    Run {
        /// FASTQ/FASTA files or folders (searched recursively); R1/R2 and lanes of a sample are merged
        inputs: Vec<PathBuf>,
        /// TSV with columns id and path (rows with the same id are one sample)
        #[arg(long)]
        sample_sheet: Option<PathBuf>,
        /// label table: first column the sample id, then one column per label set (the first is the primary one); blank = unlabelled
        #[arg(long)]
        labels: Option<PathBuf>,
        #[arg(short, long, default_value = "popsketch_out")]
        out: PathBuf,
        /// auto, assembly or reads
        #[arg(long, default_value = "auto")]
        kind: String,
        /// do not merge the read files of a sample
        #[arg(long)]
        no_group: bool,
        #[arg(short = 'k', long, default_value_t = 21)]
        ksize: u32,
        #[arg(short, long, default_value_t = 200)]
        scaled: u32,
        /// bases read per sample, shared by its files; 0 = no limit
        #[arg(long, default_value_t = 1_000_000_000)]
        max_bases: u64,
        /// read only the first file of each sample (fast mode)
        #[arg(long)]
        r1_only: bool,
        /// samples with k-mer coverage below this are flagged low_coverage
        #[arg(long, default_value_t = 4.0)]
        min_cov: f64,
        /// refuse to run on more samples than this (exact all-pairs distances grow with the square)
        #[arg(long, default_value_t = 1500)]
        max_samples: usize,
        /// what to do with flagged samples
        #[arg(long, value_enum, default_value = "exclude")]
        flagged: FlaggedArg,
        /// reference index (.gmgi): names the organism and finds mixtures (optional)
        #[arg(long)]
        index: Option<PathBuf>,
        #[arg(long, value_enum, default_value = "hash-pca")]
        ordination: OrdArg,
        /// number of components to compute
        #[arg(long, default_value_t = 5)]
        components: usize,
        /// force this many clusters
        #[arg(long)]
        k: Option<usize>,
        /// cut the clustering at this distance (1 - ANI)
        #[arg(long)]
        cut: Option<f64>,
        /// call the largest jump a clear cut only if it is at least this many noise floors
        #[arg(long, default_value_t = 3.0)]
        min_floors: f64,
        /// do not build a tree
        #[arg(long)]
        no_tree: bool,
        /// root the tree on this sample
        #[arg(long)]
        root: Option<String>,
        /// collapse tree branches shorter than this (default: the noise floor)
        #[arg(long)]
        tree_floor: Option<f64>,
        /// samples that get split-half sketches to measure the noise floor (0: fixed fallback)
        #[arg(long, default_value_t = 20)]
        floor_samples: usize,
        /// reuse sketches from an earlier run in the same output folder (only sketching is resumed)
        #[arg(long)]
        resume: bool,
        #[arg(short, long)]
        quiet: bool,
        #[arg(long)]
        threads: Option<usize>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Run {
            inputs,
            sample_sheet,
            labels,
            out,
            kind,
            no_group,
            ksize,
            scaled,
            max_bases,
            r1_only,
            min_cov,
            max_samples,
            flagged,
            index,
            ordination,
            components,
            k,
            cut,
            min_floors,
            no_tree,
            root,
            tree_floor,
            floor_samples,
            resume,
            quiet,
            threads,
        } => {
            let o = Options {
                inputs,
                sample_sheet,
                labels,
                out: out.clone(),
                kind,
                no_group,
                ksize,
                scaled,
                max_bases,
                r1_only,
                min_cov,
                max_samples,
                flagged: match flagged {
                    FlaggedArg::Exclude => FlagMode::Exclude,
                    FlaggedArg::Mark => FlagMode::Mark,
                    FlaggedArg::Decontaminate => FlagMode::Decontaminate,
                },
                index,
                ordination: match ordination {
                    OrdArg::HashPca => OrdKind::HashPca,
                    OrdArg::Pcoa => OrdKind::Pcoa,
                },
                components,
                k,
                cut,
                min_floors,
                tree: !no_tree,
                root,
                tree_floor,
                floor_samples,
                resume,
                quiet,
                threads,
            };
            let r = run(&o)?;
            if !o.quiet {
                let flagged = r
                    .rows
                    .iter()
                    .filter(|x| x.qc.as_ref().is_some_and(|q| q.flagged()))
                    .count();
                eprintln!(
                    "{} samples, {} flagged, {} in the analysis; report: {}",
                    r.rows.len(),
                    flagged,
                    r.included.len(),
                    out.join("report.html").display()
                );
            }
            Ok(())
        }
    }
}
