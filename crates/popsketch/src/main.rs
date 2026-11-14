use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use popsketch::pipeline::{run, FlagMode, Options, OrdKind, ReportFormat};
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
enum ReportArg {
    Html,
    Md,
    Pdf,
    All,
}

#[derive(Clone, Copy, ValueEnum)]
enum OrdArg {
    HashPca,
    Pcoa,
}

// the subcommands are parsed once at start-up; the size difference between them does not matter
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Cmd {
    /// Sketch the files of ONE sample into a folder of sketches (including the split halves), so the reads can be deleted:
    /// `popsketch run --from-sketches DIR` analyses the folder later
    #[command(
        after_long_help = "EXAMPLES\n  popsketch sketch --id SRR123 -o sketches/ SRR123_1.fastq.gz SRR123_2.fastq.gz\n  popsketch run --from-sketches sketches/ -o out/"
    )]
    Sketch {
        /// the sample's files (R1 and R2 are one sample)
        #[arg(required = true)]
        files: Vec<PathBuf>,
        #[arg(long)]
        id: String,
        /// folder of sketches (created; shared by all samples of one set, with the same parameters)
        #[arg(short, long)]
        out: PathBuf,
        /// reads or assembly
        #[arg(long, default_value = "reads", value_parser = ["reads", "assembly"])]
        kind: String,
        #[arg(short = 'k', long = "ksize", default_value_t = 21)]
        ksize: u32,
        #[arg(short = 's', long, default_value_t = 200)]
        scaled: u32,
        #[arg(long, default_value_t = 1_000_000_000)]
        max_bases: u64,
        #[arg(long)]
        r1_only: bool,
    },
    /// Check one sketched sample against the genome it should cover (was the run sampled evenly?). Prints one TSV line;
    /// exit 0 = ok, 75 = undersampled (the reads miss part of the genome or pile up on part of it: resample), 65 = not
    /// assessed (too little k-mer coverage, or not this genome: resampling will not help), 2 = error
    #[command(
        after_long_help = "EXAMPLES\n  popsketch check --id SRR123 -s sketches/ --genome-mb 12.25 --coverage 25\n\nThe thresholds (0.5 of the genome size, 2.5 times the expected k-mer coverage) come from one species and are UNMEASURED beyond it."
    )]
    Check {
        /// the sample's id (its sketch is `<ID>.gsk` in the folder)
        #[arg(long)]
        id: String,
        /// folder of sketches (as written by `popsketch sketch`)
        #[arg(short, long)]
        sketches: PathBuf,
        /// the genome size in Mb the sample should cover
        #[arg(long)]
        genome_mb: f64,
        /// the x the kept reads should give if sampled evenly (bases kept / genome size); omit to check the size only
        #[arg(long)]
        coverage: Option<f64>,
        /// fail when the estimated size is below this fraction of the genome
        #[arg(long, default_value_t = 0.5)]
        min_size_ratio: f64,
        /// fail when the k-mer coverage is above this many times the expected coverage
        #[arg(long, default_value_t = 2.5)]
        max_cov_ratio: f64,
        /// a failing sample with less k-mer coverage than this is not assessed (exit 65), not undersampled
        #[arg(long, default_value_t = 3.0)]
        min_kmer_cov: f64,
        /// a sample whose estimated size is more than this many times the genome is not assessed (exit 65): not this genome
        #[arg(long, default_value_t = 5.0)]
        max_size_ratio: f64,
        /// print the header line first
        #[arg(long)]
        header: bool,
    },
    /// Write a small simulated dataset with planted problems and a swapped label, to try `popsketch run` on
    #[command(
        after_long_help = "EXAMPLES\n  popsketch demo -o demo/\n  popsketch run demo/reads -o demo/out --labels demo/labels.tsv --scaled 5 --max-bases 0"
    )]
    Demo {
        /// folder to write (reads/ and labels.tsv are created in it)
        #[arg(short, long, default_value = "popsketch_demo")]
        out: PathBuf,
        /// bases per simulated genome
        #[arg(long, default_value_t = 120_000)]
        genome_len: usize,
        /// k-mer coverage of the clean samples
        #[arg(long, default_value_t = 12.0)]
        coverage: f64,
    },
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
        #[arg(long, default_value = "auto", value_parser = ["auto", "assembly", "reads"])]
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
        /// stop reading a sample once its k-mer coverage (from the reads so far) reaches this; --max-bases stays as a ceiling
        #[arg(long)]
        target_cov: Option<f64>,
        /// size flag when |robust z| of the estimated size exceeds this
        #[arg(long, default_value_t = 5.0)]
        size_z: f64,
        /// floor of the size spread, as a fraction of the cohort median size
        #[arg(long, default_value_t = 0.03)]
        size_floor: f64,
        /// abundance_shoulder flag when its robust z exceeds this
        #[arg(long, default_value_t = 5.0)]
        shoulder_z: f64,
        /// skip the ordination, DAPC and PC comparison (quick threshold sweeps)
        #[arg(long)]
        no_pca: bool,
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
        /// PCs used by DAPC on the primary label column (0: no DAPC)
        #[arg(long, default_value_t = 10)]
        dapc_pcs: usize,
        /// label column whose groups (for example the study) the permutation tests stay inside
        #[arg(long, value_name = "COLUMN")]
        strata: Option<String>,
        /// permutations for the PERMANOVA/PERMDISP test of each label column (0: no test)
        #[arg(long, default_value_t = 999)]
        permutations: usize,
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
        /// largest presence matrix (GiB) for the fast distance method; above it the slow pairwise merge is used
        #[arg(long, default_value_t = 4.0)]
        bitset_max_gib: f64,
        /// collapse near-identical samples to one representative each before the ordination, tree, clusters and label tests;
        /// the threshold of a pair is set by the samples' own split-half distances: `mean` of the two (default when given),
        /// `max` of the two, or the species `median`. Needs split-half sketches (see --floor-samples)
        #[arg(long, value_name = "RULE", num_args = 0..=1, default_missing_value = "mean", value_parser = ["mean", "max", "median"])]
        collapse: Option<String>,
        /// visit samples in a pseudo-random order when collapsing (tests how much the result depends on the order)
        #[arg(long, hide = true)]
        collapse_seed: Option<u64>,
        /// with --collapse: merge samples only when they carry the same labels in these label columns (comma-separated, or `all`),
        /// so no representative stands for samples with different labels
        #[arg(long, value_name = "COLUMNS", value_delimiter = ',')]
        collapse_within: Vec<String>,
        /// Newick tree of the same samples (for example from SNPs) to compare with: distance correlation, and whether each cluster is a clade
        #[arg(long)]
        compare_tree: Option<PathBuf>,
        /// PC scores of the same samples (for example from SNPs; id column then PC columns) to compare with the sketch PCA
        #[arg(long)]
        compare_pcs: Option<PathBuf>,
        /// reuse sketches from an earlier run in the same output folder (only sketching is resumed)
        #[arg(long)]
        resume: bool,
        /// take the samples from a folder of sketches made by `popsketch sketch` (the inputs are then not needed)
        #[arg(long)]
        from_sketches: Option<PathBuf>,
        /// report.html (default), report.md with figures/*.svg, report.pdf, or all three
        #[arg(long, value_enum, default_value = "html")]
        report_format: ReportArg,
        #[arg(short, long)]
        quiet: bool,
        #[arg(long)]
        threads: Option<usize>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Sketch {
            files,
            id,
            out,
            kind,
            ksize,
            scaled,
            max_bases,
            r1_only,
        } => {
            use popsketch::sketch_io::{cache_init, cache_write, sketch_sample, SketchParams};
            let params = SketchParams {
                ksize,
                scaled,
                max_bases,
                r1_only,
                ..SketchParams::default()
            };
            // refuses a folder made with other parameters (it would mix sketches that cannot be compared)
            cache_init(&out, &params)?;
            // halves always: the noise floor picks its samples when the whole set is known, which is after the reads are gone
            let s = sketch_sample(&id, &kind, &files, &params, true)?;
            cache_write(&out, &s, &params)?;
            eprintln!(
                "{id}: {} reads, {} bases, {} hashes",
                s.counts.reads,
                s.counts.bases,
                s.spectrum.hashes.len()
            );
            Ok(())
        }
        Cmd::Check {
            id,
            sketches,
            genome_mb,
            coverage,
            min_size_ratio,
            max_cov_ratio,
            min_kmer_cov,
            max_size_ratio,
            header,
        } => {
            use popsketch::check::{
                check, validate, CheckParams, Verdict, EXIT_ERROR, EXIT_NOT_ASSESSED,
                EXIT_UNDERSAMPLED, TSV_HEADER,
            };
            // every failure of `check` is exit 2, not anyhow's 1: a caller branches on 0, 75 and 2
            let fail = |msg: String| -> ! {
                eprintln!("popsketch check: {msg}");
                std::process::exit(EXIT_ERROR)
            };
            let params = CheckParams {
                min_size_ratio,
                max_cov_ratio,
                min_kmer_cov,
                max_size_ratio,
            };
            if let Err(e) = validate(genome_mb, coverage, &params) {
                fail(e);
            }
            let params_text =
                std::fs::read_to_string(sketches.join("params.txt")).unwrap_or_else(|e| {
                    fail(format!(
                        "{}: no readable params.txt, is it a folder of sketches? ({e})",
                        sketches.display()
                    ))
                });
            let (ksize, scaled) = popsketch::sketch_io::ksize_scaled_from_tag(&params_text)
                .unwrap_or_else(|| {
                    fail(format!(
                        "no k= and scaled= in {}",
                        sketches.join("params.txt").display()
                    ))
                });
            let sample = match popsketch::sketch_io::cache_read_checked(&sketches, &id, false) {
                Ok(Some(s)) => s,
                Ok(None) => fail(format!("no sketch for {id} in {}", sketches.display())),
                Err(e) => fail(format!("{e:#}")),
            };
            let r = check(&sample, scaled, ksize, genome_mb, coverage, &params);
            if header {
                println!("{TSV_HEADER}");
            }
            println!("{}", r.tsv());
            match r.verdict {
                Verdict::Ok => {}
                Verdict::Undersampled => std::process::exit(EXIT_UNDERSAMPLED),
                Verdict::NotAssessed => std::process::exit(EXIT_NOT_ASSESSED),
            }
            Ok(())
        }
        Cmd::Demo {
            out,
            genome_len,
            coverage,
        } => {
            let d = popsketch::demo::generate(&out, genome_len, coverage)?;
            eprintln!(
                "wrote {} samples to {} and the labels to {}\nplanted: low coverage lowcov0 lowcov1; another organism foreign; half the content half; a second genome mixed; the label of {} says group1 but its reads are group0\ntry: popsketch run {} -o {} --labels {} --scaled 5 --max-bases 0",
                d.group_of.len() + d.problems.len(),
                d.reads_dir.display(),
                d.labels.display(),
                d.swapped,
                d.reads_dir.display(),
                out.join("out").display(),
                d.labels.display()
            );
            Ok(())
        }
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
            target_cov,
            min_cov,
            size_z,
            size_floor,
            shoulder_z,
            no_pca,
            max_samples,
            flagged,
            index,
            ordination,
            components,
            dapc_pcs,
            strata,
            permutations,
            k,
            cut,
            min_floors,
            no_tree,
            root,
            tree_floor,
            floor_samples,
            bitset_max_gib,
            collapse,
            collapse_seed,
            collapse_within,
            compare_tree,
            compare_pcs,
            resume,
            from_sketches,
            report_format,
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
                target_cov,
                min_cov,
                size_z,
                size_floor,
                shoulder_z,
                no_pca,
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
                dapc_pcs,
                strata,
                permutations,
                k,
                cut,
                min_floors,
                tree: !no_tree,
                root,
                tree_floor,
                floor_samples,
                bitset_max_bytes: (bitset_max_gib.max(0.0) * (1u64 << 30) as f64) as usize,
                collapse: collapse
                    .as_deref()
                    .and_then(popsketch::collapse::Rule::parse),
                collapse_seed,
                collapse_within,
                compare_tree,
                compare_pcs,
                resume,
                from_sketches,
                report_format: match report_format {
                    ReportArg::Html => ReportFormat::Html,
                    ReportArg::Md => ReportFormat::Md,
                    ReportArg::Pdf => ReportFormat::Pdf,
                    ReportArg::All => ReportFormat::All,
                },
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
                let report = match o.report_format {
                    ReportFormat::Html => "report.html",
                    ReportFormat::Md => "report.md",
                    ReportFormat::Pdf => "report.pdf",
                    ReportFormat::All => "report.html, report.md and report.pdf",
                };
                eprintln!(
                    "{} samples, {} flagged, {} in the analysis; report: {}/{report}",
                    r.rows.len(),
                    flagged,
                    r.included.len(),
                    out.display()
                );
            }
            Ok(())
        }
    }
}
