//! Input discovery shared by `gmg run` and tools built on the same sketching: find FASTA/FASTQ files, name samples,
//! merge the read files of one sample (R1/R2, lanes) and read a sample sheet.
//!
//! Moved out of the `gmg` binary crate (`gmg-cli/src/run.rs`) without a change of behavior, so another crate can reuse it.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

/// the path that means standard input (gzip and FASTA/FASTQ are detected from the bytes)
pub fn is_stdin(p: &Path) -> bool {
    p.as_os_str() == "-"
}

/// A file name without a compression suffix and a sequence-file extension (also `.gsk`).
pub fn strip_ext(name: &str) -> String {
    let name = COMPRESSED
        .iter()
        .find_map(|c| name.strip_suffix(c))
        .unwrap_or(name);
    for e in [
        ".fa", ".fasta", ".fna", ".fsa", ".seq", ".dna", ".fastq", ".fq", ".gsk",
    ] {
        if let Some(s) = name.strip_suffix(e) {
            return s.to_string();
        }
    }
    name.to_string()
}

pub struct Input {
    /// one file, or the files of one read sample (R1/R2, lanes) merged into one sketch; `path` is the first
    pub path: PathBuf,
    pub files: Vec<PathBuf>,
    pub id: String,
    pub kind: &'static str, // "assembly" | "reads"
    /// from a sample sheet: the label the user expects (carried to results.tsv for validation)
    pub expected: Option<String>,
    /// from a sample sheet: the reference genome id of this query, for `--exclude-level`
    pub ref_id: Option<String>,
    /// named by a sample sheet: never merged automatically
    pub explicit: bool,
}

pub const FASTA_EXT: [&str; 8] = ["fa", "fna", "fasta", "fas", "ffn", "fsa", "seq", "dna"];
pub const FASTQ_EXT: [&str; 2] = ["fq", "fastq"];

/// Compression suffixes needletail reads (it detects the format from the first bytes, the suffix only selects the file).
pub const COMPRESSED: [&str; 4] = [".gz", ".zst", ".bz2", ".xz"];

/// File extension without a trailing compression suffix, lowercased.
fn base_ext(p: &Path) -> Option<String> {
    let name = p.file_name()?.to_string_lossy().to_lowercase();
    let name = COMPRESSED
        .iter()
        .find_map(|c| name.strip_suffix(c))
        .unwrap_or(&name)
        .to_string();
    name.rsplit_once('.').map(|(_, e)| e.to_string())
}

pub fn kind_of(p: &Path, forced: &str) -> Option<&'static str> {
    let ext = base_ext(p)?;
    let auto = if FASTA_EXT.contains(&ext.as_str()) {
        Some("assembly")
    } else if FASTQ_EXT.contains(&ext.as_str()) {
        Some("reads")
    } else {
        None
    };
    match forced {
        "assembly" => auto.map(|_| "assembly"),
        "reads" => auto.map(|_| "reads"),
        _ => auto,
    }
}

fn walk(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .with_context(|| format!("reading directory {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            walk(&p, out)?;
        } else {
            out.push(p);
        }
    }
    Ok(())
}

/// Sample key of a read file: the stem without the read-direction (`_R1`, `_R2`, `_R1_001`, `.1`/`_1`/`_2`) and lane (`_L001`) parts.
pub fn sample_key(stem: &str) -> String {
    let mut parts: Vec<&str> = stem.split('_').collect();
    // drop lane tokens such as L001
    parts.retain(|t| {
        !(t.len() == 4 && t.starts_with('L') && t[1..].chars().all(|c| c.is_ascii_digit()))
    });
    // a trailing 3-digit chunk after R1/R2 (Illumina `_R1_001`)
    if parts.len() >= 2 {
        let last = parts[parts.len() - 1];
        let prev = parts[parts.len() - 2];
        if last.len() == 3
            && last.chars().all(|c| c.is_ascii_digit())
            && matches!(prev, "R1" | "R2")
        {
            parts.truncate(parts.len() - 2);
            return parts.join("_");
        }
    }
    if let Some(last) = parts.last().copied() {
        if matches!(last, "R1" | "R2" | "1" | "2") && parts.len() > 1 {
            parts.pop();
        } else if let Some(b) = last.strip_suffix(".1").or_else(|| last.strip_suffix(".2")) {
            let n = parts.len();
            parts[n - 1] = b;
        }
    }
    parts.join("_")
}

/// Merge read files of one sample into a single input (same directory, same key, kind reads). A lone file keeps its own id.
pub fn group_reads(inputs: Vec<Input>) -> Vec<Input> {
    let mut order: Vec<(PathBuf, String)> = Vec::new();
    let mut groups: HashMap<(PathBuf, String), Vec<Input>> = HashMap::new();
    let mut out: Vec<Input> = Vec::new();
    for i in inputs {
        if i.kind != "reads" || i.explicit {
            out.push(i);
            continue;
        }
        let key = (
            i.path.parent().map(|p| p.to_path_buf()).unwrap_or_default(),
            sample_key(&i.id),
        );
        if !groups.contains_key(&key) {
            order.push(key.clone());
        }
        groups.entry(key).or_default().push(i);
    }
    for key in order {
        let mut g = groups.remove(&key).unwrap();
        if g.len() == 1 {
            out.push(g.pop().unwrap());
        } else {
            let files: Vec<PathBuf> = g.iter().map(|i| i.path.clone()).collect();
            out.push(Input {
                path: files[0].clone(),
                files,
                id: key.1.clone(),
                kind: "reads",
                expected: None,
                ref_id: None,
                explicit: false,
            });
        }
    }
    out.sort_by(|a, b| a.id.cmp(&b.id));
    out
}

/// Sample sheet: TSV with a header and columns `id`, `path` (relative to the sheet) and optionally `expected`.
/// Rows with the same id are the files (R1/R2, lanes) of one sample.
pub fn read_sheet(sheet: &Path, forced: &str) -> Result<Vec<Input>> {
    let text =
        fs::read_to_string(sheet).with_context(|| "reading the --sample-sheet file".to_string())?;
    let base = sheet.parent().map(Path::to_path_buf).unwrap_or_default();
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let header: Vec<&str> = lines
        .next()
        .unwrap_or("")
        .split('\t')
        .map(str::trim)
        .collect();
    let col = |n: &str| header.iter().position(|h| h.eq_ignore_ascii_case(n));
    let (Some(ci), Some(cp)) = (col("id"), col("path")) else {
        bail!("sample sheet needs a header with columns 'id' and 'path' (and optionally 'expected'), separated by tabs");
    };
    let ce = col("expected");
    let cr = col("ref");
    let mut out: Vec<Input> = Vec::new();
    for (n, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let (Some(id), Some(path)) = (
            f.get(ci).filter(|x| !x.is_empty()),
            f.get(cp).filter(|x| !x.is_empty()),
        ) else {
            bail!("sample sheet row {}: id and path are required", n + 2);
        };
        let path = {
            let p = PathBuf::from(path);
            if p.is_absolute() {
                p
            } else {
                base.join(p)
            }
        };
        if !path.is_file() {
            bail!(
                "sample sheet row {}: file not found: {}",
                n + 2,
                path.display()
            );
        }
        let Some(kind) = kind_of(&path, forced) else {
            bail!(
                "sample sheet row {}: {} does not look like FASTA/FASTQ",
                n + 2,
                path.display()
            );
        };
        let expected = ce
            .and_then(|c| f.get(c))
            .filter(|x| !x.is_empty())
            .map(|x| x.to_string());
        let ref_id = cr
            .and_then(|c| f.get(c))
            .filter(|x| !x.is_empty())
            .map(|x| x.to_string());
        match out.iter_mut().find(|i| i.id == *id) {
            Some(i) => {
                if i.kind != kind {
                    bail!("sample sheet: '{id}' mixes assemblies and reads");
                }
                i.files.push(path);
            }
            None => out.push(Input {
                files: vec![path.clone()],
                path,
                id: id.to_string(),
                kind,
                expected,
                ref_id,
                explicit: true,
            }),
        }
    }
    Ok(out)
}

/// What to look at: files and directories (searched recursively), an optional list file, an optional sample sheet.
pub struct DiscoverArgs<'a> {
    pub inputs: &'a [PathBuf],
    /// file with one input path per line
    pub list: Option<&'a Path>,
    /// merge read files of one sample (R1/R2, lanes)
    pub group: bool,
    pub quiet: bool,
    /// TSV naming inputs (id, path, expected, ref)
    pub sample_sheet: Option<&'a Path>,
    /// `auto`, `assembly` or `reads`
    pub kind: &'a str,
}

pub fn discover(args: &DiscoverArgs) -> Result<Vec<Input>> {
    let mut files = Vec::new();
    let mut given: Vec<PathBuf> = args.inputs.to_vec();
    if let Some(l) = args.list {
        let text = fs::read_to_string(l).with_context(|| "reading the --list file".to_string())?;
        given.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(PathBuf::from),
        );
    }
    let mut stdin_given = false;
    for p in &given {
        if is_stdin(p) {
            if stdin_given {
                bail!("standard input (-) can be given only once");
            }
            stdin_given = true;
            files.push(p.clone());
        } else if p.is_dir() {
            walk(p, &mut files)?;
        } else if p.is_file() {
            files.push(p.clone());
        } else {
            bail!("input not found: {}", p.display());
        }
    }
    let mut inputs: Vec<Input> = Vec::new();
    let mut seen: HashMap<String, PathBuf> = HashMap::new();
    for p in files {
        // directories are filtered by extension; explicitly named files must also look like sequence files;
        // standard input has no name: `--kind reads` says it is FASTQ/reads, anything else is an assembly (FASTQ is refused later)
        let kind = if is_stdin(&p) {
            if args.kind == "reads" {
                "reads"
            } else {
                "assembly"
            }
        } else if let Some(k) = kind_of(&p, args.kind) {
            k
        } else {
            continue;
        };
        let id = if is_stdin(&p) {
            "stdin".to_string()
        } else {
            strip_ext(
                &p.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            )
        };
        if let Some(prev) = seen.insert(id.clone(), p.clone()) {
            bail!(
                "two inputs would get the id '{id}': {} and {}; rename one",
                prev.display(),
                p.display()
            );
        }
        inputs.push(Input {
            files: vec![p.clone()],
            path: p,
            id,
            kind,
            expected: None,
            ref_id: None,
            explicit: false,
        });
    }
    if args.group {
        inputs = group_reads(inputs);
        let merged: Vec<&Input> = inputs.iter().filter(|i| i.files.len() > 1).collect();
        if !merged.is_empty() && !args.quiet {
            eprintln!(
                "merged {} read files into {} samples (R1/R2 and lanes; --no-group keeps them separate)",
                merged.iter().map(|i| i.files.len()).sum::<usize>(),
                merged.len()
            );
        }
    }
    if let Some(sheet) = args.sample_sheet {
        for e in read_sheet(sheet, args.kind)? {
            if let Some(prev) = inputs.iter().find(|i| i.id == e.id) {
                bail!(
                    "sample sheet id '{}' is also used by {}; ids must be unique",
                    e.id,
                    prev.path.display()
                );
            }
            inputs.push(e);
        }
    }
    if inputs.is_empty() {
        bail!("no FASTA/FASTQ files found in the given inputs (extensions: fa fna fasta fas ffn fsa seq dna fq fastq, optionally .gz .zst .bz2 .xz)");
    }
    Ok(inputs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_key_drops_read_direction_and_lane() {
        assert_eq!(sample_key("S1_R1_001"), "S1");
        assert_eq!(sample_key("S1_R2_001"), "S1");
        assert_eq!(sample_key("S1_S5_L001_R1_001"), "S1_S5");
        assert_eq!(sample_key("run_1"), "run");
        assert_eq!(sample_key("run_2"), "run");
        assert_eq!(sample_key("SRR123.1"), "SRR123");
        assert_eq!(sample_key("SRR123_1"), "SRR123");
        assert_eq!(sample_key("sample"), "sample");
    }

    #[test]
    fn strip_ext_handles_compression() {
        assert_eq!(strip_ext("a.fa.gz"), "a");
        assert_eq!(strip_ext("b_R1.fastq.zst"), "b_R1");
        assert_eq!(strip_ext("c.gsk"), "c");
        assert_eq!(strip_ext("d.txt"), "d.txt");
    }

    #[test]
    fn kind_follows_the_extension() {
        assert_eq!(kind_of(Path::new("x.fna.gz"), "auto"), Some("assembly"));
        assert_eq!(kind_of(Path::new("x.FASTQ"), "auto"), Some("reads"));
        assert_eq!(kind_of(Path::new("x.txt"), "auto"), None);
        assert_eq!(kind_of(Path::new("x.fq"), "assembly"), Some("assembly"));
    }

    #[test]
    fn read_pairs_and_lanes_merge_into_one_sample() {
        let dir = std::env::temp_dir().join(format!("gmg_input_test_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        for f in [
            "A_S1_L001_R1_001.fastq.gz",
            "A_S1_L001_R2_001.fastq.gz",
            "A_S1_L002_R1_001.fastq.gz",
            "A_S1_L002_R2_001.fastq.gz",
            "solo.fastq",
            "g.fa",
        ] {
            fs::write(dir.join(f), b"").unwrap();
        }
        let inputs = vec![dir.clone()];
        let args = DiscoverArgs {
            inputs: &inputs,
            list: None,
            group: true,
            quiet: true,
            sample_sheet: None,
            kind: "auto",
        };
        let got = discover(&args).unwrap();
        let ids: Vec<(&str, usize, &str)> = got
            .iter()
            .map(|i| (i.id.as_str(), i.files.len(), i.kind))
            .collect();
        assert_eq!(
            ids,
            vec![
                ("A_S1", 4, "reads"),
                ("g", 1, "assembly"),
                ("solo", 1, "reads")
            ]
        );
        let args = DiscoverArgs {
            group: false,
            ..args
        };
        assert_eq!(discover(&args).unwrap().len(), 6);
        fs::remove_dir_all(&dir).ok();
    }
}
