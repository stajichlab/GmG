//! Plain-language failure codes for inputs `gmg run` cannot process (#55).
use anyhow::Error;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::Path;

pub struct InputError {
    pub code: &'static str,
    /// one sentence, no paths
    pub message: String,
    /// what to do about it
    pub action: &'static str,
}

fn e(code: &'static str, message: &str, action: &'static str) -> InputError {
    InputError {
        code,
        message: message.to_string(),
        action,
    }
}

/// Look at the file itself (size, gzip integrity, first character) and at the parser's error to say what is wrong.
/// `bases` is the number of bases read when the file parsed but no usable k-mer window was found.
pub fn classify(path: &Path, size: u64, err: &Error) -> InputError {
    let text = format!("{err:#}");
    if text.contains("internal error") {
        return e(
            "internal_error",
            "GmG hit an unexpected problem with this file.",
            "Please report it (https://github.com/stajichlab/GmG/issues) with the file if you can share it.",
        );
    }
    if size == 0 {
        return e(
            "empty_file",
            "The file is empty (0 bytes).",
            "Re-download or re-copy it.",
        );
    }
    let mut f = match File::open(path) {
        Ok(f) => f,
        Err(_) => {
            return e(
                "unreadable",
                "The file could not be opened.",
                "Check that it exists and that you have read permission.",
            )
        }
    };
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic).unwrap_or(0);
    let gz = n == 2 && magic == [0x1f, 0x8b];
    let mut head = Vec::new();
    let rd = File::open(path).and_then(|f| {
        if gz {
            MultiGzDecoder::new(f).take(4096).read_to_end(&mut head)
        } else {
            f.take(4096).read_to_end(&mut head)
        }
    });
    if gz && rd.is_err() {
        return e(
            "bad_gzip",
            "The gzip file is truncated or corrupt.",
            "Re-download it, or check it with `gzip -t`.",
        );
    }
    if rd.is_err() {
        return e(
            "unreadable",
            "The file could not be read.",
            "Check the disk and file permissions.",
        );
    }
    match head.iter().find(|b| !b.is_ascii_whitespace()) {
        None => {
            return e(
                "no_sequence_records",
                "The file has no sequence records (only blank content).",
                "Check that it is the right file.",
            )
        }
        Some(b'>') => {
            // the whole (small) file is a header with no sequence lines
            let has_seq = head.split(|b| *b == b'\n').any(|l| {
                !l.is_empty() && !l.starts_with(b">") && !l.iter().all(|b| b.is_ascii_whitespace())
            });
            if !has_seq && head.len() < 4096 {
                return e(
                    "no_sequence_records",
                    "The file has a header but no sequence.",
                    "Check that the download or conversion finished.",
                );
            }
        }
        Some(b'@') => {}
        Some(_) => return e(
            "not_sequence_file",
            "The file does not look like FASTA or FASTQ.",
            "GmG reads FASTA (starting with '>') or FASTQ (starting with '@'), optionally gzipped.",
        ),
    }
    if text.contains("no sequence records") {
        return e(
            "no_sequence_records",
            "The file has a header but no sequence.",
            "Check that the download or conversion finished.",
        );
    }
    if text.contains("too short") {
        return e(
            "too_short",
            "The input has fewer bases than one k-mer window.",
            "Provide a real genome or read set.",
        );
    }
    e(
        "parse_error",
        "The file could not be parsed as FASTA or FASTQ.",
        "Check that it is complete and well-formed.",
    )
}
