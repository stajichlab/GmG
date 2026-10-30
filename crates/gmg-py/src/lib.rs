//! Python bindings: the same `gmg-core` code as `gmg run` / the web app, so a record from `Index.classify` equals the
//! corresponding line of `results.jsonl` (the tests in `tests/` check it against the CLI).
use flate2::read::MultiGzDecoder;
use gmg_core::classify::{
    bytes_error, classify, error_record, hits_tsv, results_tsv, Context, Item, SCHEMA_VERSION,
};
use gmg_core::{guess, guess_json, sketch_fastx, sketch_fastx_stats, RefIndex, Thresholds};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::io::Read;

fn to_py(py: Python<'_>, json_text: &str) -> PyResult<PyObject> {
    Ok(py
        .import("json")?
        .call_method1("loads", (json_text,))?
        .unbind())
}

/// Bytes of a file, decompressed when they start with the gzip magic.
fn read_maybe_gz(path: &str) -> PyResult<Vec<u8>> {
    let raw = std::fs::read(path).map_err(|e| PyIOError::new_err(format!("{path}: {e}")))?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .map_err(|e| PyIOError::new_err(format!("{path}: bad gzip: {e}")))?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

/// A reference index with the thresholds and context needed to classify.
#[pyclass]
struct Index {
    ix: RefIndex,
    th: Thresholds,
    ctx: Context,
    th_name: String,
}

#[pymethods]
impl Index {
    /// `Index(path)`: load an index file written by `gmg build` / `gmg convert` (any format version).
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let bytes = std::fs::read(path).map_err(|e| PyIOError::new_err(format!("{path}: {e}")))?;
        let ix = RefIndex::read_owned(bytes)
            .map_err(|e| PyValueError::new_err(format!("{path}: not a GmG index: {e}")))?;
        let ctx = Context::new(&ix);
        Ok(Index {
            ix,
            th: Thresholds::default(),
            ctx,
            th_name: "known".into(),
        })
    }

    /// Choose a built-in threshold set: "known" (default), "permissive" or "novel".
    fn use_thresholds(&mut self, name: &str) -> PyResult<()> {
        let text = match name {
            "known" => None,
            "permissive" => Some(include_str!(
                "../../../data/thresholds/bfd_s10000.permissive.json"
            )),
            "novel" => Some(include_str!(
                "../../../data/thresholds/bfd_s10000.novel.json"
            )),
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown threshold set '{other}' (known, permissive, novel)"
                )))
            }
        };
        self.th = match text {
            None => Thresholds::default(),
            Some(t) => Thresholds::from_json_str(t).map_err(PyValueError::new_err)?,
        };
        self.th_name = name.to_string();
        Ok(())
    }

    #[getter]
    fn ksize(&self) -> u32 {
        self.ix.ksize
    }
    #[getter]
    fn scaled(&self) -> u32 {
        self.ix.scaled
    }
    #[getter]
    fn n_genomes(&self) -> usize {
        self.ix.genomes.len()
    }
    #[getter]
    fn thresholds(&self) -> String {
        self.th_name.clone()
    }

    /// The manifest of the index as a dict (None for indexes built before manifests).
    fn manifest(&self, py: Python<'_>) -> PyResult<PyObject> {
        match &self.ix.manifest {
            Some(m) => to_py(py, m),
            None => Ok(py.None()),
        }
    }

    /// Sketch FASTA/FASTQ text (bytes, plain) at this index's k and scaled: sorted unique hashes.
    #[pyo3(signature = (data, min_abund=1, max_records=0))]
    fn sketch(&self, data: &[u8], min_abund: u32, max_records: u64) -> Vec<u64> {
        sketch_fastx(
            data,
            self.ix.ksize as usize,
            self.ix.scaled,
            self.ix.seed,
            min_abund,
            max_records,
        )
    }

    /// The JSON of `gmg guess` (as a dict) from hashes made with `sketch`.
    #[pyo3(signature = (hashes, top=5))]
    fn guess_hashes(&self, py: Python<'_>, hashes: Vec<u64>, top: usize) -> PyResult<PyObject> {
        let mut h = hashes;
        h.sort_unstable();
        h.dedup();
        let g = guess(&self.ix, &h, &|_| false, &self.th);
        to_py(py, &guess_json(&self.ix, &g, "", top).to_string())
    }

    /// Classify a FASTA/FASTQ file (plain or gzip) and return its record: the line `gmg run` writes to `results.jsonl`
    /// (a failure is a record with status "error", a code and a message, not an exception). `kind`: "assembly", "reads"
    /// or "auto" (FASTQ when the first byte is "@"); `min_abund` 0 = default (1 for assemblies, 2 for reads);
    /// `max_records` 0 = all.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (path, kind="auto", min_abund=0, max_records=0, top=5, id=None))]
    fn classify(
        &self,
        py: Python<'_>,
        path: &str,
        kind: &str,
        min_abund: u32,
        max_records: u64,
        top: usize,
        id: Option<&str>,
    ) -> PyResult<PyObject> {
        let bytes = read_maybe_gz(path)?;
        let name = std::path::Path::new(path)
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.to_string());
        let id_s = id.map(str::to_string).unwrap_or_else(|| {
            let base = name.strip_suffix(".gz").unwrap_or(&name);
            base.rsplit_once('.').map_or(base, |(a, _)| a).to_string()
        });
        let rec = self.classify_bytes(&id_s, path, &bytes, kind, min_abund, max_records, top);
        to_py(py, &rec)
    }

    /// Same, for sequence text already in memory (bytes, not gzipped).
    #[allow(clippy::too_many_arguments)]
    #[pyo3(name = "classify_bytes", signature = (data, name="input", kind="auto", min_abund=0, max_records=0, top=5))]
    fn classify_bytes_py(
        &self,
        py: Python<'_>,
        data: &[u8],
        name: &str,
        kind: &str,
        min_abund: u32,
        max_records: u64,
        top: usize,
    ) -> PyResult<PyObject> {
        let rec = self.classify_bytes(name, name, data, kind, min_abund, max_records, top);
        to_py(py, &rec)
    }
}

impl Index {
    #[allow(clippy::too_many_arguments)]
    fn classify_bytes(
        &self,
        id: &str,
        path: &str,
        bytes: &[u8],
        kind: &str,
        min_abund: u32,
        max_records: u64,
        top: usize,
    ) -> String {
        let kind = match kind {
            "reads" => "reads",
            "assembly" => "assembly",
            _ => {
                if bytes.first() == Some(&b'@') {
                    "reads"
                } else {
                    "assembly"
                }
            }
        };
        let item = Item {
            id,
            path,
            size_bytes: bytes.len() as u64,
            kind,
        };
        let reads = kind == "reads";
        let min_abund = if min_abund == 0 {
            if reads {
                2
            } else {
                1
            }
        } else {
            min_abund
        };
        let fail = |reason: &str| {
            let (code, msg, action) = bytes_error(bytes, reason);
            error_record(&item, code, msg, action).to_string()
        };
        let Some((sk, st)) =
            sketch_fastx_stats(bytes, self.ix.ksize, self.ix.scaled, reads, max_records)
        else {
            return fail("");
        };
        match classify(
            &self.ix,
            &self.th,
            &self.ctx,
            &item,
            sk,
            &st,
            min_abund,
            top,
            &|_| false,
        ) {
            Ok((rec, _)) => rec.to_string(),
            Err(msg) => fail(&msg),
        }
    }
}

/// `results.tsv` text from a list of records (dicts from `Index.classify`).
#[pyfunction]
fn results_tsv_text(py: Python<'_>, records: PyObject) -> PyResult<String> {
    tsv(py, records, false)
}

/// `hits.tsv` text from a list of records.
#[pyfunction]
fn hits_tsv_text(py: Python<'_>, records: PyObject) -> PyResult<String> {
    tsv(py, records, true)
}

fn tsv(py: Python<'_>, records: PyObject, hits: bool) -> PyResult<String> {
    let text: String = py
        .import("json")?
        .call_method1("dumps", (records,))?
        .extract()?;
    let v: Vec<serde_json::Value> =
        serde_json::from_str(&text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (header, rows) = if hits { hits_tsv(&v) } else { results_tsv(&v) };
    let mut out = header.join("\t");
    out.push('\n');
    for r in rows {
        out.push_str(&r.join("\t"));
        out.push('\n');
    }
    Ok(out)
}

#[pymodule]
fn _gmg(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Index>()?;
    m.add_function(wrap_pyfunction!(results_tsv_text, m)?)?;
    m.add_function(wrap_pyfunction!(hits_tsv_text, m)?)?;
    m.add("SCHEMA_VERSION", SCHEMA_VERSION)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
