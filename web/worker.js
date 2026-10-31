// Web Worker: loads the index, classifies files and renders the outputs. All heavy work happens here, off the UI thread.
// Messages in:  {type:'index', url} | {type:'index', file} | {type:'thresholds', name}
//               {type:'classify', files:[File], opts:{maxReads, minAbund, top}} | {type:'cancel'}
// Messages out: {type:'index-progress', loaded, total} {type:'index-ready', info} {type:'progress', i, n, name}
//               {type:'record', record} {type:'done', ...outputs} {type:'error', message}
import init, { Index, make_results_tsv, make_hits_tsv, render_report } from "./pkg/gmg_wasm.js";

const CACHE = "gmg-index-v1";
let ix = null;
let thresholds = "known";
let cancelled = false;
const ready = init();

const post = (m) => self.postMessage(m);
const stem = (name) => name.replace(/\.gz$/i, "").replace(/\.(fa|fna|fasta|fas|ffn|fq|fastq)$/i, "");

async function fetchIndex(url) {
  const cache = await caches.open(CACHE);
  let resp = await cache.match(url);
  if (resp) {
    post({ type: "index-progress", loaded: 1, total: 1, cached: true });
    return new Uint8Array(await resp.arrayBuffer());
  }
  resp = await fetch(url);
  if (!resp.ok) throw new Error(`could not fetch the index (${resp.status})`);
  const total = Number(resp.headers.get("content-length")) || 0;
  const copy = resp.clone(); // stored in the cache while the original is read for progress
  const reader = resp.body.getReader();
  const chunks = [];
  let loaded = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    loaded += value.length;
    post({ type: "index-progress", loaded, total });
  }
  cache.put(url, copy).catch(() => {}); // offline use next time; a failure (quota, private mode) only costs the cache
  const bytes = new Uint8Array(loaded);
  let off = 0;
  for (const c of chunks) { bytes.set(c, off); off += c.length; }
  return bytes;
}

function setIndex(bytes) {
  if (ix) ix.free();
  ix = new Index(bytes);
  ix.use_thresholds(thresholds);
  let manifest = null;
  try { manifest = JSON.parse(ix.manifest() || "null"); } catch { /* none */ }
  post({ type: "index-ready", info: { n_genomes: ix.n_genomes(), k: ix.ksize(), scaled: ix.scaled(), manifest } });
}

// Read a whole file, or for FASTQ stop after `maxReads` records (so a 2 GB read file is not read to the end).
async function readFile(file, maxReads) {
  const gz = /\.gz$/i.test(file.name);
  let stream = file.stream();
  if (gz) stream = stream.pipeThrough(new DecompressionStream("gzip"));
  const reader = stream.getReader();
  const chunks = [];
  let lines = 0, total = 0, fastq = null;
  const maxLines = maxReads > 0 ? maxReads * 4 + 4 : Infinity;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    total += value.length;
    if (fastq === null && value.length) fastq = value[0] === 0x40; // '@'
    if (fastq && Number.isFinite(maxLines)) {
      for (let i = 0; i < value.length; i++) if (value[i] === 10) lines++;
      if (lines >= maxLines) { await reader.cancel(); break; }
    }
  }
  const bytes = new Uint8Array(total);
  let off = 0;
  for (const c of chunks) { bytes.set(c, off); off += c.length; }
  return bytes;
}

async function classifyAll(files, opts) {
  cancelled = false;
  const t0 = performance.now();
  const started = new Date().toISOString().replace(/\.\d+Z$/, "Z");
  ix.clear_session();
  const records = [];
  const used = new Set();
  for (let i = 0; i < files.length; i++) {
    if (cancelled) break;
    const f = files[i];
    post({ type: "progress", i, n: files.length, name: f.name });
    let id = stem(f.name);
    for (let k = 2; used.has(id); k++) id = `${stem(f.name)}_${k}`;
    used.add(id);
    let rec;
    try {
      const bytes = await readFile(f, opts.maxReads);
      rec = JSON.parse(ix.classify(id, f.name, bytes, "auto", opts.minAbund || 0, BigInt(opts.maxReads || 0), opts.top || 5));
    } catch (e) {
      rec = JSON.parse(JSON.stringify({ schema_version: 8, id, path: f.name, size_bytes: f.size, kind: "assembly", status: "error", best_rank: null,
        error_code: "unreadable", message: "The file could not be read.", error_action: String(e && e.message || e), qc: null, qc_flags: [], guess: null }));
    }
    records.push(rec);
    post({ type: "record", record: rec });
    await new Promise((r) => setTimeout(r, 0)); // let a cancel message through
  }
  const jsonl = records.map((r) => JSON.stringify(r)).join("\n") + "\n";
  const rel = ix.relationships(jsonl);
  const run = ix.run_json(started, (performance.now() - t0) / 1000);
  post({
    type: "done", cancelled, n: records.length,
    results_tsv: make_results_tsv(jsonl), hits_tsv: make_hits_tsv(jsonl), results_jsonl: jsonl, run_json: run,
    report_html: records.length ? render_report(jsonl, run, rel) : "",
  });
}

self.onmessage = async (e) => {
  const m = e.data;
  try {
    await ready;
    if (m.type === "index") {
      const bytes = m.file ? new Uint8Array(await m.file.arrayBuffer()) : await fetchIndex(m.url);
      setIndex(bytes);
    } else if (m.type === "thresholds") {
      thresholds = m.name;
      if (ix) ix.use_thresholds(m.name);
    } else if (m.type === "classify") {
      if (!ix) throw new Error("load a reference index first");
      await classifyAll(m.files, m.opts || {});
    } else if (m.type === "cancel") {
      cancelled = true;
    }
  } catch (err) {
    post({ type: "error", message: String((err && err.message) || err) });
  }
};
