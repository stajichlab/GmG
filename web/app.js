// GmG web app: UI only. The worker (worker.js) does the sketching and classification with the same Rust code as the CLI.
const $ = (id) => document.getElementById(id);
const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });

const SEQ_EXT = /\.(fa|fna|fasta|fas|ffn|fq|fastq)(\.gz)?$/i;
const RANKS = ["phylum", "class", "order", "family", "genus", "species"];
const STATUS_LABEL = { ok: "classified", low_confidence: "low confidence", no_match: "no match", error: "error" };

let queue = [];           // File objects waiting to be classified
let records = [];         // result records so far
let outputs = null;       // files produced at the end of a run
let sortKey = "id", sortDir = 1;
let indexReady = false, running = false;

// ---------------------------------------------------------------- reference index
async function loadManifest() {
  try {
    const r = await fetch("index/manifest.json", { cache: "no-cache" });
    if (!r.ok) return;
    const list = await r.json();
    const sel = $("index-select");
    sel.innerHTML = "";
    for (const e of list) {
      const o = document.createElement("option");
      o.value = e.file; o.textContent = `${e.title} (${e.size_mb} MB)`; sel.append(o);
    }
  } catch { /* no bundled index: the file picker still works */ }
}
function showProgress(bar, frac) { bar.hidden = false; bar.firstElementChild.style.width = `${Math.round(frac * 100)}%`; }
function setStatus(el, text, err) { el.textContent = text; el.classList.toggle("err", !!err); }

$("index-load").onclick = () => {
  const file = $("index-select").value;
  if (!file) return setStatus($("index-status"), "No bundled index is available here: use “Use a .gmgi file…”.", true);
  setStatus($("index-status"), "Loading the index…");
  worker.postMessage({ type: "index", url: `index/${file}` });
};
$("index-file").onchange = (e) => {
  const f = e.target.files[0];
  if (!f) return;
  setStatus($("index-status"), `Reading ${f.name}…`);
  worker.postMessage({ type: "index", file: f });
};

function showIndex(info) {
  indexReady = true;
  const m = info.manifest || {};
  const facts = [["Genomes", info.n_genomes.toLocaleString()], ["k-mer size / scaled", `${info.k} / ${info.scaled}`]];
  if (m.taxonomy_version) facts.push(["Taxonomy", m.taxonomy_version]);
  if (m.built_utc) facts.push(["Built", m.built_utc.slice(0, 10)]);
  if (m.note) facts.push(["Note", m.note]);
  const dl = $("index-facts");
  dl.replaceChildren(...facts.flatMap(([k, v]) => { const a = document.createElement("dt"), b = document.createElement("dd"); a.textContent = k; b.textContent = v; return [a, b]; }));
  dl.hidden = false;
  $("index-bar").hidden = true;
  setStatus($("index-status"), "Index ready. Files are compared with it on this computer.");
  updateRun();
}

// ---------------------------------------------------------------- choosing files
const meta = new WeakMap(); // File -> { label, sha } for the list under the drop zone
const HASH_MAX = 256e6;      // crypto.subtle needs the whole file in memory, so larger files are listed without a hash
const fmtSize = (n) => (n >= 1e9 ? `${(n / 1e9).toFixed(1)} GB` : n >= 1e6 ? `${(n / 1e6).toFixed(1)} MB` : `${Math.max(1, Math.round(n / 1e3))} kB`);
let hashing = Promise.resolve();
function hashFile(f) {
  const m = meta.get(f);
  if (f.size > HASH_MAX) { m.sha = "too large to hash here"; return; }
  hashing = hashing.then(async () => {
    try {
      const d = await crypto.subtle.digest("SHA-256", await f.arrayBuffer());
      m.sha = [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
    } catch { m.sha = "hash unavailable"; }
    renderQueue();
  });
}
function renderQueue() {
  const ul = $("queue-list");
  ul.hidden = queue.length === 0;
  ul.replaceChildren(...queue.map((f) => {
    const m = meta.get(f), li = document.createElement("li");
    const name = document.createElement("span"); name.className = "qname"; name.textContent = f.name;
    const info = document.createElement("span"); info.className = "qinfo";
    const sha = m.sha ? (m.sha.length === 64 ? `sha256 ${m.sha.slice(0, 12)}…` : m.sha) : "hashing…";
    info.textContent = [m.label, fmtSize(f.size), sha].filter(Boolean).join(" · ");
    if (m.sha && m.sha.length === 64) info.title = `sha256 ${m.sha}`;
    const x = document.createElement("button"); x.type = "button"; x.className = "qx"; x.textContent = "✕"; x.title = `Remove ${f.name}`; x.setAttribute("aria-label", `Remove ${f.name} from the list`);
    x.onclick = () => { queue = queue.filter((q) => q !== f); queueStatus(); };
    li.append(name, info, x);
    return li;
  }));
}
function queueStatus(extra) {
  setStatus($("queue-status"), queue.length ? `${queue.length} file${queue.length === 1 ? "" : "s"} ready${extra || ""}. Remove any with ✕.` : `No bite yet: no files chosen${extra || ""}.`);
  renderQueue(); updateRun();
}
function addFiles(files, label) {
  const keep = [...files].filter((f) => SEQ_EXT.test(f.name));
  const same = (a, b) => a.name === b.name && a.size === b.size && a.lastModified === b.lastModified;
  const fresh = keep.filter((f, i) => !queue.some((q) => same(q, f)) && keep.findIndex((g) => same(g, f)) === i);
  const dup = keep.length - fresh.length, skipped = files.length - keep.length;
  for (const f of fresh) { meta.set(f, { label: label || "", sha: "" }); hashFile(f); }
  queue = queue.concat(fresh);
  const notes = [skipped ? `${skipped} skipped: not FASTA/FASTQ` : "", dup ? `${dup} already in the list` : ""].filter(Boolean).join("; ");
  queueStatus(notes ? ` (${notes})` : "");
}
function updateRun() { $("run").disabled = running || !indexReady || queue.length === 0; }

$("pick-files").onchange = (e) => { addFiles(e.target.files); e.target.value = ""; };
$("pick-dir").onchange = (e) => { addFiles(e.target.files); e.target.value = ""; };
$("drop").onkeydown = (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); $("pick-files").click(); } };
for (const ev of ["dragenter", "dragover"]) $("drop").addEventListener(ev, (e) => { e.preventDefault(); $("drop").classList.add("over"); });
for (const ev of ["dragleave", "drop"]) $("drop").addEventListener(ev, (e) => { e.preventDefault(); $("drop").classList.remove("over"); });
$("drop").addEventListener("drop", async (e) => {
  const items = [...(e.dataTransfer.items || [])].map((i) => i.webkitGetAsEntry && i.webkitGetAsEntry()).filter(Boolean);
  if (!items.length) return addFiles(e.dataTransfer.files);
  const files = [];
  const walk = async (entry) => {
    if (entry.isFile) files.push(await new Promise((res, rej) => entry.file(res, rej)));
    else if (entry.isDirectory) {
      const reader = entry.createReader();
      for (;;) { // readEntries returns batches until it returns an empty one
        const batch = await new Promise((res, rej) => reader.readEntries(res, rej));
        if (!batch.length) break;
        for (const c of batch) await walk(c);
      }
    }
  };
  for (const it of items) await walk(it);
  addFiles(files);
});

// ---------------------------------------------------------------- example genomes
// examples/examples.tsv (shared with scripts/fetch_examples.sh): id, accession, name, expected, note, url. NCBI's FTP server sends
// CORS headers, so the browser downloads the gzipped FASTA directly; the worker already reads .gz.
async function loadExamples() {
  try {
    const r = await fetch("examples/examples.tsv", { cache: "no-cache" });
    if (!r.ok) return;
    const [head, ...lines] = (await r.text()).trim().split("\n").map((l) => l.split("\t"));
    const rows = lines.map((c) => Object.fromEntries(head.map((h, i) => [h, c[i]])));
    $("example-buttons").replaceChildren(...rows.map((ex) => {
      const b = document.createElement("button");
      b.type = "button"; b.className = "btn"; b.textContent = ex.name; b.title = `${ex.accession}: ${ex.note}`;
      b.onclick = () => addExample(ex, b);
      return b;
    }));
    $("examples").hidden = false;
  } catch { /* examples are optional */ }
}
// Any NCBI assembly: the Datasets API (CORS-enabled) gives the assembly name, which with the accession fixes the FTP path.
const ACC = /^GC[AF]_\d{9}(\.\d+)?$/i;
async function addAccession(raw) {
  const acc = raw.trim().toUpperCase(), btn = $("acc-go"), box = $("acc");
  if (!ACC.test(acc)) return setStatus($("queue-status"), `"${raw.trim()}" is not an NCBI assembly accession. Expected GCF_ or GCA_, nine digits, and optionally a version: GCF_000146045.2.`, true);
  btn.disabled = box.disabled = true;
  setStatus($("queue-status"), `Looking up ${acc} at NCBI…`);
  try {
    const r = await fetch(`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/dataset_report`);
    if (!r.ok) throw new Error(`NCBI answered ${r.status}`);
    const rep = ((await r.json()).reports || [])[0];
    if (!rep) throw new Error("no such assembly");
    const id = rep.accession, d = id.slice(4, 13).match(/.{3}/g).join("/");
    const dir = `${id}_${rep.assembly_info.assembly_name.replace(/[^A-Za-z0-9._-]/g, "_")}`;
    const url = `https://ftp.ncbi.nlm.nih.gov/genomes/all/${id.slice(0, 3)}/${d}/${dir}/${dir}_genomic.fna.gz`;
    const mb = Math.round((rep.assembly_stats && rep.assembly_stats.total_sequence_length || 0) / 1e6);
    setStatus($("queue-status"), `Downloading ${rep.organism.organism_name} (${id}${mb ? `, ${mb} Mb` : ""}) from NCBI…`);
    const g = await fetch(url);
    if (!g.ok) throw new Error(`NCBI answered ${g.status} for the genome file`);
    addFiles([new File([await g.blob()], `${id}.fna.gz`)], rep.organism.organism_name);
    box.value = "";
  } catch (err) {
    setStatus($("queue-status"), `Could not fetch ${acc}: ${err.message}. Check the accession and your connection.`, true);
  } finally { btn.disabled = box.disabled = false; }
}
$("acc-go").onclick = () => addAccession($("acc").value);
const useEg = () => { $("acc").value = $("acc-eg").textContent; $("acc").focus(); };
$("acc-eg").onclick = useEg;
$("acc-eg").onkeydown = (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); useEg(); } };
$("acc").onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); addAccession($("acc").value); } };

async function addExample(ex, btn) {
  btn.disabled = true;
  setStatus($("queue-status"), `Downloading ${ex.name} (${ex.accession}) from NCBI…`);
  try {
    const r = await fetch(ex.url);
    if (!r.ok) throw new Error(`NCBI answered ${r.status}`);
    addFiles([new File([await r.blob()], `${ex.accession}.fna.gz`)], ex.name);
  } catch (err) {
    setStatus($("queue-status"), `Could not download ${ex.name}: ${err.message}. Check your connection, or use scripts/fetch_examples.sh.`, true);
  } finally { btn.disabled = false; }
}

// ---------------------------------------------------------------- running
$("th").onchange = () => worker.postMessage({ type: "thresholds", name: $("th").value });
$("run").onclick = () => {
  records = []; outputs = null;
  $("results").hidden = false;
  render();
  running = true; updateRun();
  $("cancel").hidden = false;
  const opts = { maxReads: Number($("max-reads").value) || 0, minAbund: Number($("min-abund").value) || 0, top: 5 };
  worker.postMessage({ type: "classify", files: queue, opts });
  queue = [];
  queueStatus();
};
$("cancel").onclick = () => worker.postMessage({ type: "cancel" });

worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === "index-progress") {
    if (m.cached) setStatus($("index-status"), "Loading the index from this browser's cache…");
    else if (m.total) { showProgress($("index-bar"), m.loaded / m.total); setStatus($("index-status"), `Downloading the index… ${(m.loaded / 1e6).toFixed(0)} of ${(m.total / 1e6).toFixed(0)} MB`); }
  } else if (m.type === "index-ready") {
    showIndex(m.info);
  } else if (m.type === "progress") {
    showProgress($("run-bar"), m.i / m.n);
    setStatus($("run-status"), `Classifying ${m.i + 1} of ${m.n}: ${m.name}`);
  } else if (m.type === "record") {
    records.push(m.record);
    render();
  } else if (m.type === "done") {
    outputs = m;
    running = false; updateRun();
    $("cancel").hidden = true; $("run-bar").hidden = true;
    setStatus($("run-status"), m.cancelled ? `Cancelled after ${m.n} files.` : `Done: ${m.n} file${m.n === 1 ? "" : "s"}.`);
    if (m.report_html) $("report").srcdoc = m.report_html;
    render();
  } else if (m.type === "error") {
    running = false; updateRun();
    $("cancel").hidden = true; $("run-bar").hidden = true;
    setStatus($("index-bar").hidden && !indexReady ? $("index-status") : $("run-status"), m.message, true);
  }
};

// ---------------------------------------------------------------- results table
const call = (r) => {
  const c = r.guess && r.guess.calls;
  if (!c) return "";
  for (let i = 5; i >= 0; i--) if (c[RANKS[i]].taxon) return c[RANKS[i]].taxon;
  const t = r.guess.tentative_genus;
  return t ? `${t.taxon} (tentative genus)` : "";
};
const phylum = (r) => (r.guess && r.guess.calls && r.guess.calls.phylum.taxon) || "";
const ani = (r) => (r.guess && r.guess.top ? r.guess.top.ani : null);
const why = (r) => {
  if (r.status === "error") return `${r.error_code}: ${r.message} ${r.error_action || ""}`.trim();
  const w = r.guess && r.guess.withheld;
  return w ? `${w.rank} withheld: ${w.detail}` : "";
};
const COLS = [
  ["id", "Input", (r) => r.id],
  ["status", "Status", (r) => r.status],
  ["rank", "Deepest rank", (r) => (r.best_rank || "")],
  ["call", "Call", call],
  ["phylum", "Phylum", phylum],
  ["ani", "ANI", ani],
];

function render() {
  const counts = {};
  for (const r of records) counts[r.status] = (counts[r.status] || 0) + 1;
  $("chips").replaceChildren(...Object.entries(counts).map(([k, v]) => { const s = document.createElement("span"); s.className = "chip"; const b = document.createElement("b"); b.textContent = v; s.append(b, STATUS_LABEL[k] || k); return s; }));
  const q = $("filter").value.trim().toLowerCase(), sf = $("status-filter").value;
  const get = Object.fromEntries(COLS.map(([k, , f]) => [k, f]));
  let rows = records.filter((r) => (!sf || r.status === sf) && (!q || `${r.id} ${call(r)} ${phylum(r)}`.toLowerCase().includes(q)));
  rows = rows.sort((a, b) => {
    const x = get[sortKey](a), y = get[sortKey](b);
    if (x === y) return 0;
    if (x === null || x === "") return 1;
    if (y === null || y === "") return -1;
    return (x < y ? -1 : 1) * sortDir;
  });
  const thead = $("table").tHead;
  thead.replaceChildren(Object.assign(document.createElement("tr"), {}));
  for (const [k, label] of COLS) {
    const th = document.createElement("th");
    th.textContent = label; th.tabIndex = 0; th.scope = "col";
    if (k === sortKey) th.setAttribute("aria-sort", sortDir > 0 ? "ascending" : "descending");
    const go = () => { if (sortKey === k) sortDir = -sortDir; else { sortKey = k; sortDir = 1; } render(); };
    th.onclick = go; th.onkeydown = (e) => { if (e.key === "Enter") go(); };
    thead.rows[0].append(th);
  }
  const th = document.createElement("th"); th.textContent = "Flags and why"; thead.rows[0].append(th);
  const body = $("table").tBodies[0];
  body.replaceChildren(...rows.map((r) => {
    const tr = document.createElement("tr");
    for (const [k, , f] of COLS) {
      const td = document.createElement("td");
      const v = f(r);
      if (k === "status") { const s = document.createElement("span"); s.className = `st ${r.status}`; s.textContent = STATUS_LABEL[r.status] || r.status; td.append(s); }
      else if (k === "ani") {
        td.className = "num"; td.textContent = v === null ? "" : v.toFixed(4);
        if (v !== null) { const s = document.createElement("span"); s.className = "sound"; s.style.setProperty("--v", Math.min(1, Math.max(0, (v - 0.7) / 0.3)).toFixed(3)); td.append(s); }
      }
      else td.textContent = v || "";
      if (k === "call" && r.guess && r.guess.calls) {
        // depth profile: one tick per rank with a call, phylum first
        const p = document.createElement("span"); p.className = "prof"; p.title = "ranks called, phylum to species";
        for (const rk of RANKS) { const i = document.createElement("i"); if (r.guess.calls[rk].taxon) i.className = "on"; p.append(i); }
        td.append(p);
      }
      tr.append(td);
    }
    const td = document.createElement("td");
    for (const fl of r.qc_flags || []) { const s = document.createElement("span"); s.className = "flag"; s.textContent = fl; td.append(s); }
    const w = why(r);
    if (w) { const s = document.createElement("div"); s.className = "why"; s.textContent = w; td.append(s); }
    tr.append(td);
    return tr;
  }));
}
$("filter").oninput = render;
$("status-filter").onchange = render;

// ---------------------------------------------------------------- downloads
function save(name, text, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name; document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}
for (const b of document.querySelectorAll("[data-dl]")) b.onclick = () => {
  if (!outputs) return;
  const name = b.dataset.dl;
  if (name === "results.tsv") save(name, outputs.results_tsv, "text/tab-separated-values");
  else if (name === "hits.tsv") save(name, outputs.hits_tsv, "text/tab-separated-values");
  else save(name, outputs.report_html, "text/html");
};
$("print").onclick = () => {
  if (!outputs || !outputs.report_html) return;
  // the report is shown in a sandboxed frame; printing happens from its own tab so the browser's "save as PDF" gets the whole page
  const url = URL.createObjectURL(new Blob([outputs.report_html], { type: "text/html" }));
  const w = window.open(url, "_blank");
  if (!w) return setStatus($("run-status"), "Allow pop-ups to print, or download report.html and print it from there.", true);
  w.addEventListener("load", () => w.print());
  setTimeout(() => URL.revokeObjectURL(url), 60000);
};

loadManifest();
loadExamples();
