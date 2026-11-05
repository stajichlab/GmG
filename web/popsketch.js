// popsketch page: choose FASTQ files (and an optional label table), group them into samples, run the worker, show the report.
const $ = (id) => document.getElementById(id);
let files = [];
let labelsText = null;
let treeText = null;
let pcsText = null;
let worker = null;

// sample id and mate: sampleA_R1.fastq.gz / sampleA_1.fq -> {id: sampleA, mate: 1}
function parse(name) {
  const base = name.replace(/\.gz$/i, "").replace(/\.(fq|fastq)$/i, "");
  const m = base.match(/^(.*?)(?:[._-]R?|R)([12])(?:_001)?$/i); // a bare trailing digit without a separator is part of the name (g0s1)
  return m && m[1] ? { id: m[1], mate: Number(m[2]) } : { id: base, mate: 1 };
}

function group(list) {
  const by = new Map();
  for (const f of list) {
    if (!/\.(fq|fastq)(\.gz)?$/i.test(f.name)) continue;
    const { id, mate } = parse(f.name);
    if (!by.has(id)) by.set(id, []);
    by.get(id)[mate === 2 ? 1 : 0] = f;
  }
  return [...by].map(([id, fs]) => ({ id, files: fs.filter(Boolean) })).sort((a, b) => a.id.localeCompare(b.id));
}

function refresh() {
  const samples = group(files);
  $("queue-status").textContent = samples.length
    ? `${samples.length} samples (${samples.reduce((n, s) => n + s.files.length, 0)} files).` + (samples.length > 200 ? " More than about 200 samples is slow in a browser; consider the command line." : "")
    : "No FASTQ files chosen.";
  $("run").disabled = samples.length < 4;
  return samples;
}

$("pick-files").addEventListener("change", (e) => { files = [...e.target.files]; refresh(); });
$("pick-dir").addEventListener("change", (e) => { files = [...e.target.files]; refresh(); });
$("pick-labels").addEventListener("change", async (e) => {
  const f = e.target.files[0];
  labelsText = f ? await f.text() : null;
  $("labels-name").textContent = f ? f.name : "";
});

$("pick-tree").addEventListener("change", async (e) => { const f = e.target.files[0]; treeText = f ? await f.text() : null; });
$("pick-pcs").addEventListener("change", async (e) => { const f = e.target.files[0]; pcsText = f ? await f.text() : null; });

function download(name, text, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function show(r) {
  $("results").hidden = false;
  const s = r.summary;
  const gap = s.clusters && s.clusters.gap ? s.clusters.gap.verdict : "n/a";
  $("summary-line").textContent = `${s.samples.total} samples, ${s.samples.in_analysis} in the analysis; cut verdict: ${gap}.`;
  $("report").srcdoc = r.html;
  const row = $("downloads");
  row.replaceChildren();
  const add = (name, text, type = "text/plain") => {
    const b = document.createElement("button");
    b.className = "btn"; b.type = "button"; b.textContent = name;
    b.addEventListener("click", () => download(name, text, type));
    row.append(b);
  };
  add("report.html", r.html, "text/html");
  add("popsketch.json", JSON.stringify(r.summary, null, 1), "application/json");
  add("qc.tsv", ["id\tstatus\tflags", ...r.qc.map((q) => `${q.id}\t${q.status}\t${q.flags.map((f) => f.code).join(",")}`)].join("\n") + "\n");
  if (r.clusters) add("clusters.tsv", ["id\tcluster", ...r.clusters.map(([id, k]) => `${id}\t${k}`)].join("\n") + "\n");
  if (r.pca) add("pca.tsv", ["id\t" + r.pca[0].pcs.map((_, i) => `PC${i + 1}`).join("\t"), ...r.pca.map((p) => `${p.id}\t${p.pcs.join("\t")}`)].join("\n") + "\n");
  if (r.newick) add("tree.nwk", r.newick + "\n");
}

$("run").addEventListener("click", () => {
  const samples = refresh();
  if (worker) worker.terminate();
  worker = new Worker("popsketch-worker.js", { type: "module" });
  $("run").disabled = true;
  $("run-bar").hidden = false;
  worker.onmessage = (e) => {
    const m = e.data;
    if (m.type === "progress") {
      $("run-status").textContent = m.id === "analysing" ? "Analysing…" : `Sketching ${m.i + 1} of ${m.n}: ${m.id}`;
      $("run-bar").firstElementChild.style.width = `${(100 * m.i) / m.n}%`;
    } else if (m.type === "done") {
      $("run-status").textContent = "Done.";
      $("run-bar").hidden = true;
      $("run").disabled = false;
      show(m.result);
    } else if (m.type === "error") {
      $("run-status").textContent = `Error: ${m.message}`;
      $("run-bar").hidden = true;
      $("run").disabled = false;
    }
  };
  worker.postMessage({
    type: "run", samples, labels: labelsText, compareTree: treeText, comparePcs: pcsText,
    opts: { maxReads: Number($("max-reads").value), scaled: Number($("scaled").value), flagged: $("flagged").value },
  });
});
refresh();
