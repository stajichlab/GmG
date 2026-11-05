// Web Worker for popsketch: reads each sample's first reads, sketches them, then runs the analysis. One thread, no files uploaded.
// In:  {type:'run', samples:[{id, files:[File]}], labels:string|null, compareTree:string|null, comparePcs:string|null, opts:{maxReads, scaled, flagged}}
// Out: {type:'progress', i, n, id} {type:'done', result} {type:'error', message}
import init, { Session } from "./pkg-ps/popsketch_wasm.js";

const ready = init();
const post = (m) => self.postMessage(m);

// Read a file (gunzipped when it ends in .gz) and stop after `maxReads` FASTQ records.
async function readHead(file, maxReads) {
  let stream = file.stream();
  if (/\.gz$/i.test(file.name)) stream = stream.pipeThrough(new DecompressionStream("gzip"));
  const reader = stream.getReader();
  const chunks = [];
  let lines = 0, total = 0, fastq = null;
  const maxLines = maxReads * 4;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    total += value.length;
    if (fastq === null && value.length) fastq = value[0] === 0x40;
    if (fastq) {
      for (let i = 0; i < value.length; i++) if (value[i] === 10) lines++;
      if (lines >= maxLines) { await reader.cancel(); break; }
    }
  }
  const bytes = new Uint8Array(total);
  let off = 0;
  for (const c of chunks) { bytes.set(c, off); off += c.length; }
  return bytes;
}

self.onmessage = async (e) => {
  const m = e.data;
  if (m.type !== "run") return;
  try {
    await ready;
    const { maxReads, scaled, flagged } = m.opts;
    const s = new Session(21, scaled, 0n, false);
    const ids = m.samples.map((x) => x.id);
    const halves = new Set(JSON.parse(s.halves_for(JSON.stringify(ids), 20)));
    for (let i = 0; i < m.samples.length; i++) {
      const sm = m.samples[i];
      post({ type: "progress", i, n: m.samples.length, id: sm.id });
      const r1 = await readHead(sm.files[0], maxReads);
      const r2 = sm.files[1] ? await readHead(sm.files[1], maxReads) : new Uint8Array(0);
      s.add_sample(sm.id, "reads", r1, r2, halves.has(sm.id));
    }
    post({ type: "progress", i: m.samples.length, n: m.samples.length, id: "analysing" });
    await new Promise((r) => setTimeout(r));
    const result = JSON.parse(s.analyze(m.labels || undefined, JSON.stringify({ flagged }), m.compareTree || undefined, m.comparePcs || undefined));
    s.free();
    post({ type: "done", result });
  } catch (err) {
    post({ type: "error", message: String((err && err.message) || err) });
  }
};
