// Web Worker for popsketch: reads each sample's first reads, sketches them, then runs the analysis. One thread, no files uploaded.
// In:  {type:'run', samples:[{id, files:[File]}] or bundle:{k, scaled, sketches:[{id, meta, whole, a, b}]} (finished sketches), labels:string|null, compareTree:string|null, comparePcs:string|null, opts:{maxReads, scaled, flagged, strata, permutations, collapse, collapseWithin}}
// Out: {type:'progress', i, n, id} {type:'note', message} {type:'done', result} {type:'error', message}
import init, { Session } from "./pkg-ps/popsketch_wasm.js";
import { readHead } from "./readhead.js";

const ready = init();
const post = (m) => self.postMessage(m);

self.onmessage = async (e) => {
  const m = e.data;
  if (m.type !== "run") return;
  try {
    await ready;
    const { maxReads, scaled, flagged, strata, permutations, collapse, collapseWithin } = m.opts;
    const aopts = JSON.stringify({ flagged, ...(strata ? { strata } : {}), ...(Number.isFinite(permutations) ? { permutations } : {}), ...(collapse ? { collapse } : {}), ...(collapse && collapseWithin ? { collapse_within: collapseWithin } : {}) });
    // a bundle of finished sketches sets k and the scale itself; the sketches are not made here
    const s = m.bundle ? new Session(m.bundle.k, m.bundle.scaled, 0n, false) : new Session(21, scaled, 0n, false);
    if (m.index) {
      post({ type: "progress", i: 0, n: m.samples.length, id: "loading the index" });
      s.set_index(m.index);
    }
    if (m.bundle) {
      const n = m.bundle.sketches.length;
      for (let i = 0; i < n; i++) {
        const sk = m.bundle.sketches[i];
        post({ type: "progress", i, n, id: sk.id });
        s.add_sketch(sk.meta, sk.whole, sk.a, sk.b);
      }
      post({ type: "progress", i: n, n, id: "analysing" });
      await new Promise((r) => setTimeout(r));
      const result = JSON.parse(s.analyze(m.labels || undefined, aopts, m.compareTree || undefined, m.comparePcs || undefined));
      s.free();
      post({ type: "done", result });
      return;
    }
    const ids = m.samples.map((x) => x.id);
    const halves = new Set(JSON.parse(s.halves_for(JSON.stringify(ids), 20)));
    for (let i = 0; i < m.samples.length; i++) {
      const sm = m.samples[i];
      post({ type: "progress", i, n: m.samples.length, id: sm.id });
      const note = (message) => post({ type: "note", message });
      const r1 = await readHead(sm.files[0], maxReads, { onNote: note });
      const r2 = sm.files[1] ? await readHead(sm.files[1], maxReads, { onNote: note }) : new Uint8Array(0);
      s.add_sample(sm.id, "reads", r1, r2, halves.has(sm.id));
    }
    post({ type: "progress", i: m.samples.length, n: m.samples.length, id: "analysing" });
    await new Promise((r) => setTimeout(r));
    const result = JSON.parse(s.analyze(m.labels || undefined, aopts, m.compareTree || undefined, m.comparePcs || undefined));
    s.free();
    post({ type: "done", result });
  } catch (err) {
    post({ type: "error", message: String((err && err.message) || err) });
  }
};
