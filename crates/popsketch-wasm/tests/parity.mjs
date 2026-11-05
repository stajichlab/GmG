// Parity: the same simulated reads through the native `popsketch run` and through the wasm Session must give the same
// audit, noise floor, clusters, PCA and label check. Run with `pixi run popsketch-wasm-test`
// (needs target/release/popsketch and crates/popsketch-wasm/pkg-node).
import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import assert from "node:assert/strict";

const require = createRequire(import.meta.url);
const here = path.dirname(new URL(import.meta.url).pathname);
const { Session } = require(path.join(here, "../pkg-node/popsketch_wasm.js"));
const bin = path.join(here, "../../../target/release/popsketch");

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "popsketch-parity-"));
execFileSync(bin, ["demo", "-o", tmp], { stdio: ["ignore", "ignore", "inherit"] });
const out = path.join(tmp, "out");
execFileSync(
  bin,
  ["run", path.join(tmp, "reads"), "--labels", path.join(tmp, "labels.tsv"), "-o", out, "--scaled", "5", "--max-bases", "0", "--floor-samples", "20", "--threads", "2", "--quiet"],
  { stdio: ["ignore", "ignore", "inherit"] },
);
const native = JSON.parse(fs.readFileSync(path.join(out, "popsketch.json"), "utf8"));
const nativeClusters = new Map(fs.readFileSync(path.join(out, "clusters.tsv"), "utf8").trim().split("\n").slice(1).map((l) => l.split("\t")).map((c) => [c[0], c[1]]));

const files = fs.readdirSync(path.join(tmp, "reads")).filter((f) => f.endsWith(".fastq")).sort();
const ids = files.map((f) => f.replace(/\.fastq$/, ""));
const s = new Session(21, 5, 0n, false);
const withHalves = JSON.parse(s.halves_for(JSON.stringify(ids), 20));
files.forEach((f, i) => s.add_sample(ids[i], "reads", new Uint8Array(fs.readFileSync(path.join(tmp, "reads", f))), new Uint8Array(0), withHalves.includes(ids[i])));
const res = JSON.parse(s.analyze(fs.readFileSync(path.join(tmp, "labels.tsv"), "utf8"), "{}"));
const w = res.summary;
// numbers agree to rounding (the native code sums in parallel, the wasm code on one thread); everything else is exact
const same = (a, b, what) => {
  const walk = (x, y, at) => {
    if (typeof x === "number" && typeof y === "number") {
      assert.ok(Math.abs(x - y) <= 1e-9 * Math.max(1, Math.abs(x), Math.abs(y)), `${what}${at}: ${x} vs ${y}`);
    } else if (x && y && typeof x === "object") {
      assert.deepEqual(Object.keys(x).sort(), Object.keys(y).sort(), `${what}${at}: keys`);
      for (const k of Object.keys(x)) walk(x[k], y[k], `${at}.${k}`);
    } else assert.deepEqual(x, y, `${what}${at}`);
  };
  walk(a, b, "");
};

same(w.samples, native.samples, "sample counts by status");
same(w.flags, native.flags, "flags");
assert.equal(w.noise_floor.source, native.noise_floor.source);
assert.ok(Math.abs(w.noise_floor.value - native.noise_floor.value) < 1e-12, `noise floor ${w.noise_floor.value} vs ${native.noise_floor.value}`);
same(w.clusters, native.clusters, "cluster verdict");
same(w.labels, native.labels, "label check");
same(w.ordination, native.ordination, "ordination");
assert.equal(w.tree.branches_collapsed, native.tree.branches_collapsed);
// cluster membership: same partition (labels may be numbered differently, so compare as a partition)
const wasm = new Map(res.clusters.map(([id, k]) => [id, String(k)]));
assert.equal(wasm.size, nativeClusters.size);
const pairs = (m) => { const o = [...m.keys()].sort(); const r = new Set(); for (const a of o) for (const b of o) if (a < b && m.get(a) === m.get(b)) r.add(a + "|" + b); return r; };
assert.deepEqual([...pairs(wasm)].sort(), [...pairs(nativeClusters)].sort(), "cluster partition");
assert.ok(res.html.includes("1. Audit") && res.newick && res.newick.endsWith(";"), "report and tree present");
assert.equal(res.pca.length, wasm.size);
console.log(`popsketch wasm parity OK: ${ids.length} samples, ${wasm.size} analysed, verdict ${JSON.stringify(w.clusters.gap.verdict)}`);
fs.rmSync(tmp, { recursive: true, force: true });
