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

const out3 = path.join(tmp, "out3");
execFileSync(bin, ["run", path.join(tmp, "reads"), "--labels", path.join(tmp, "labels.tsv"), "-o", out3, "--scaled", "5", "--max-bases", "0", "--floor-samples", "20", "--threads", "2", "--quiet", "--report-format", "md"], { stdio: ["ignore", "ignore", "inherit"] });
const files = fs.readdirSync(path.join(tmp, "reads")).filter((f) => f.endsWith(".fastq")).sort();
const ids = files.map((f) => f.replace(/\.fastq$/, ""));
const s = new Session(21, 5, 0n, false);
const withHalves = JSON.parse(s.halves_for(JSON.stringify(ids), 20));
files.forEach((f, i) => s.add_sample(ids[i], "reads", new Uint8Array(fs.readFileSync(path.join(tmp, "reads", f))), new Uint8Array(0), withHalves.includes(ids[i])));
const res = JSON.parse(s.analyze(fs.readFileSync(path.join(tmp, "labels.tsv"), "utf8"), "{}", undefined, undefined));
const res3 = res;
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
// comparison with an external tree and PCs: the first run's own tree and PCs stand in for SNP results
const tree = fs.readFileSync(path.join(out, "tree.nwk"), "utf8");
// the PC table the comparison reads is id plus PC columns only
const pcs = fs.readFileSync(path.join(out, "pca.tsv"), "utf8").trim().split("\n").map((l) => l.split("\t")).map((c) => [c[0], ...c.slice(2, -1)].join("\t")).join("\n") + "\n";
fs.writeFileSync(path.join(tmp, "snp_pcs.tsv"), pcs);
const out2 = path.join(tmp, "out2");
execFileSync(bin, ["run", path.join(tmp, "reads"), "--labels", path.join(tmp, "labels.tsv"), "-o", out2, "--scaled", "5", "--max-bases", "0", "--floor-samples", "20", "--threads", "2", "--quiet", "--compare-tree", path.join(out, "tree.nwk"), "--compare-pcs", path.join(tmp, "snp_pcs.tsv")], { stdio: ["ignore", "ignore", "inherit"] });
const native2 = JSON.parse(fs.readFileSync(path.join(out2, "popsketch.json"), "utf8"));
const s2 = new Session(21, 5, 0n, false);
files.forEach((f, i) => s2.add_sample(ids[i], "reads", new Uint8Array(fs.readFileSync(path.join(tmp, "reads", f))), new Uint8Array(0), withHalves.includes(ids[i])));
const w2 = JSON.parse(s2.analyze(fs.readFileSync(path.join(tmp, "labels.tsv"), "utf8"), "{}", tree, pcs)).summary;
assert.ok(native2.compare.tree && native2.compare.pcs, "native comparison present");
same(w2.compare, native2.compare, "comparison with an external tree and PCs");
// Markdown report: identical text and figures
const nativeMd = fs.readFileSync(path.join(out3, "report.md"), "utf8");
assert.equal(res3.markdown.md, nativeMd, "report.md");
for (const [file, svg] of res3.markdown.figures) assert.equal(svg, fs.readFileSync(path.join(out3, "figures", file), "utf8"), `figure ${file}`);

// index screen: a small index of three of the sets (as genomes), through `gmg build`, loaded in both
const gmg = path.join(here, "../../../target/release/gmg");
const ixdir = path.join(tmp, "ix");
fs.mkdirSync(ixdir);
const meta = ["id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tphylum_taxid\tclass_taxid\torder_taxid\tfamily_taxid\tgenus_taxid\tspecies_taxid\tlabel"];
const refs = [["g0s1", "Alpha"], ["g1s1", "Beta"], ["g2s1", "Gamma"]];
refs.forEach(([id, genus], i) => {
  fs.copyFileSync(path.join(tmp, "reads", id + ".fastq"), path.join(ixdir, id + ".fastq"));
  meta.push([id, 100 + i, "Ascomycota", "C", "O", "F", genus, genus + " sp", 4890, 0, 0, 0, 200 + i, 100 + i, genus].join("\t"));
});
fs.writeFileSync(path.join(ixdir, "meta.tsv"), meta.join("\n") + "\n");
execFileSync(gmg, ["sketch", "-k", "21", "-s", "10", "-o", path.join(ixdir, "sk"), ...refs.map(([id]) => path.join(ixdir, id + ".fastq"))], { stdio: ["ignore", "ignore", "inherit"] });
execFileSync(gmg, ["build", "--meta", path.join(ixdir, "meta.tsv"), "--sketches", path.join(ixdir, "sk"), "-o", path.join(ixdir, "raw.gmgi")], { stdio: ["ignore", "ignore", "inherit"] });
execFileSync(gmg, ["convert", path.join(ixdir, "raw.gmgi"), "-o", path.join(ixdir, "ix.gmgi")], { stdio: ["ignore", "ignore", "inherit"] });
const out4 = path.join(tmp, "out4");
execFileSync(bin, ["run", path.join(tmp, "reads"), "-o", out4, "--scaled", "5", "--max-bases", "0", "--floor-samples", "20", "--threads", "2", "--quiet", "--index", path.join(ixdir, "ix.gmgi")], { stdio: ["ignore", "ignore", "inherit"] });
const native4 = JSON.parse(fs.readFileSync(path.join(out4, "popsketch.json"), "utf8"));
const s4 = new Session(21, 5, 0n, false);
s4.set_index(new Uint8Array(fs.readFileSync(path.join(ixdir, "ix.gmgi"))));
files.forEach((f, i) => s4.add_sample(ids[i], "reads", new Uint8Array(fs.readFileSync(path.join(tmp, "reads", f))), new Uint8Array(0), withHalves.includes(ids[i])));
const w4 = JSON.parse(s4.analyze(undefined, "{}", undefined, undefined));
same(w4.summary.flags, native4.flags, "flags with an index");
same(w4.summary.samples, native4.samples, "samples with an index");
assert.ok(native4.options.index, "native used the index");
assert.ok(!w4.warnings.some((x) => x.includes("index screen skipped")), "wasm did not skip the index screen");
const nativeIdx = fs.readFileSync(path.join(out4, "qc.tsv"), "utf8").trim().split("\n").slice(1).map((l) => l.split("\t")).map((c) => [c[0], c[16]]);
const wasmIdx = new Map(w4.qc.map((q) => [q.id, q]));
assert.equal(wasmIdx.size, nativeIdx.length, "qc rows with an index");
console.log(`index screen: ${nativeIdx.filter((r) => r[1] && r[1] !== "none").length} samples with an organism call in the native run`);
console.log(`popsketch wasm parity OK: ${ids.length} samples, ${wasm.size} analysed, verdict ${JSON.stringify(w.clusters.gap.verdict)}`);
fs.rmSync(tmp, { recursive: true, force: true });
