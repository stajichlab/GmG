// Parity: 20+ inputs classified in wasm must give the same records and TSVs as `gmg run` (issue #24 acceptance).
// Run with `pixi run wasm-test` (needs target/release/gmg and crates/gmg-wasm/pkg-node).
import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import assert from "node:assert/strict";

const require = createRequire(import.meta.url);
const here = path.dirname(new URL(import.meta.url).pathname);
const { Index, make_results_tsv, make_hits_tsv, render_report } = require(path.join(here, "../pkg-node/gmg_wasm.js"));
const gmg = path.join(here, "../../../target/release/gmg");

let seed = 777;
const rnd = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 2 ** 32);
const randGenome = (n) => Array.from({ length: n }, () => "ACGT"[Math.floor(rnd() * 4)]).join("");
const mutate = (s, rate) => [...s].map((c) => (rnd() < rate ? "ACGT"[Math.floor(rnd() * 4)] : c)).join("");
const fasta = (name, seq) => `>${name}\n${seq.match(/.{1,70}/g).join("\n")}\n`;
const fastq = (g, n, len) => Array.from({ length: n }, (_, i) => {
  const a = Math.floor(rnd() * (g.length - len));
  return `@r${i}\n${g.slice(a, a + len)}\n+\n${"I".repeat(len)}\n`;
}).join("");

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "gmg-parity-"));
const NG = 12;
const genomes = Array.from({ length: NG }, () => randGenome(60_000));
const genera = ["Aspergillus", "Fusarium", "Agaricus", "Candida"];
const meta = ["id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tphylum_taxid\tclass_taxid\torder_taxid\tfamily_taxid\tgenus_taxid\tspecies_taxid\tlabel"];
genomes.forEach((g, i) => {
  const genus = genera[i % 4];
  const phylum = i % 4 === 2 ? "Basidiomycota" : "Ascomycota";
  fs.writeFileSync(path.join(tmp, `g${i}.fa`), fasta(`g${i}`, g));
  meta.push([`g${i}`, 1000 + i, phylum, `C${i % 4}`, `O${i % 4}`, `F${i % 4}`, genus, `${genus} sp${i}`, phylum === "Ascomycota" ? 4890 : 5204, 0, 0, 0, 2000 + (i % 4), 1000 + i, `label ${i}`].join("\t"));
});
fs.writeFileSync(path.join(tmp, "meta.tsv"), meta.join("\n") + "\n");
const run = (...args) => execFileSync(gmg, args, { encoding: "utf8", stdio: ["ignore", "pipe", "inherit"] });
run("sketch", "-k", "21", "-s", "10", "-o", path.join(tmp, "sk"), ...genomes.map((_, i) => path.join(tmp, `g${i}.fa`)));
run("build", "--meta", path.join(tmp, "meta.tsv"), "--sketches", path.join(tmp, "sk"), "-o", path.join(tmp, "raw.gmgi"));
run("convert", path.join(tmp, "raw.gmgi"), "-o", path.join(tmp, "ix.gmgi"));

// the inputs: 14 assemblies (exact copies, mutants at several rates, divergent, unrelated, a tiny piece), 2 read sets, 4 bad files
const inp = path.join(tmp, "in");
fs.mkdirSync(inp);
const files = {};
genomes.slice(0, 8).forEach((g, i) => (files[`exact_${i}.fa`] = fasta(`e${i}`, g)));
[0.005, 0.02, 0.05, 0.10, 0.15].forEach((r, i) => (files[`mutant_${i}.fa`] = fasta(`m${i}`, mutate(genomes[i + 1], r))));
files["unrelated.fa"] = fasta("u", randGenome(60_000));
files["piece.fa"] = fasta("p", genomes[3].slice(0, 5_000));
files["reads_a.fastq"] = fastq(genomes[2], 3000, 100);
files["reads_b.fastq"] = fastq(genomes[5], 3000, 100);
files["empty.fa"] = "";
files["header_only.fa"] = ">only_a_header\n";
files["notes.fa"] = "this is a text file, not a genome\n";
files["tiny.fa"] = ">t\nACGTACGT\n";
for (const [n, t] of Object.entries(files)) fs.writeFileSync(path.join(inp, n), t);
assert.ok(Object.keys(files).length >= 20);

const out = path.join(tmp, "out");
run("run", inp, "--index", path.join(tmp, "ix.gmgi"), "-o", out, "--no-report", "-q");
// results.jsonl is written as inputs finish (any order); results.tsv, hits.tsv and the report use the input order (sorted by id)
const stem = (n) => n.replace(/\.(fa|fastq)$/, "");
const names = Object.keys(files).sort((a, b) => (stem(a) < stem(b) ? -1 : 1));
const byId = Object.fromEntries(fs.readFileSync(path.join(out, "results.jsonl"), "utf8").trim().split("\n").map((l) => JSON.parse(l)).map((r) => [r.id, r]));
const cliRecs = names.map((n) => byId[stem(n)]);
assert.ok(cliRecs.every(Boolean), "the CLI produced a record for every input");

const ix = new Index(new Uint8Array(fs.readFileSync(path.join(tmp, "ix.gmgi"))));
const wasmRecs = names.map((n) => JSON.parse(ix.classify(stem(n), n, new Uint8Array(Buffer.from(files[n])), "auto", 0, 0n, 5)));
// wasm and native libm may differ in the last bit of a float (powf): compare numbers to 1e-12 relative
const close = (a, b, where) => {
  if (typeof a === "number" && typeof b === "number") {
    assert.ok(a === b || Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(a), Math.abs(b)), `${where}: ${a} vs ${b}`);
  } else if (a && b && typeof a === "object" && typeof b === "object") {
    assert.deepEqual(Object.keys(a).sort(), Object.keys(b).sort(), `${where}: keys`);
    for (const k of Object.keys(a)) close(a[k], b[k], `${where}.${k}`);
  } else {
    assert.equal(a, b, where);
  }
};
cliRecs.forEach((r, i) => {
  r.path = path.basename(r.path);
  if (r.guess) r.guess.query = path.basename(r.guess.query);
  close(wasmRecs[i], r, `record ${r.id}`);
});
// the TSVs are byte for byte the same apart from the path column
const jsonl = (recs) => recs.map((r) => JSON.stringify(r)).join("\n") + "\n";
const norm = (t) => t.trim().split("\n").map((l) => l.split("\t").filter((_, i) => i !== 1).join("\t")).join("\n");
assert.equal(norm(make_results_tsv(jsonl(wasmRecs))), norm(fs.readFileSync(path.join(out, "results.tsv"), "utf8")), "results.tsv differs");
assert.equal(make_hits_tsv(jsonl(wasmRecs)).trim(), fs.readFileSync(path.join(out, "hits.tsv"), "utf8").trim(), "hits.tsv differs");
// relationships among the assemblies match the CLI's relationships.json
const relCli = JSON.parse(fs.readFileSync(path.join(out, "relationships.json"), "utf8"));
const relWasm = JSON.parse(ix.relationships(jsonl(wasmRecs)));
assert.deepEqual(relWasm, relCli, "relationships differ");
// the statuses cover the cases that matter
const st = Object.fromEntries(wasmRecs.map((r) => [r.id, r.status]));
assert.equal(st.exact_0, "ok");
assert.equal(st.unrelated, "no_match");
assert.equal(st.empty, "error");
assert.equal(wasmRecs.find((r) => r.id === "empty").error_code, "empty_file");
assert.equal(wasmRecs.find((r) => r.id === "notes").error_code, "not_sequence_file");
assert.equal(wasmRecs.find((r) => r.id === "header_only").error_code, "no_sequence_records");
assert.equal(wasmRecs.find((r) => r.id === "tiny").error_code, "too_short");
assert.equal(wasmRecs.find((r) => r.id === "reads_a").kind, "reads");
// taxids and labels are in the records
assert.equal(wasmRecs.find((r) => r.id === "exact_0").guess.calls.species.taxid, 1000);
assert.equal(wasmRecs.find((r) => r.id === "exact_0").guess.top.label, "label 0");
// a report renders from what the browser has
const html = render_report(jsonl(wasmRecs), ix.run_json("2026-10-29T00:00:00Z", 1.5), JSON.stringify(relWasm));
assert.ok(html.includes("Classification report") && html.includes("Relationships among inputs"));
// threshold sets are selectable and recorded
ix.use_thresholds("novel");
assert.equal(JSON.parse(ix.run_json("", 0)).thresholds.name, "novel");
assert.throws(() => ix.use_thresholds("nope"));
ix.free();
fs.rmSync(tmp, { recursive: true });
console.log(`gmg-wasm parity OK: ${names.length} inputs give identical records, results.tsv, hits.tsv and relationships in wasm and the CLI`);
