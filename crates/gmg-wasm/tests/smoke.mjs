// Parity test: the wasm build must give the same answer as the native CLI on the same query.
// Builds a tiny synthetic reference (6 random genomes + one 2%-mutated copy of genome 0 as the query), runs both.
// Run with `pixi run wasm-test` (needs target/release/gmg and crates/gmg-wasm/pkg-node).
import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import assert from "node:assert/strict";

const require = createRequire(import.meta.url);
const here = path.dirname(new URL(import.meta.url).pathname);
const { Index } = require(path.join(here, "../pkg-node/gmg_wasm.js"));
const gmg = path.join(here, "../../../target/release/gmg");

let seed = 12345;
const rnd = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 2 ** 32);
const randGenome = (n) => Array.from({ length: n }, () => "ACGT"[Math.floor(rnd() * 4)]).join("");
const mutate = (s, rate) => [...s].map((c) => (rnd() < rate ? "ACGT"[Math.floor(rnd() * 4)] : c)).join("");
const fasta = (name, seq) => `>${name}\n${seq.match(/.{1,70}/g).join("\n")}\n`;

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "gmg-smoke-"));
const genomes = Array.from({ length: 6 }, () => randGenome(200_000));
const lineage = [
  ["Ascomycota", "Eurotiomycetes", "Eurotiales", "Aspergillaceae", "Aspergillus", "Aspergillus alpha"],
  ["Ascomycota", "Eurotiomycetes", "Eurotiales", "Aspergillaceae", "Aspergillus", "Aspergillus beta"],
  ["Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae", "Fusarium", "Fusarium gamma"],
  ["Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae", "Fusarium", "Fusarium delta"],
  ["Basidiomycota", "Agaricomycetes", "Agaricales", "Agaricaceae", "Agaricus", "Agaricus eps"],
  ["Basidiomycota", "Agaricomycetes", "Agaricales", "Agaricaceae", "Agaricus", "Agaricus zeta"],
];
const meta = ["id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies"];
genomes.forEach((g, i) => {
  fs.writeFileSync(path.join(tmp, `g${i}.fa`), fasta(`g${i}`, g));
  meta.push([`g${i}`, 100 + i, ...lineage[i]].join("\t"));
});
fs.writeFileSync(path.join(tmp, "meta.tsv"), meta.join("\n") + "\n");
const queryFa = fasta("query", mutate(genomes[2], 0.02));
fs.writeFileSync(path.join(tmp, "query.fa"), queryFa);

const run = (...args) => execFileSync(gmg, args, { encoding: "utf8", stdio: ["ignore", "pipe", "inherit"] });
run("sketch", "-k", "21", "-s", "10", "-o", path.join(tmp, "sk"), ...genomes.map((_, i) => path.join(tmp, `g${i}.fa`)));
run("build", "--meta", path.join(tmp, "meta.tsv"), "--sketches", path.join(tmp, "sk"), "-o", path.join(tmp, "ix.gmgi"));
run("convert", path.join(tmp, "ix.gmgi"), "-o", path.join(tmp, "ix.v1"));
const cli = JSON.parse(run("guess", "--index", path.join(tmp, "ix.gmgi"), path.join(tmp, "query.fa")));

for (const file of ["ix.gmgi", "ix.v1"]) {   // both index formats
  const ix = new Index(new Uint8Array(fs.readFileSync(path.join(tmp, file))));
  assert.equal(ix.ksize(), 21);
  assert.equal(ix.scaled(), 10);
  assert.equal(ix.n_genomes(), 6);
  const wasm = JSON.parse(ix.guess_sequence(new Uint8Array(Buffer.from(queryFa)), 1, 0n, 5));
  assert.equal(wasm.query_hashes, cli.query_hashes, `${file}: query hashes differ`);
  assert.deepEqual(wasm.top, cli.top, `${file}: top hit differs`);
  assert.deepEqual(wasm.calls, cli.calls, `${file}: calls differ`);
  assert.equal(wasm.top.id, "g2");
  // a 2%-mutated copy has ANI ~0.98, below the calibrated species threshold (0.994): the species call is withheld,
  // the raw vote still names the right species, and the genus call is made
  assert.equal(wasm.calls.species.taxon, null);
  assert.equal(wasm.calls.species.vote, "Fusarium gamma");
  assert.equal(wasm.calls.genus.taxon, "Fusarium");
  ix.free();
}
// bad input must throw, not abort
assert.throws(() => new Index(new Uint8Array([1, 2, 3])));
fs.rmSync(tmp, { recursive: true });
console.log("gmg-wasm smoke test OK: wasm and native CLI agree (v0 and v1 index)");
