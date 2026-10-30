// #49: WASM against the native CLI on the REAL golden-panel genomes (the CI parity test uses synthetic ones).
// usage: node scripts/wasm_golden.mjs PKG_NODE_DIR GMG_BIN INDEX.gmgi PANEL_CACHE_DIR
// Classifies every <id>.fna.gz of the cache in wasm and with `gmg run`, and requires equal records (numbers to 1e-12 relative: libm may
// differ in the last bit between wasm and native).
import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import zlib from "node:zlib";
import assert from "node:assert/strict";

const [pkg, gmg, indexPath, cache] = process.argv.slice(2);
const require = createRequire(import.meta.url);
const { Index } = require(path.resolve(pkg, "gmg_wasm.js"));
const files = fs.readdirSync(cache).filter((f) => f.endsWith(".fna.gz")).sort();
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "gmg-wasm-golden-"));
const inp = path.join(tmp, "in");
fs.mkdirSync(inp);
for (const f of files) fs.symlinkSync(path.resolve(cache, f), path.join(inp, f));
execFileSync(gmg, ["-q", "run", inp, "--index", indexPath, "-o", path.join(tmp, "out"), "--no-report"], { stdio: ["ignore", "inherit", "inherit"] });
const native = Object.fromEntries(
  fs.readFileSync(path.join(tmp, "out", "results.jsonl"), "utf8").trim().split("\n").map((l) => JSON.parse(l)).map((r) => [r.id, r]),
);
const ix = new Index(new Uint8Array(fs.readFileSync(indexPath)));
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
const skip = new Set(["path", "size_bytes", "files", "query"]); // path and size differ by construction (the wasm input is the decompressed text)
const strip = (o) => {
  if (Array.isArray(o)) o.forEach(strip);
  else if (o && typeof o === "object") {
    for (const k of Object.keys(o)) {
      if (skip.has(k) && typeof o[k] !== "object") delete o[k];
      else strip(o[k]);
    }
  }
};
let n = 0;
for (const f of files) {
  const id = f.replace(/\.fna\.gz$/, "");
  const seq = new Uint8Array(zlib.gunzipSync(fs.readFileSync(path.join(cache, f))));
  const w = JSON.parse(ix.classify(id, f, seq, "auto", 0, 0n, 5));
  const c = native[id];
  strip(w);
  strip(c);
  close(w, c, id);
  n++;
}
console.log(`wasm == native on ${n} real genomes`);
