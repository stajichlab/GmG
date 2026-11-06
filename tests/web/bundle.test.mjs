// node --test tests/web/bundle.test.mjs : the bundle packer (scripts/popsketch_bundle.py) and readBundle / fetchBundle (web/bundle.js).
// The sketch files here are placeholders: readBundle checks the bundle's structure, popsketch-wasm parses the sketches (parity.mjs).
import { test } from "node:test";
import assert from "node:assert/strict";
import { execFileSync, spawnSync } from "node:child_process";
import { mkdtempSync, mkdirSync, writeFileSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { readBundle, fetchBundle } from "../../web/bundle.js";

const root = join(import.meta.dirname, "../..");
const packer = join(root, "scripts/popsketch_bundle.py");
const dir = mkdtempSync(join(tmpdir(), "psbundle-"));
const sk = join(dir, "sk");
mkdirSync(sk);
writeFileSync(join(sk, "params.txt"), "k=21;scaled=200;max_bases=1000000000;r1_only=false;target_cov=none");
const ids = ["SRR1", "SRR2", "SRR3", "SRR4"];
for (const id of ids) {
  writeFileSync(join(sk, `${id}.json`), JSON.stringify({ id, kind: "reads", reads: 10, bases: 1000, non_acgt: 0, files: [] }));
  writeFileSync(join(sk, `${id}.gsk`), Buffer.from(`whole ${id}`));
  if (id !== "SRR4") { // SRR4 has no halves
    writeFileSync(join(sk, `${id}.half_a.gsk`), Buffer.from(`a ${id}`));
    writeFileSync(join(sk, `${id}.half_b.gsk`), Buffer.from(`b ${id}`));
  }
}
writeFileSync(join(dir, "labels.tsv"), "id\tcountry\nSRR1\tUSA\nSRR2\tChina\n");
const pack = (out, ...extra) => spawnSync("python3", [packer, sk, "-o", join(dir, out), ...extra], { encoding: "utf8" });
const asFile = (path, name) => new File([readFileSync(path)], name);

test("a packed bundle reads back with its parameters, labels and the bytes of every sketch file", async () => {
  const r = pack("b.tar.gz", "--labels", join(dir, "labels.tsv"), "--note", "test set");
  assert.equal(r.status, 0, r.stderr);
  assert.match(r.stderr, /4 samples \(3 with split halves\), k=21 scaled=200/); // all four are within the default 20
  const b = await readBundle(asFile(join(dir, "b.tar.gz"), "b.tar.gz"));
  assert.equal(b.k, 21);
  assert.equal(b.scaled, 200);
  assert.equal(b.note, "test set");
  assert.match(b.labels, /^id\tcountry\nSRR1\tUSA/);
  assert.deepEqual(b.sketches.map((s) => s.id), ids);
  const s1 = b.sketches[0], s4 = b.sketches[3];
  assert.equal(new TextDecoder().decode(s1.whole), "whole SRR1");
  assert.equal(new TextDecoder().decode(s1.a), "a SRR1");
  assert.equal(new TextDecoder().decode(s1.b), "b SRR1");
  assert.equal(s4.a.length + s4.b.length, 0); // no halves: empty arrays, which add_sketch reads as none
  assert.equal(JSON.parse(new TextDecoder().decode(s1.meta)).id, "SRR1");
});

test("packing the same sketches twice gives the same members in the same order", () => {
  pack("x1.tar.gz");
  const r = pack("x2.tar.gz");
  assert.equal(r.status, 0);
  // (the manifest records the creation time, so the archives are not byte-identical)
  const list = (f) => execFileSync("tar", ["-tzf", join(dir, f)]).toString();
  assert.equal(list("x1.tar.gz"), list("x2.tar.gz"));
});

test("a sketch without its .gsk is refused by the packer, a non-bundle and a missing member by the reader", async () => {
  writeFileSync(join(sk, "SRR5.json"), JSON.stringify({ id: "SRR5", kind: "reads", reads: 1, bases: 1 }));
  const r = pack("bad.tar.gz");
  assert.notEqual(r.status, 0);
  assert.match(r.stderr, /SRR5.*without SRR5\.gsk/);
  rmSync(join(sk, "SRR5.json"));
  // a tar that is not a bundle
  execFileSync("tar", ["-czf", join(dir, "plain.tar.gz"), "labels.tsv"], { cwd: dir });
  await assert.rejects(readBundle(asFile(join(dir, "plain.tar.gz"), "plain.tar.gz")), /not a popsketch bundle/);
  // a manifest that names a sample whose files are not in the archive
  const bd = join(dir, "broken");
  mkdirSync(bd);
  writeFileSync(join(bd, "popsketch-bundle.json"), JSON.stringify({ format: 1, k: 21, scaled: 200, samples: ["SRR1"] }));
  execFileSync("tar", ["-czf", join(dir, "broken.tar.gz"), "popsketch-bundle.json"], { cwd: bd });
  await assert.rejects(readBundle(asFile(join(dir, "broken.tar.gz"), "broken.tar.gz")), /SRR1\.json is listed in the bundle but missing/);
  writeFileSync(join(bd, "popsketch-bundle.json"), JSON.stringify({ format: 2, k: 21, scaled: 200, samples: ["SRR1"] }));
  execFileSync("tar", ["-czf", join(dir, "v2.tar.gz"), "popsketch-bundle.json"], { cwd: bd });
  await assert.rejects(readBundle(asFile(join(dir, "v2.tar.gz"), "v2.tar.gz")), /bundle format 2/);
});

test("fetchBundle reads the kind of archive from the bytes, and reports an HTTP error", async () => {
  pack("f.tar.gz");
  const data = readFileSync(join(dir, "f.tar.gz"));
  const f = await fetchBundle("https://example.org/download?id=7", async () => new Response(data));
  assert.match(f.name, /\.tar\.gz$/);
  assert.equal((await readBundle(f)).sketches.length, 4);
  await assert.rejects(fetchBundle("https://example.org/x", async () => new Response("no", { status: 404, statusText: "Not Found" })), /404 Not Found/);
});
