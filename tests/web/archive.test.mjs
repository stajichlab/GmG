// node --test tests/web/archive.test.mjs  (Node 20+: Blob, File, DecompressionStream and Response are built in)
import { test } from "node:test";
import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { mkdtempSync, mkdirSync, writeFileSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { gzipSync } from "node:zlib";
import { expandArchive, expandArchiveInWorker, KEEP } from "../../web/archive.js";

const dir = mkdtempSync(join(tmpdir(), "gmg-arch-"));
const src = join(dir, "src");
mkdirSync(join(src, "a"), { recursive: true });
mkdirSync(join(src, "b"), { recursive: true });
const big = ">x\n" + "ACGT".repeat(40000) + "\n"; // compresses, spans several 512-byte blocks and stream chunks
writeFileSync(join(src, "a", "genome.fa"), big);
writeFileSync(join(src, "b", "genome.fa"), ">y\nTTTTGGGGCCCCAAAA\n");
writeFileSync(join(src, "reads.fastq"), "@r\nACGT\n+\nIIII\n");
writeFileSync(join(src, "notes.txt"), "hello\n");
const keep = (n) => /\.(fa|fna|fasta|fas|ffn|fsa|seq|dna|fq|fastq)(\.gz)?$/i.test(n);
const asFile = (path, name) => new File([readFileSync(path)], name);
const sh = (cmd, args, cwd = src) => execFileSync(cmd, args, { cwd });

sh("zip", ["-q", "-r", join(dir, "d.zip"), "."]);
sh("zip", ["-q", "-r", "-0", join(dir, "s.zip"), "."]); // stored
sh("tar", ["-cf", join(dir, "t.tar"), "."]);
sh("tar", ["-czf", join(dir, "t.tar.gz"), "."]);

for (const [file, label] of [["d.zip", "zip (deflate)"], ["s.zip", "zip (stored)"], ["t.tar", "tar"], ["t.tar.gz", "tar.gz"]]) {
  test(label, async () => {
    const { files, skipped } = await expandArchive(asFile(join(dir, file), file), keep);
    assert.equal(files.length, 3);
    assert.equal(skipped.notSeq, 1); // notes.txt
    const names = files.map((f) => f.name).sort();
    assert.deepEqual(names, ["a_genome.fa", "b_genome.fa", "reads.fastq"]); // the two genome.fa get their folder
    const a = files.find((f) => f.name === "a_genome.fa");
    assert.equal(await a.text(), big);
    assert.match(a.innerPath, /a\/genome\.fa$/);
    assert.equal(await files.find((f) => f.name === "reads.fastq").text(), "@r\nACGT\n+\nIIII\n");
  });
}

test("a long path (GNU/pax long name) survives", async () => {
  const deep = join(src, "very".repeat(20), "deeper".repeat(20));
  mkdirSync(deep, { recursive: true });
  writeFileSync(join(deep, "long_named_genome.fna"), ">z\nACGTACGT\n");
  sh("tar", ["-cf", join(dir, "long.tar"), "."]);
  const { files } = await expandArchive(asFile(join(dir, "long.tar"), "long.tar"), keep);
  const f = files.find((x) => x.name === "long_named_genome.fna");
  assert.ok(f && f.innerPath.length > 100, "long path kept");
  assert.equal(await f.text(), ">z\nACGTACGT\n");
});

test("a damaged archive gives a readable error", async () => {
  await assert.rejects(expandArchive(new File([new Uint8Array(100)], "bad.zip"), keep), /bad\.zip: not a zip file/);
  const t = readFileSync(join(dir, "t.tar")).subarray(0, 700);
  await assert.rejects(expandArchive(new File([t], "cut.tar"), keep), /cut\.tar: the archive ends inside a file/);
});

test("an empty tar and a zip of nothing useful return no files", async () => {
  const { files } = await expandArchive(new File([new Uint8Array(1024)], "empty.tar"), keep);
  assert.equal(files.length, 0);
});

test("a .tar.gz of several gzip members, and one with padding after the data, unpack (Chrome's DecompressionStream refuses both)", async () => {
  const tar = readFileSync(join(dir, "t.tar"));
  const half = Math.floor(tar.length / 2);
  const multi = Buffer.concat([gzipSync(tar.subarray(0, half)), gzipSync(tar.subarray(half))]);
  const padded = Buffer.concat([gzipSync(tar), Buffer.alloc(1536)]);
  for (const [name, bytes] of [["multi.tar.gz", multi], ["padded.tgz", padded]]) {
    const { files } = await expandArchive(new File([bytes], name), keep);
    assert.deepEqual(files.map((f) => f.name).sort(), ["genome.fa", "genome.fa", "reads.fastq"].map((n, i) => (i < 2 ? ["a_genome.fa", "b_genome.fa"][i] : n)).sort(), name);
    assert.equal(await files.find((f) => f.name === "reads.fastq").text(), "@r\nACGT\n+\nIIII\n");
  }
});

test("a damaged .tar.gz is an error that names the file: a flipped CRC, and a file cut short", async () => {
  const good = readFileSync(join(dir, "t.tar.gz"));
  const badCrc = Buffer.from(good); badCrc[badCrc.length - 8] ^= 1;
  await assert.rejects(expandArchive(new File([badCrc], "crc.tar.gz"), keep), /crc\.tar\.gz: gzip checksum mismatch/);
  await assert.rejects(expandArchive(new File([good.subarray(0, good.length - 20)], "cut.tar.gz"), keep), /cut\.tar\.gz: .*(ends early|ends inside)/);
});

test("the kind of file to keep can be named, and expandArchiveInWorker works without a Worker (node)", async () => {
  assert.ok(KEEP.seq.test("x/genome.fa.gz") && !KEEP.seq.test("notes.txt") && KEEP.bundle.test("a.gsk") && !KEEP.bundle.test("a.fa"));
  const { files } = await expandArchiveInWorker(asFile(join(dir, "t.tar.gz"), "t.tar.gz"), "seq");
  assert.ok(files.some((f) => f.name === "reads.fastq"));
  const viaKind = await expandArchive(asFile(join(dir, "t.tar"), "t.tar"), "seq");
  assert.equal(viaKind.files.length, 3);
});
