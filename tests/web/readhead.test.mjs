// node --test tests/web/readhead.test.mjs : readHead / readWhole (web/readhead.js) and openDecoded (web/gunzip-browser.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import zlib from "node:zlib";
import { readHead, readWhole } from "../../web/readhead.js";
import { openDecoded } from "../../web/gunzip-browser.js";

// a File-like object whose stream() delivers the bytes in chunks of a chosen size, so chunk boundaries can be controlled
function fileOf(name, bytes, chunk = 1 << 16) {
  const buf = Buffer.from(bytes);
  return { name, size: buf.length, stream: () => new ReadableStream({ start(c) { for (let i = 0; i < buf.length; i += chunk) c.enqueue(new Uint8Array(buf.subarray(i, i + chunk))); c.close(); } }) };
}
const fq = (n, tag = "r") => Array.from({ length: n }, (_, i) => `@${tag}${i}\nACGTACGTAC\n+\nIIIIIIIIII\n`).join("");
const text = (u8) => Buffer.from(u8).toString();
const lines = (s) => s.split("\n").length - 1;

test("exactly maxReads records are returned, at any chunk size, plain and gzip", async () => {
  const data = fq(2500);
  for (const [name, bytes] of [["a.fastq", data], ["a.fastq.gz", zlib.gzipSync(data)]]) {
    for (const chunk of [1, 7, 100, 4096, 1 << 20]) {
      for (const reads of [1, 2, 999, 1000, 2500, 5000]) {
        const got = text(await readHead(fileOf(name, bytes, chunk), reads));
        assert.equal(lines(got), 4 * Math.min(reads, 2500), `${name} chunk ${chunk} reads ${reads}`);
        assert.equal(got, data.slice(0, got.length));
      }
    }
  }
});

test("a read limit that is not a positive whole number is refused, not turned into a few reads", async () => {
  for (const bad of [0, -1, 1.5, NaN, "", undefined]) await assert.rejects(readHead(fileOf("a.fastq", fq(10)), bad), /whole number of 1 or more/);
});

test("gzip is recognised by its first bytes, not by the file name", async () => {
  const data = fq(50);
  assert.equal(text(await readHead(fileOf("renamed.fastq", zlib.gzipSync(data), 1), 50)), data); // gzip called .fastq, 1-byte chunks
  assert.equal(text(await readHead(fileOf("plain.fastq.gz", data, 1), 50)), data); // plain text called .gz
  const s = await openDecoded(fileOf("x", Buffer.from([0x1f]), 1)); // a one-byte file: not gzip, passed through
  const r = s.getReader();
  assert.deepEqual([...(await r.read()).value], [0x1f]);
});

test("a multi-member file is read through its members", async () => {
  const gz = (s) => zlib.gzipSync(s);
  const bytes = Buffer.concat([gz(fq(100, "a")), gz(fq(100, "b")), gz(fq(100, "c"))]);
  const got = text(await readHead(fileOf("m.fastq.gz", bytes, 333), 250));
  assert.equal(lines(got), 1000);
  assert.ok(got.includes("@a99\n") && got.includes("@b0\n") && got.includes("@c49\n") && !got.includes("@c50\n"));
});

test("an error beyond the read limit is never seen, and one inside the data that is read is", async () => {
  const good = zlib.gzipSync(fq(5000));
  const damaged = Buffer.from(good);
  damaged[damaged.length - 6] ^= 0xff; // the CRC at the end of the member: beyond a limit of 100 reads
  const got = text(await readHead(fileOf("a.fastq.gz", damaged, 1 << 20), 100));
  assert.equal(lines(got), 400);
  await assert.rejects(readHead(fileOf("a.fastq.gz", damaged, 1 << 20), 6000), /checksum mismatch/); // the file ends first: read to the end, checked
  // (a limit equal to the number of reads in the file stops at the last record, before the trailer: not checked, as documented)
  assert.equal(lines(text(await readHead(fileOf("a.fastq.gz", damaged, 1 << 20), 5000))), 20000);
  const truncated = good.subarray(0, 600);
  assert.equal(lines(text(await readHead(fileOf("t.fastq.gz", truncated, 1 << 20), 20))), 80); // enough reads before the cut
  await assert.rejects(readHead(fileOf("t.fastq.gz", truncated, 1 << 20), 6000), /ends early/);
});

test("the bytes limit cuts a FASTQ file at the end of a record, with a note; a FASTA file that large is an error", async () => {
  const notes = [];
  const got = text(await readHead(fileOf("big.fastq", fq(1000), 500), 1000, { maxBytes: 10000, onNote: (m) => notes.push(m) }));
  assert.ok(got.length <= 10500 && lines(got) % 4 === 0 && lines(got) > 0, `${got.length} bytes, ${lines(got)} lines`);
  assert.match(notes.join(), /stopped at \d+ reads/);
  await assert.rejects(readHead(fileOf("g.fa", ">x\n" + "ACGT".repeat(5000) + "\n", 700), 10, { maxBytes: 5000 }), /too large for the browser/);
});

test("a line with no line break for too long is not FASTQ", async () => {
  await assert.rejects(readHead(fileOf("n.fastq", "@" + "A".repeat(10000), 1000), 10, { maxLine: 4000 }), /no line break/);
});

test("a FASTA file is read whole, and notes about bytes after the last gzip member reach onNote", async () => {
  const fa = ">c1\n" + "ACGT".repeat(500) + "\n>c2\nTTTT\n";
  assert.equal(text(await readHead(fileOf("g.fa", fa, 100), 10)), fa);
  const notes = [];
  const padded = Buffer.concat([zlib.gzipSync(fq(10)), Buffer.alloc(100)]);
  assert.equal(lines(text(await readHead(fileOf("p.fastq.gz", padded, 64), 1000, { onNote: (m) => notes.push(m) }))), 40);
  assert.match(notes.join(), /100 bytes after the last gzip member were ignored/);
});

test("an empty file says so, whatever its name", async () => {
  await assert.rejects(readHead(fileOf("e.fastq.gz", Buffer.alloc(0)), 10), /e\.fastq\.gz: empty file/);
  await assert.rejects(readHead(fileOf("e.fastq", Buffer.alloc(0)), 10), /empty file/);
});

test("readWhole checks the CRC of a whole file", async () => {
  const data = fq(300);
  assert.equal(text(await readWhole(fileOf("a.fastq.gz", zlib.gzipSync(data), 77))), data);
  const damaged = Buffer.from(zlib.gzipSync(data));
  damaged[damaged.length - 7] ^= 1;
  await assert.rejects(readWhole(fileOf("a.fastq.gz", damaged, 77)), /checksum mismatch/);
});
