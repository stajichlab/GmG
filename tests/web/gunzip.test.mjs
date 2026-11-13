// node tests/web/gunzip.test.mjs (and `node --test tests/web/`): the multi-member gzip decoder against zlib, with the chunkings and
// damaged inputs that break a decoder. Every group below was checked to fail when the part of gunzip.js it covers is broken.
import assert from "node:assert/strict";
import zlib from "node:zlib";
import { makeGunzip, parseHeader, crc32 } from "../../web/gunzip.js";

const pako = (await import("../../web/vendor/pako_inflate.mjs")).default;
const gunzip = makeGunzip(pako);

function source(bytes, chunk) {
  return new ReadableStream({
    start(c) {
      for (let i = 0; i < bytes.length; i += chunk) c.enqueue(bytes.subarray(i, i + chunk));
      if (bytes.length === 0) c.enqueue(new Uint8Array(0));
      c.close();
    },
  });
}
async function decode(bytes, chunk, opts) {
  const out = [];
  for await (const v of source(bytes, chunk).pipeThrough(gunzip(opts))) out.push(Buffer.from(v));
  return Buffer.concat(out).toString();
}
async function fails(bytes, chunk, pattern, opts) {
  await assert.rejects(decode(bytes, chunk, opts), pattern, `chunk ${chunk}`);
}

const fq = (n, tag) => Array.from({ length: n }, (_, i) => `@${tag}${i}\nACGTACGTAC\n+\nIIIIIIIIII\n`).join("");
const A = fq(500, "a"), B = fq(300, "b"), C = fq(1, "c");
const gz = (s) => zlib.gzipSync(Buffer.from(s));
let n = 0;
const ok = (name) => { n++; console.log("ok:", name); };
const CHUNKS = [1, 3, 7, 64, 1 << 16, 1 << 20];

// ---------------------------------------------------------------------------------------------------- CRC-32
{
  const data = Buffer.from(fq(2000, "x"));
  assert.equal(crc32(0, data), zlib.crc32(data));
  assert.equal(crc32(0, Buffer.alloc(0)), 0);
  for (const cut of [1, 7, 8, 9, 100, 4093]) assert.equal(crc32(crc32(0, data.subarray(0, cut)), data.subarray(cut)), zlib.crc32(data), `split at ${cut}`);
  for (let len = 0; len < 40; len++) assert.equal(crc32(0, data.subarray(3, 3 + len)), zlib.crc32(data.subarray(3, 3 + len)), `length ${len}`);
  ok("crc32 equals zlib.crc32, whole, split and at every short length and alignment");
}

// ---------------------------------------------------------------------------------------------------- valid input
for (const chunk of CHUNKS) {
  assert.equal(await decode(gz(A), chunk), A);
  assert.equal(await decode(Buffer.concat([gz(A), gz(B), gz(C)]), chunk), A + B + C);
  assert.equal(await decode(Buffer.concat([gz(A), gz(""), gz(B)]), chunk), A + B); // an empty member in the middle
}
ok("single and multi-member files, an empty member, at chunk sizes 1 byte to 1 MiB");

// a header with every optional field, the extra field holding zero bytes (a reader that looks for the first zero would stop in it)
const member = (body, { extra = Buffer.from([0, 1, 0, 2, 0, 3]), name = "name.fq", comment = "a comment" } = {}) => {
  const raw = zlib.deflateRawSync(Buffer.from(body));
  const t = Buffer.alloc(8);
  t.writeUInt32LE(zlib.crc32(Buffer.from(body)), 0);
  t.writeUInt32LE(Buffer.byteLength(body) >>> 0, 4);
  const xl = Buffer.alloc(2); xl.writeUInt16LE(extra.length);
  const hdr = Buffer.concat([Buffer.from([0x1f, 0x8b, 8, 4 | 8 | 16 | 2, 0, 0, 0, 0, 0, 3]), xl, extra, Buffer.from(`${name}\0`), Buffer.from(`${comment}\0`), Buffer.from([0xab, 0xcd])]);
  return Buffer.concat([hdr, raw, t]);
};
const fancy = member(A);
for (const chunk of [1, 2, 5, 11, 4096]) assert.equal(await decode(Buffer.concat([fancy, gz(B)]), chunk), A + B);
assert.equal(parseHeader(fancy), 10 + 2 + 6 + 8 + 10 + 2);
ok("a header with FEXTRA (holding zero bytes), FNAME, FCOMMENT and FHCRC, split anywhere, with a second member after it");

// bgzip: the BGZF layout (FEXTRA with the BC subfield, 64 KiB blocks, an empty EOF block), which has zero bytes in its extra field
function bgzf(text, block = 4096) {
  const buf = Buffer.from(text), blocks = [];
  for (let i = 0; i <= buf.length; i += block) {
    const part = buf.subarray(i, Math.min(i + block, buf.length));
    const raw = zlib.deflateRawSync(part);
    const bsize = 18 + raw.length + 8 - 1;
    const hdr = Buffer.from([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0, bsize & 255, bsize >> 8]);
    const t = Buffer.alloc(8); t.writeUInt32LE(zlib.crc32(part), 0); t.writeUInt32LE(part.length, 4);
    blocks.push(Buffer.concat([hdr, raw, t]));
  }
  return Buffer.concat(blocks);
}
const bg = bgzf(A + B, 997);
for (const chunk of [1, 13, 997, 5000, 1 << 20]) assert.equal(await decode(bg, chunk), A + B);
ok("a bgzf file of many small members with the BC extra field, at several chunk sizes");

// bytes after the last member
const warned = [];
for (const chunk of [1, 5, 4096]) {
  warned.length = 0;
  assert.equal(await decode(Buffer.concat([gz(A), Buffer.alloc(37)]), chunk, { onWarn: (m) => warned.push(m) }), A);
  assert.match(warned.join(), /37 bytes after the last gzip member were ignored/);
  warned.length = 0;
  assert.equal(await decode(Buffer.concat([gz(A), gz(B), Buffer.from("trailing text")]), chunk, { onWarn: (m) => warned.push(m) }), A + B);
  assert.match(warned.join(), /13 bytes/);
}
assert.equal(await decode(gz(A), 64, { onWarn: () => assert.fail("no warning for a clean file") }), A);
ok("zero padding and text after the last member are ignored, with a warning that counts them");

// ---------------------------------------------------------------------------------------------------- damaged input
{
  const good = gz(A);
  // a wrong CRC, a wrong length
  const badCrc = Buffer.from(good); badCrc[badCrc.length - 8] ^= 1;
  const badLen = Buffer.from(good); badLen[badLen.length - 4] ^= 1;
  for (const chunk of CHUNKS) { await fails(badCrc, chunk, /checksum mismatch.*damaged/); await fails(badLen, chunk, /length mismatch.*damaged/); }
  // a flipped bit in a stored (uncompressed) block decodes without a deflate error, so only the CRC can see it
  const stored = zlib.gzipSync(Buffer.from(A), { level: 0 });
  const flipped = Buffer.from(stored); flipped[stored.length / 2 | 0] ^= 0x04;
  await fails(flipped, 4096, /checksum mismatch/);
  // the same file with the check off decodes to wrong bytes without any error, which is what `crc: false` accepts
  assert.notEqual(await decode(flipped, 4096, { crc: false }), A);
  // a damaged second member is found, and the first is delivered before the error
  const two = Buffer.concat([gz(A), badCrc]);
  const got = [];
  await assert.rejects((async () => { for await (const v of source(two, 1 << 20).pipeThrough(gunzip())) got.push(Buffer.from(v)); })(), /member 2/);
  assert.equal(Buffer.concat(got).toString(), A + A, "everything decoded before the error is delivered, including the damaged member's own bytes");
  ok("a bad CRC, a bad length and a flipped bit in a stored block are errors (and are not with crc: false); the error names the member");
}

// every cut of a multi-member file inside a trailer or a header is an error, not a shorter file
{
  const m1 = gz(A), m2 = gz(B), file = Buffer.concat([m1, m2]);
  for (let cut = m1.length - 12; cut < m1.length + 14; cut++) {
    if (cut === m1.length) { assert.equal(await decode(file.subarray(0, cut), 7), A, "cut exactly at the end of member 1 is a complete file"); continue; }
    await fails(file.subarray(0, cut), 7, /ends early/);
    await fails(file.subarray(0, cut), 1 << 20, /ends early/);
  }
  await fails(file.subarray(0, file.length - 1), 4096, /ends early/); // inside the last trailer
  await fails(file.subarray(0, 5), 4096, /ends early/); // inside the first header
  ok("a file cut inside a deflate body, a trailer or the next member's header is an error; cut exactly at a member end it is complete");
}

// what is decoded before the input ends early or turns corrupt is delivered (pako keeps up to 16 KiB until its buffer fills or the member ends)
{
  const big = fq(5000, "z"), g = gz(big);
  for (const chunk of [1 << 20, 97]) {
    const got = [];
    await assert.rejects((async () => { for await (const v of source(g.subarray(0, 700), chunk).pipeThrough(gunzip())) got.push(Buffer.from(v)); })(), /ends early/);
    const text = Buffer.concat(got).toString();
    assert.ok(text.length > 1000 && big.startsWith(text), `chunk ${chunk}: ${text.length} bytes delivered before the truncation error`);
  }
  ok("the bytes decoded before a truncation are delivered, then the error");
}

// not gzip, and the pieces of a header
await fails(Buffer.from("this is not gzip at all"), 8, /not a gzip file/);
await fails(Buffer.from([0x1f]), 1, /ends early/);
await fails(Buffer.from([0x1e, 0x8b]), 1, /not a gzip file/);
await fails(Buffer.from([0x1f, 0x8a, 8, 0, 0, 0, 0, 0, 0, 3]), 1, /not a gzip file/);
await fails(Buffer.from([0x1f, 0x8b, 7, 0, 0, 0, 0, 0, 0, 3]), 3, /unknown compression method/);
await fails(Buffer.from([0x1f, 0x8b, 8, 0xe0, 0, 0, 0, 0, 0, 3]), 3, /reserved flags/);
await fails(Buffer.concat([Buffer.from([0x1f, 0x8b, 8, 8, 0, 0, 0, 0, 0, 3]), Buffer.alloc(70000, 65)]), 4096, /header too long/); // FNAME never ends
assert.throws(() => parseHeader(Buffer.from([0x1f, 0x8b, 9])), /unknown compression method/);
assert.equal(parseHeader(Buffer.from([0x1f, 0x8b, 8])), -1);
ok("not gzip, a short or odd header, an unknown method, reserved flags and a header that never ends are errors");

// corrupt deflate data
{
  const bad = Buffer.from(gz(A));
  for (const i of [20, 21, 22]) bad[i] ^= 0xff;
  await assert.rejects(decode(bad, 64), /corrupt|checksum|ends early/);
  ok("corrupt deflate data is an error");
}

// empty input
await fails(Buffer.alloc(0), 8, /empty file/);
ok("an empty file is 'empty file'");

// an error is raised after the bytes decoded from its chunk have been delivered, so a reader that stops first never sees it
{
  const bad = Buffer.concat([gz(A), Buffer.from(gz(B)).fill(0, 20, 30)]); // member 2 damaged
  const reader = source(bad, 1 << 20).pipeThrough(gunzip()).getReader();
  let text = "";
  const first = await reader.read(); // the bytes of member 1 are there
  text += Buffer.from(first.value).toString();
  assert.ok(text.startsWith("@a0\n"));
  await reader.cancel(); // enough: the damage later in the file is never reported
  ok("the bytes decoded before a damaged member are delivered, and a reader that cancels never sees the error");
}

console.log(`all ${n} gunzip test groups passed`);
