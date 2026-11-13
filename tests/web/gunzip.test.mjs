// node tests/web/gunzip.test.mjs (and `node --test tests/web/`): the multi-member gzip decoder against zlib, with the chunkings and oddities that break a decoder.
import assert from "node:assert/strict";
import zlib from "node:zlib";
import { makeGunzip, parseHeader } from "../../web/gunzip.js";

const pako = (await import("../../web/vendor/pako_inflate.mjs")).default;
const gunzip = makeGunzip(pako);

async function decode(bytes, chunk) {
  const rs = new ReadableStream({
    start(c) {
      for (let i = 0; i < bytes.length; i += chunk) c.enqueue(bytes.subarray(i, i + chunk));
      c.close();
    },
  });
  const out = [];
  for await (const v of rs.pipeThrough(gunzip())) out.push(Buffer.from(v));
  return Buffer.concat(out).toString();
}

const fq = (n, tag) => Array.from({ length: n }, (_, i) => `@${tag}${i}\nACGTACGTAC\n+\nIIIIIIIIII\n`).join("");
const A = fq(500, "a"), B = fq(300, "b"), C = fq(1, "c");
const gz = (s, o) => zlib.gzipSync(Buffer.from(s), o);
let n = 0;
const ok = (name) => { n++; console.log("ok:", name); };

for (const chunk of [1, 3, 7, 64, 1 << 16]) {
  assert.equal(await decode(gz(A), chunk), A);
  assert.equal(await decode(Buffer.concat([gz(A), gz(B), gz(C)]), chunk), A + B + C);
  // trailing zero padding (tape or some pipelines) and trailing junk are ignored, as gzip does
  assert.equal(await decode(Buffer.concat([gz(A), Buffer.alloc(37)]), chunk), A);
  assert.equal(await decode(Buffer.concat([gz(A), gz(B), Buffer.from("trailing text")]), chunk), A + B);
  // an empty member in the middle
  assert.equal(await decode(Buffer.concat([gz(A), gz(""), gz(B)]), chunk), A + B);
}
ok("single, multi-member, padding, junk and an empty member, at chunk sizes 1 to 64 KiB");

// a header with every optional field (FEXTRA, FNAME, FCOMMENT, FHCRC), split anywhere
const fancy = (() => {
  const body = zlib.deflateRawSync(Buffer.from(A));
  const crc = Buffer.alloc(8); // the trailer is not checked
  const hdr = Buffer.concat([
    Buffer.from([0x1f, 0x8b, 8, 4 | 8 | 16 | 2, 0, 0, 0, 0, 0, 3]),
    Buffer.from([5, 0, 1, 2, 3, 4, 5]), Buffer.from("name.fq\0"), Buffer.from("a comment\0"), Buffer.from([0xab, 0xcd]),
  ]);
  return Buffer.concat([hdr, body, crc]);
})();
for (const chunk of [1, 2, 5, 11, 4096]) assert.equal(await decode(Buffer.concat([fancy, gz(B)]), chunk), A + B);
assert.equal(parseHeader(fancy), 10 + 7 + 8 + 10 + 2);
ok("a header with FEXTRA, FNAME, FCOMMENT and FHCRC, with a second member after it");

// bgzip-like: many small members
const many = Buffer.concat(Array.from({ length: 200 }, (_, i) => gz(fq(2, "r" + i))));
assert.equal(await decode(many, 1000), Array.from({ length: 200 }, (_, i) => fq(2, "r" + i)).join(""));
ok("200 small members");

// errors: not gzip, truncated, corrupt
await assert.rejects(decode(Buffer.from("this is not gzip at all"), 8), /not a gzip file/);
await assert.rejects(decode(gz(A).subarray(0, 40), 8), /ends early/);
const bad = Buffer.from(gz(A)); bad[20] ^= 0xff; bad[21] ^= 0xff; bad[22] ^= 0xff;
await assert.rejects(decode(bad, 64), /corrupt|ends early/);
assert.equal(await decode(Buffer.alloc(0), 8), "");
ok("not gzip, truncated and corrupt data are errors; an empty file is empty");

// what the native decoder does with the same multi-member bytes, for the record
let native = "ok";
try {
  const rs = new ReadableStream({ start(c) { c.enqueue(Buffer.concat([gz(A), gz(B)])); c.close(); } }).pipeThrough(new DecompressionStream("gzip"));
  const out = []; for await (const v of rs) out.push(Buffer.from(v));
  native = Buffer.concat(out).toString() === A + B ? "reads every member" : "stops after the first member";
} catch (e) { native = `throws: ${e.message}`; }
console.log(`(node's DecompressionStream on two members: ${native})`);
console.log(`all ${n} gunzip tests passed`);
