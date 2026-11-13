// A gzip decoder for the browser that reads every member of a file, as `gzip -d` does.
//
// Why not `new DecompressionStream("gzip")`: Chrome's decoder stops with "Junk found after end of compressed data" when anything
// follows the first gzip member. That is what a file made of several members looks like (bgzip, `cat a.gz b.gz`, the windows
// that `popseq fetch` concatenates), and what a .gz with padding after the data looks like. Here the header of each member is
// parsed, its deflate data is inflated with pako (vendor/pako_inflate.min.js, raw mode), the 8 trailer bytes are skipped and the
// next member is read. Bytes after the last member that do not start a gzip member (zeros, junk) are ignored, as gzip does with a
// warning. The CRC and length trailer are not checked.
//
//   const gunzip = makeGunzip(pako);                    // pako: the pako_inflate object (self.pako in the browser)
//   const text = readable.pipeThrough(gunzip());        // a TransformStream: Uint8Array in, Uint8Array out
//
// `parseHeader` and the state machine are plain code with no browser API apart from TransformStream, so test_gunzip.mjs runs
// them in node.

// Length of the gzip header at the start of `d`, -1 if `d` ends before the header does; throws if it is not a gzip header.
export function parseHeader(d) {
  if (d.length < 10) {
    if (d.length >= 1 && d[0] !== 0x1f) throw new Error("not a gzip file");
    if (d.length >= 2 && d[1] !== 0x8b) throw new Error("not a gzip file");
    return -1;
  }
  if (d[0] !== 0x1f || d[1] !== 0x8b) throw new Error("not a gzip file");
  if (d[2] !== 8) throw new Error("gzip with an unknown compression method");
  const flg = d[3];
  let p = 10;
  if (flg & 4) { // FEXTRA: 2-byte length, then that many bytes
    if (d.length < p + 2) return -1;
    p += 2 + (d[p] | (d[p + 1] << 8));
    if (d.length < p) return -1;
  }
  for (const bit of [8, 16]) { // FNAME, FCOMMENT: zero-terminated strings
    if (flg & bit) {
      const z = d.indexOf(0, p);
      if (z < 0) return -1;
      p = z + 1;
    }
  }
  if (flg & 2) { // FHCRC
    p += 2;
    if (d.length < p) return -1;
  }
  return p;
}

function concat(a, b) {
  if (a.length === 0) return b;
  const out = new Uint8Array(a.length + b.length);
  out.set(a, 0);
  out.set(b, a.length);
  return out;
}

export function makeGunzip(pako) {
  return function gunzip() {
    let state = "header"; // header | body | trailer | done
    let members = 0;
    let pending = new Uint8Array(0);
    let inf = null;
    let trailerLeft = 0;

    function run(data, ctl) {
      while (data.length > 0 && state !== "done") {
        if (state === "header") {
          data = concat(pending, data);
          pending = new Uint8Array(0);
          if (members > 0) { // after the first member, anything that is not a gzip header is trailing junk
            if (data[0] !== 0x1f || (data.length > 1 && data[1] !== 0x8b)) { state = "done"; return; }
          }
          const h = parseHeader(data);
          if (h < 0) { pending = data; return; }
          data = data.subarray(h);
          inf = new pako.Inflate({ raw: true });
          inf.onData = (c) => ctl.enqueue(new Uint8Array(c));
          state = "body";
        } else if (state === "body") {
          inf.push(data, false);
          if (inf.err) throw new Error(`gzip data is corrupt: ${inf.msg || inf.err}`);
          if (!inf.ended) return; // all of `data` was consumed, the member goes on in the next chunk
          data = data.subarray(data.length - inf.strm.avail_in);
          inf = null;
          state = "trailer";
          trailerLeft = 8;
        } else { // trailer: CRC32 and ISIZE
          const n = Math.min(trailerLeft, data.length);
          trailerLeft -= n;
          data = data.subarray(n);
          if (trailerLeft === 0) { members++; state = "header"; }
        }
      }
    }

    return new TransformStream({
      transform(chunk, ctl) { run(chunk, ctl); },
      flush() {
        if (state === "body" || (state === "header" && members === 0 && pending.length > 0)) {
          throw new Error("gzip file ends early (truncated or still downloading)");
        }
      },
    });
  };
}
