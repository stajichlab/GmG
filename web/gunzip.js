// A gzip decoder for the browser that reads every member of a file, as `gzip -d` does, and checks what it can.
//
// Why not `new DecompressionStream("gzip")`: Chrome's decoder stops with "Junk found after end of compressed data" when anything
// follows the first gzip member. That is what a file of several members looks like (bgzip, `cat a.gz b.gz`, the windows
// that `popseq fetch` concatenates), and what a .gz with padding after the data looks like. Here the header of each member is
// parsed, its deflate data is inflated with pako (vendor/pako_inflate.mjs, raw mode), the 8 trailer bytes are read and the next
// member follows.
//
//   const gunzip = makeGunzip(pako);                              // pako: the pako_inflate object
//   const bytes = readable.pipeThrough(gunzip({ crc: true, onWarn: console.warn }));   // a TransformStream: bytes in, bytes out
//
// What is checked
// * With `crc` (default true) the CRC-32 and the length that end each member are compared with the bytes decoded; a mismatch is an
//   error ("the file is damaged"). Without it (about 40% faster) a damaged download that still decodes gives wrong bytes with no
//   error, so use `crc: false` only where the file is not read to its end anyway: a reader that stops after the first N reads never
//   gets to the trailer of its first member, so there is nothing to compare, and the caller should say so.
// * A file that ends inside a member, inside a trailer, or inside the header of a next member is an error, not a shorter file.
// * An empty input is the error "empty file". A header with reserved flag bits, an unknown compression method or more than 64 KiB
//   of header is not gzip.
// * Bytes after the last member that do not start a gzip member (zero padding, junk) are ignored, as `gzip -d` ignores them with a
//   warning; `onWarn` is called with the number of bytes at the end of the stream.
// * An error found while decoding a chunk is raised on the next chunk or at the end, after the bytes decoded from that chunk have been
//   delivered: a reader that has enough and cancels never sees an error that lies beyond what it read.

// ---------------------------------------------------------------------------------------------------------- CRC-32
// Slicing-by-8: eight tables, eight bytes per step. About 4 times a one-table loop.
const T = (() => {
  const t = new Int32Array(8 * 256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c;
  }
  for (let n = 0; n < 256; n++) {
    let c = t[n];
    for (let k = 1; k < 8; k++) {
      c = t[c & 0xff] ^ (c >>> 8);
      t[k * 256 + n] = c;
    }
  }
  return t;
})();

/** CRC-32 of `buf` continuing from `crc` (0 to start): `crc32(crc32(0, a), b) === crc32(0, a + b)`, as zlib.crc32. */
export function crc32(crc, buf) {
  let c = ~crc;
  let i = 0;
  const n = buf.length;
  for (; i + 8 <= n; i += 8) {
    const a = c ^ (buf[i] | (buf[i + 1] << 8) | (buf[i + 2] << 16) | (buf[i + 3] << 24));
    c = T[1792 + (a & 0xff)] ^ T[1536 + ((a >>> 8) & 0xff)] ^ T[1280 + ((a >>> 16) & 0xff)] ^ T[1024 + (a >>> 24)] ^
      T[768 + buf[i + 4]] ^ T[512 + buf[i + 5]] ^ T[256 + buf[i + 6]] ^ T[buf[i + 7]];
  }
  for (; i < n; i++) c = T[(c ^ buf[i]) & 0xff] ^ (c >>> 8);
  return ~c >>> 0;
}

// ---------------------------------------------------------------------------------------------------------- header
const MAX_HEADER = 65536;

// Length of the gzip header at the start of `d`, -1 if `d` ends before the header does; throws if it is not a gzip header.
export function parseHeader(d) {
  if (d.length >= 1 && d[0] !== 0x1f) throw new Error("not a gzip file");
  if (d.length >= 2 && d[1] !== 0x8b) throw new Error("not a gzip file");
  if (d.length >= 3 && d[2] !== 8) throw new Error("gzip with an unknown compression method");
  if (d.length >= 4 && d[3] & 0xe0) throw new Error("gzip header with reserved flags: not a gzip file");
  if (d.length < 10) return -1;
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

// ---------------------------------------------------------------------------------------------------------- decoder
export function makeGunzip(pako, defaults = {}) {
  return function gunzip(opts = {}) {
    const check = opts.crc ?? defaults.crc ?? true;
    const onWarn = opts.onWarn ?? defaults.onWarn ?? (() => {});
    let state = "header"; // header | body | trailer | done
    let members = 0;
    let totalIn = 0;
    let junk = 0;
    let pending = new Uint8Array(0); // an incomplete header
    let inf = null;
    let crc = 0;
    let len = 0;
    const trailer = new Uint8Array(8);
    let trailerHave = 0;
    let failed = null;

    // pako hands out output only when its 16 KiB buffer is full or the member ends; what is decoded and still in the buffer when
    // the input ends early is valid data (a clean prefix), delivered here before the error is raised
    function emitPartial(ctl) {
      const st = inf && inf.strm;
      if (st && st.output && st.next_out > 0) {
        const part = st.output.subarray(0, st.next_out);
        if (check) crc = crc32(crc, part);
        len += part.length;
        ctl.enqueue(part.slice());
        st.next_out = 0;
        st.avail_out = st.output.length;
      }
    }

    function run(data, ctl) {
      while (data.length > 0) {
        if (state === "done") { junk += data.length; return; }
        if (state === "header") {
          data = concat(pending, data);
          pending = new Uint8Array(0);
          if (members > 0) { // after the first member, anything that is not a gzip header is trailing junk
            if (data[0] !== 0x1f || (data.length > 1 && data[1] !== 0x8b)) { state = "done"; junk += data.length; return; }
          }
          const h = parseHeader(data);
          if (h < 0) {
            if (data.length > MAX_HEADER) throw new Error("gzip header too long: not a gzip file");
            pending = data;
            return;
          }
          data = data.subarray(h);
          inf = new pako.Inflate({ raw: true });
          inf.onData = (c) => {
            if (check) crc = crc32(crc, c);
            len += c.length;
            ctl.enqueue(c); // pako hands over a buffer it does not reuse
          };
          crc = 0;
          len = 0;
          state = "body";
        } else if (state === "body") {
          inf.push(data, false);
          // no emitPartial here: after damage the bytes still in pako's buffer may be garbage decoded before the error showed
          if (inf.err) throw new Error(`gzip data is corrupt: ${inf.msg || inf.err}`);
          if (!inf.ended) return; // all of `data` was consumed, the member goes on in the next chunk
          data = data.subarray(data.length - inf.strm.avail_in);
          inf = null;
          state = "trailer";
          trailerHave = 0;
        } else { // trailer: CRC-32 and ISIZE, both little-endian
          const n = Math.min(8 - trailerHave, data.length);
          trailer.set(data.subarray(0, n), trailerHave);
          trailerHave += n;
          data = data.subarray(n);
          if (trailerHave === 8) {
            const want = (trailer[0] | (trailer[1] << 8) | (trailer[2] << 16) | (trailer[3] << 24)) >>> 0;
            const size = (trailer[4] | (trailer[5] << 8) | (trailer[6] << 16) | (trailer[7] << 24)) >>> 0;
            members++;
            state = "header";
            if (check && (crc !== want || (len >>> 0) !== size)) {
              throw new Error(
                crc !== want
                  ? `gzip checksum mismatch in member ${members}: the file is damaged`
                  : `gzip length mismatch in member ${members}: the file is damaged`,
              );
            }
          }
        }
      }
    }

    return new TransformStream({
      transform(chunk, ctl) {
        if (failed) throw failed;
        totalIn += chunk.length;
        try {
          run(chunk, ctl);
        } catch (e) {
          failed = e; // raised on the next chunk or at the end: what was decoded from this chunk is delivered first
        }
      },
      async flush(ctl) {
        // Erroring the readable side discards the chunks still queued in it, so an error found at the end waits until the reader
        // has taken what was decoded (the queue is empty when desiredSize is no longer negative); a reader that stops reading
        // and cancels leaves this waiting, which costs nothing.
        const fail = async (e) => {
          while (ctl.desiredSize !== null && ctl.desiredSize < 0) await new Promise((r) => setTimeout(r, 0));
          throw e;
        };
        if (failed) await fail(failed);
        if (totalIn === 0) await fail(new Error("empty file"));
        const early = new Error("gzip file ends early (truncated or still downloading)");
        if (state === "body") emitPartial(ctl);
        if (state === "body" || state === "trailer") await fail(early);
        if (state === "header" && pending.length > 0) {
          // part of a header: of the first member, or of one that follows (a file cut in the middle of a member boundary);
          // anything else after the last member is junk
          if (members === 0 || (pending[0] === 0x1f && (pending.length === 1 || pending[1] === 0x8b))) await fail(early);
          junk += pending.length;
        }
        if (junk > 0) {
          onWarn(`${junk} byte${junk === 1 ? "" : "s"} after the last gzip member ${junk === 1 ? "was" : "were"} ignored (padding or junk)`);
        }
      },
    });
  };
}
