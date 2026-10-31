// Unpack .zip, .tar, .tar.gz and .tgz in the browser (no library: tar is 512-byte headers, zip a central directory, and
// the browser inflates natively with DecompressionStream). Entries come back as File objects held in memory, named by their
// path inside the archive; `keep(path)` picks the ones wanted. Folders, links and other entry types are ignored.

export const ARCHIVE = /\.(zip|tar|tar\.gz|tgz)$/i;
export const MAX_UNPACKED = 3e9; // bytes of wanted entries held at once; the whole archive is read, only wanted entries are kept

const dec = new TextDecoder();
const cut = (b) => { const i = b.indexOf(0); return dec.decode(i < 0 ? b : b.subarray(0, i)); };

// ---------------------------------------------------------------- tar
// ustar and GNU/pax long names; sizes in octal or base-256 (>8 GB).
function tarNumber(b) {
  if (b[0] & 0x80) { let n = b[0] & 0x7f; for (let i = 1; i < b.length; i++) n = n * 256 + b[i]; return n; }
  return parseInt(cut(b).trim() || "0", 8);
}

async function untar(stream, keep, out, skipped) {
  const reader = stream.getReader();
  let buf = new Uint8Array(0), done = false;
  const fill = async () => {
    const r = await reader.read();
    if (r.done) { done = true; return; }
    const n = new Uint8Array(buf.length + r.value.length); n.set(buf); n.set(r.value, buf.length); buf = n;
  };
  const take = async (n) => { while (buf.length < n && !done) await fill(); const k = Math.min(n, buf.length); const o = buf.subarray(0, k); buf = buf.subarray(k); return o; };
  // read `size` bytes as chunks (so a large entry is not copied through the growing buffer), then the padding to 512
  const body = async (size, wanted) => {
    const chunks = []; let left = size;
    while (left > 0) {
      if (!buf.length && !done) await fill();
      if (!buf.length) throw new Error("the archive ends inside a file");
      const k = Math.min(left, buf.length);
      if (wanted) chunks.push(buf.slice(0, k));
      buf = buf.subarray(k); left -= k;
    }
    await take((512 - (size % 512)) % 512);
    return chunks;
  };
  let longName = null, paxPath = null, held = 0;
  for (;;) {
    const h = await take(512);
    if (h.length === 0 || h.every((x) => x === 0)) break;
    if (h.length < 512) throw new Error("the archive ends inside a file");
    const size = tarNumber(h.subarray(124, 136)), type = String.fromCharCode(h[156] || 48);
    let name = cut(h.subarray(0, 100));
    if (cut(h.subarray(257, 262)) === "ustar") { const p = cut(h.subarray(345, 500)); if (p) name = `${p}/${name}`; }
    if (type === "L") { longName = cut(b2(await body(size, true))); continue; }
    if (type === "x") {
      const m = /(?:^|\n)\d+ path=([^\n]*)\n/.exec(dec.decode(b2(await body(size, true))));
      if (m) paxPath = m[1];
      continue;
    }
    if (type === "g") { await body(size, false); continue; }
    if (longName) { name = longName; longName = null; }
    if (paxPath) { name = paxPath; paxPath = null; }
    name = name.replace(/^\.\//, "");
    const regular = type === "0" || type === "7";
    const wanted = regular && keep(name);
    if (regular && !wanted) skipped.notSeq++;
    if (wanted && (held += size) > MAX_UNPACKED) throw new Error(`more than ${MAX_UNPACKED / 1e9} GB of sequence files inside the archive: split it or unpack it on your computer`);
    const chunks = await body(size, wanted);
    if (wanted) out.push(new File(chunks, name));
  }
}
const b2 = (chunks) => { let n = 0; for (const c of chunks) n += c.length; const o = new Uint8Array(n); let i = 0; for (const c of chunks) { o.set(c, i); i += c.length; } return o; };

// ---------------------------------------------------------------- zip
async function unzip(file, keep, out, skipped) {
  const tail = new Uint8Array(await file.slice(Math.max(0, file.size - 65557)).arrayBuffer());
  const tv = new DataView(tail.buffer);
  let e = tail.length - 22;
  while (e >= 0 && tv.getUint32(e, true) !== 0x06054b50) e--;
  if (e < 0) throw new Error("not a zip file (no end-of-archive record)");
  let count = tv.getUint16(e + 10, true), cdSize = tv.getUint32(e + 12, true), cdOff = tv.getUint32(e + 16, true);
  if (cdOff === 0xffffffff || count === 0xffff) { // zip64: the locator sits just before the end record
    const z = new DataView(await file.slice(file.size - tail.length + e - 20, file.size - tail.length + e).arrayBuffer());
    if (z.getUint32(0, true) !== 0x07064b50) throw new Error("zip64 locator missing");
    const r = new DataView(await file.slice(Number(z.getBigUint64(8, true)), Number(z.getBigUint64(8, true)) + 56).arrayBuffer());
    count = Number(r.getBigUint64(32, true)); cdSize = Number(r.getBigUint64(40, true)); cdOff = Number(r.getBigUint64(48, true));
  }
  const cd = new Uint8Array(await file.slice(cdOff, cdOff + cdSize).arrayBuffer()), v = new DataView(cd.buffer);
  let p = 0, held = 0;
  for (let i = 0; i < count && p + 46 <= cd.length; i++) {
    if (v.getUint32(p, true) !== 0x02014b50) throw new Error("damaged zip directory");
    const flags = v.getUint16(p + 8, true), method = v.getUint16(p + 10, true);
    let csize = v.getUint32(p + 20, true), usize = v.getUint32(p + 24, true), off = v.getUint32(p + 42, true);
    const nl = v.getUint16(p + 28, true), xl = v.getUint16(p + 30, true), cl = v.getUint16(p + 32, true);
    const name = dec.decode(cd.subarray(p + 46, p + 46 + nl));
    if (csize === 0xffffffff || usize === 0xffffffff || off === 0xffffffff) { // zip64 extra field: only the saturated fields, in this order
      const x = new DataView(cd.buffer, cd.byteOffset + p + 46 + nl, xl);
      for (let q = 0; q + 4 <= xl;) {
        const id = x.getUint16(q, true), len = x.getUint16(q + 2, true);
        if (id === 1) {
          let r = q + 4;
          if (usize === 0xffffffff) { usize = Number(x.getBigUint64(r, true)); r += 8; }
          if (csize === 0xffffffff) { csize = Number(x.getBigUint64(r, true)); r += 8; }
          if (off === 0xffffffff) { off = Number(x.getBigUint64(r, true)); }
        }
        q += 4 + len;
      }
    }
    p += 46 + nl + xl + cl;
    if (name.endsWith("/")) continue;
    if (!keep(name)) { skipped.notSeq++; continue; }
    if (flags & 1) { skipped.unsupported.push(`${name} (encrypted)`); continue; }
    if (method !== 0 && method !== 8) { skipped.unsupported.push(`${name} (compression method ${method})`); continue; }
    if ((held += usize) > MAX_UNPACKED) throw new Error(`more than ${MAX_UNPACKED / 1e9} GB of sequence files inside the archive: split it or unpack it on your computer`);
    const lh = new DataView(await file.slice(off, off + 30).arrayBuffer());
    if (lh.getUint32(0, true) !== 0x04034b50) throw new Error(`damaged zip entry ${name}`);
    const start = off + 30 + lh.getUint16(26, true) + lh.getUint16(28, true);
    const raw = file.slice(start, start + csize);
    const blob = method === 0 ? raw : await new Response(raw.stream().pipeThrough(new DecompressionStream("deflate-raw"))).blob();
    out.push(new File([blob], name));
  }
}

// ---------------------------------------------------------------- entry point
/** Unpack `file`; resolves to { files: File[], skipped: { notSeq, unsupported: string[] } }. Throws a readable Error on a damaged archive. */
export async function expandArchive(file, keep) {
  const out = [], skipped = { notSeq: 0, unsupported: [] };
  try {
    if (/\.zip$/i.test(file.name)) await unzip(file, keep, out, skipped);
    else if (/\.(tar\.gz|tgz)$/i.test(file.name)) await untar(file.stream().pipeThrough(new DecompressionStream("gzip")), keep, out, skipped);
    else await untar(file.stream(), keep, out, skipped);
  } catch (err) {
    throw new Error(`${file.name}: ${err.message}`);
  }
  // the worker names a result after the file, so entries called genome.fa in different folders get their last two path parts;
  // the full path inside the archive stays in `innerPath` for the file list
  const base = out.map((f) => f.name.split("/").pop());
  return {
    files: out.map((f, i) => {
      const dup = base.filter((b) => b === base[i]).length > 1;
      const n = dup ? f.name.split("/").filter(Boolean).slice(-2).join("_") : base[i];
      const g = new File([f], n);
      Object.defineProperty(g, "innerPath", { value: f.name });
      return g;
    }),
    skipped,
  };
}
