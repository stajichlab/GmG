// A popsketch sketch bundle (scripts/popsketch_bundle.py): a .tar.gz / .tar / .zip with popsketch-bundle.json (format 1: k, scaled, the sample
// ids), an optional labels.tsv, and per sample <id>.json, <id>.gsk and, when present, <id>.half_a.gsk and <id>.half_b.gsk.
// readBundle checks that everything the manifest names is there and returns the bytes for popsketch-wasm's Session.add_sketch.
import { expandArchiveInWorker } from "./archive.js";

export const BUNDLE_FORMAT = 1;

export async function readBundle(file) {
  const { files } = await expandArchiveInWorker(file, "bundle");
  const by = new Map(files.map((f) => [f.name, f]));
  const mf = by.get("popsketch-bundle.json");
  if (!mf) throw new Error(`${file.name}: not a popsketch bundle (no popsketch-bundle.json)`);
  const manifest = JSON.parse(await mf.text());
  if (manifest.format !== BUNDLE_FORMAT) throw new Error(`${file.name}: bundle format ${manifest.format}, this page reads format ${BUNDLE_FORMAT}`);
  if (!Array.isArray(manifest.samples) || !manifest.samples.length) throw new Error(`${file.name}: the bundle lists no samples`);
  const bytes = async (name, need) => {
    const f = by.get(name);
    if (!f) {
      if (need) throw new Error(`${file.name}: ${name} is listed in the bundle but missing`);
      return new Uint8Array(0);
    }
    return new Uint8Array(await f.arrayBuffer());
  };
  const sketches = [];
  for (const id of manifest.samples) {
    sketches.push({ id, meta: await bytes(`${id}.json`, true), whole: await bytes(`${id}.gsk`, true), a: await bytes(`${id}.half_a.gsk`, false), b: await bytes(`${id}.half_b.gsk`, false) });
  }
  const labels = by.get("labels.tsv") ? await by.get("labels.tsv").text() : null;
  return { k: manifest.k, scaled: manifest.scaled, note: manifest.note || "", labels, sketches };
}

// Fetch a bundle from a URL (the server must allow cross-origin requests) as a File that readBundle can open. The kind of archive is read
// from the first bytes, not the URL, so a link without a file extension works.
// The download is limited to `maxBytes` (a real bundle is a few MB to a few hundred MB): a link on the page that is opened
// automatically (`?data=URL`) must not be able to make the tab download without end, which highly compressible data would allow.
export const MAX_BUNDLE_BYTES = 1 * 2 ** 30;
export async function fetchBundle(url, fetchFn = fetch, maxBytes = MAX_BUNDLE_BYTES) {
  const r = await fetchFn(url);
  if (!r.ok) throw new Error(`${url}: ${r.status} ${r.statusText}`);
  const tooBig = () => new Error(`${url}: the bundle is larger than ${(maxBytes / 2 ** 20).toFixed(0)} MiB, which this page does not open`);
  const declared = Number(r.headers?.get?.("content-length"));
  if (Number.isFinite(declared) && declared > maxBytes) { r.body?.cancel?.(); throw tooBig(); }
  let blob;
  if (r.body?.getReader) { // count what arrives, whatever the headers said (content-length can be absent or the compressed size)
    const reader = r.body.getReader();
    const parts = [];
    let n = 0;
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      n += value.length;
      if (n > maxBytes) { await reader.cancel(); throw tooBig(); }
      parts.push(value);
    }
    blob = new Blob(parts);
  } else {
    blob = await r.blob();
    if (blob.size > maxBytes) throw tooBig();
  }
  const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer());
  const ext = head[0] === 0x1f && head[1] === 0x8b ? "tar.gz" : head[0] === 0x50 && head[1] === 0x4b ? "zip" : "tar";
  return new File([blob], `bundle.${ext}`);
}
