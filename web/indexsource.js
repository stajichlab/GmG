// The reference indexes bundled with the site (index/manifest.json, written by scripts/build_site.sh) and how to fetch one for popsketch.
// fetch and the Cache API are parameters so that node can test it.

export const INDEX_CACHE = "popsketch-index-v1";
export const MAX_INDEX_BYTES = 1 * 2 ** 30; // an index larger than this is not fetched, whatever the manifest says

/** The manifest entries [{file, title, size_mb, bytes}], or [] when there is none (a site without bundled indexes, or offline). */
export async function loadManifest(fetchFn = fetch, url = "index/manifest.json") {
  try {
    const r = await fetchFn(url, { cache: "no-cache" });
    if (!r.ok) return [];
    const list = await r.json();
    return Array.isArray(list) ? list.filter((e) => e && typeof e.file === "string") : [];
  } catch {
    return [];
  }
}

/**
 * The bytes of a bundled index. A copy in the Cache API is used when its size is the manifest's (a stale or half-written copy is
 * dropped), otherwise the file is downloaded with `onProgress(loaded, total)`, checked against the manifest size, and cached only
 * when complete. Returns { bytes, cached }.
 * The size from the manifest is the one to trust: content-length counts the compressed bytes when the server gzips the file.
 */
export async function fetchIndexBytes(entry, { fetchFn = fetch, cacheStorage = globalThis.caches, onProgress = () => {}, base = "index/" } = {}) {
  const url = `${base}${entry.file}`;
  const want = Number(entry.bytes) || 0;
  if (want > MAX_INDEX_BYTES) throw new Error(`${entry.file}: ${(want / 2 ** 20).toFixed(0)} MiB is more than this page will fetch`);
  let cache = null;
  try { cache = cacheStorage ? await cacheStorage.open(INDEX_CACHE) : null; } catch { cache = null; /* private mode, no storage */ }
  if (cache) {
    const hit = await cache.match(url).catch(() => undefined);
    if (hit) {
      const bytes = new Uint8Array(await hit.arrayBuffer());
      if (!want || bytes.length === want) { onProgress(bytes.length, bytes.length); return { bytes, cached: true }; }
      await cache.delete(url).catch(() => {}); // a copy of another size: stale or damaged
    }
  }
  const r = await fetchFn(url);
  if (!r.ok) throw new Error(`could not fetch ${entry.file} (${r.status})`);
  const limit = want || MAX_INDEX_BYTES;
  const chunks = [];
  let loaded = 0;
  if (r.body?.getReader) {
    const reader = r.body.getReader();
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      loaded += value.length;
      if (loaded > limit * 1.001) { await reader.cancel(); throw new Error(`${entry.file}: more bytes than the ${(limit / 2 ** 20).toFixed(0)} MiB expected`); }
      chunks.push(value);
      onProgress(loaded, want);
    }
  } else {
    const b = new Uint8Array(await r.arrayBuffer());
    chunks.push(b);
    loaded = b.length;
  }
  if (want && loaded !== want) throw new Error(`${entry.file}: the download is incomplete or changed (${loaded} of ${want} bytes); try again`);
  const bytes = new Uint8Array(loaded);
  let off = 0;
  for (const c of chunks) { bytes.set(c, off); off += c.length; }
  if (cache) await cache.put(url, new Response(bytes.slice(), { headers: { "content-type": "application/octet-stream" } })).catch(() => {});
  return { bytes, cached: false };
}

/** True when the Cache API can hold the download (then the page need not keep the bytes in memory between runs). */
export async function cacheWorks(cacheStorage = globalThis.caches) {
  try { return !!(cacheStorage && (await cacheStorage.open(INDEX_CACHE))); } catch { return false; }
}
