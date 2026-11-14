// node --test tests/web/indexsource.test.mjs : the bundled reference indexes of the popsketch page (web/indexsource.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import { loadManifest, fetchIndexBytes, cacheWorks, INDEX_CACHE, MAX_INDEX_BYTES } from "../../web/indexsource.js";

const bytesOf = (n, v = 7) => new Uint8Array(n).fill(v);
const entry = (n, file = "bfd-reps-test.gmgi") => ({ file, title: "t", size_mb: 1, bytes: n });

// a Cache API in memory
function fakeCaches({ failPut = false, failOpen = false } = {}) {
  const store = new Map();
  const log = [];
  return {
    store, log,
    async open(name) {
      if (failOpen) throw new Error("no storage");
      assert.equal(name, INDEX_CACHE);
      return {
        async match(u) { const b = store.get(u); return b ? new Response(b.slice()) : undefined; },
        async put(u, r) { if (failPut) throw new Error("quota"); log.push(`put ${u}`); store.set(u, new Uint8Array(await r.arrayBuffer())); },
        async delete(u) { log.push(`delete ${u}`); return store.delete(u); },
      };
    },
  };
}
// a fetch that serves `data` in chunks and counts the requests
function fakeFetch(data, { status = 200, chunk = 1000 } = {}) {
  const f = async () => {
    f.calls++;
    if (status !== 200) return new Response("no", { status });
    return new Response(new ReadableStream({ start(c) { for (let i = 0; i < data.length; i += chunk) c.enqueue(data.subarray(i, i + chunk)); c.close(); } }));
  };
  f.calls = 0;
  return f;
}

test("loadManifest returns the entries, and [] when there is none, it is not JSON or it is not a list", async () => {
  const list = [{ file: "a.gmgi", title: "A", size_mb: 1, bytes: 1 }, { nofile: true }, null];
  assert.deepEqual((await loadManifest(async () => new Response(JSON.stringify(list)))).map((e) => e.file), ["a.gmgi"]);
  assert.deepEqual(await loadManifest(async () => new Response("no", { status: 404 })), []);
  assert.deepEqual(await loadManifest(async () => new Response("<html>")), []);
  assert.deepEqual(await loadManifest(async () => new Response("{}")), []);
  assert.deepEqual(await loadManifest(async () => { throw new TypeError("offline"); }), []);
});

test("a download reports progress, is checked against the manifest size and is cached; the second time comes from the cache", async () => {
  const data = bytesOf(4500), caches = fakeCaches(), f = fakeFetch(data), seen = [];
  const a = await fetchIndexBytes(entry(4500), { fetchFn: f, cacheStorage: caches, onProgress: (l, t) => seen.push([l, t]) });
  assert.equal(a.cached, false);
  assert.deepEqual(a.bytes, data);
  assert.deepEqual(seen.at(-1), [4500, 4500]);
  assert.ok(seen.length >= 4 && seen.every(([l], i) => i === 0 || l > seen[i - 1][0]), "progress rises");
  assert.deepEqual(caches.log, ["put index/bfd-reps-test.gmgi"]);
  const b = await fetchIndexBytes(entry(4500), { fetchFn: f, cacheStorage: caches });
  assert.equal(b.cached, true);
  assert.deepEqual(b.bytes, data);
  assert.equal(f.calls, 1, "not fetched again");
});

test("a cached copy of another size is dropped and fetched again", async () => {
  const caches = fakeCaches();
  caches.store.set("index/bfd-reps-test.gmgi", bytesOf(100, 1)); // stale or half-written
  const f = fakeFetch(bytesOf(2000, 9));
  const r = await fetchIndexBytes(entry(2000), { fetchFn: f, cacheStorage: caches });
  assert.equal(r.cached, false);
  assert.equal(r.bytes[0], 9);
  assert.equal(f.calls, 1);
  assert.ok(caches.log.includes("delete index/bfd-reps-test.gmgi"));
  assert.equal(caches.store.get("index/bfd-reps-test.gmgi").length, 2000);
});

test("an incomplete or changed download is an error and is not cached", async () => {
  const caches = fakeCaches();
  await assert.rejects(fetchIndexBytes(entry(5000), { fetchFn: fakeFetch(bytesOf(3000)), cacheStorage: caches }), /incomplete or changed \(3000 of 5000/);
  assert.equal(caches.store.size, 0);
  await assert.rejects(fetchIndexBytes(entry(1000), { fetchFn: fakeFetch(bytesOf(1000), { status: 404 }), cacheStorage: caches }), /404/);
});

test("a download that goes on past the expected size is cut, an index the manifest says is too large is not fetched", async () => {
  let read = 0;
  const endless = async () => new Response(new ReadableStream({ pull(c) { read += 4096; c.enqueue(new Uint8Array(4096)); } }));
  await assert.rejects(fetchIndexBytes(entry(20000), { fetchFn: endless, cacheStorage: fakeCaches() }), /more bytes than/);
  assert.ok(read < 20000 + 3 * 4096, `stopped after ${read} bytes`);
  const f = fakeFetch(bytesOf(10));
  await assert.rejects(fetchIndexBytes(entry(MAX_INDEX_BYTES + 1), { fetchFn: f, cacheStorage: fakeCaches() }), /more than this page will fetch/);
  assert.equal(f.calls, 0);
});

test("without a usable cache (private mode, no Cache API, a failing put) the download still works", async () => {
  const data = bytesOf(3000, 3);
  for (const cacheStorage of [undefined, fakeCaches({ failOpen: true }), fakeCaches({ failPut: true })]) {
    const r = await fetchIndexBytes(entry(3000), { fetchFn: fakeFetch(data), cacheStorage });
    assert.deepEqual(r.bytes, data);
  }
  assert.equal(await cacheWorks(undefined), false);
  assert.equal(await cacheWorks(fakeCaches({ failOpen: true })), false);
  assert.equal(await cacheWorks(fakeCaches()), true);
});

test("a manifest entry without a size is still downloaded, within the global limit, and not size-checked", async () => {
  const r = await fetchIndexBytes({ file: "x.gmgi" }, { fetchFn: fakeFetch(bytesOf(1234)), cacheStorage: fakeCaches() });
  assert.equal(r.bytes.length, 1234);
});
