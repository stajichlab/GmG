// The multi-member gunzip (gunzip.js) with pako (vendor/pako_inflate.mjs, the unmodified pako 1.0.11 inflate build wrapped as
// an ES module so that it behaves the same in a page, a worker and node), and a reader that decides by the first bytes.
import pako from "./vendor/pako_inflate.mjs";
import { makeGunzip } from "./gunzip.js";

export const gunzip = makeGunzip(pako);

/**
 * The bytes of `file`, decoded when it is gzip. Gzip is recognised by its first two bytes (1f 8b), not by the file name, so a
 * renamed file works and a plain file called .gz does not fail. Two bytes are collected before deciding, since a chunk can hold
 * one. `opts` as for gunzip(): { crc, onWarn }; a plain file is passed through unchanged.
 */
export async function openDecoded(file, opts = {}) {
  const reader = file.stream().getReader();
  let head = new Uint8Array(0);
  while (head.length < 2) {
    const { done, value } = await reader.read();
    if (done) break;
    const n = new Uint8Array(head.length + value.length);
    n.set(head);
    n.set(value, head.length);
    head = n;
  }
  const isGz = head.length >= 2 && head[0] === 0x1f && head[1] === 0x8b;
  let first = head.length > 0 ? head : null;
  const raw = new ReadableStream({
    async pull(ctl) {
      if (first) { ctl.enqueue(first); first = null; return; }
      const { done, value } = await reader.read();
      if (done) ctl.close();
      else ctl.enqueue(value);
    },
    cancel(reason) { return reader.cancel(reason); },
  });
  return isGz ? raw.pipeThrough(gunzip(opts)) : raw;
}
