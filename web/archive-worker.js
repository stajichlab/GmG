// Worker for expandArchiveInWorker (archive.js): unpacks one archive off the page's thread.
// In: {file: File, keep: "seq" | "bundle"}   Out: {ok: true, files: [{file, innerPath}], skipped} or {ok: false, message}
import { expandArchive } from "./archive.js";

self.onmessage = async (e) => {
  try {
    const r = await expandArchive(e.data.file, e.data.keep);
    self.postMessage({ ok: true, files: r.files.map((f) => ({ file: f, innerPath: f.innerPath })), skipped: r.skipped });
  } catch (err) {
    self.postMessage({ ok: false, message: String((err && err.message) || err) });
  }
};
