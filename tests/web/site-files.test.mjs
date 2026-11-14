// node --test tests/web/site-files.test.mjs : everything a page or worker in web/ imports is a file in web/ that scripts/build_site.sh
// copies. The deployed site once lacked four modules added after the file list in the script was written; nothing ran on Run.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readdirSync, readFileSync, existsSync, statSync } from "node:fs";
import { join, dirname, resolve, relative } from "node:path";

const root = join(import.meta.dirname, "../..");
const web = join(root, "web");
const built = new Set(["pkg", "pkg-ps"]); // directories build_site.sh creates from the wasm builds

const walk = (d) => readdirSync(d).flatMap((n) => (statSync(join(d, n)).isDirectory() ? walk(join(d, n)) : [join(d, n)]));
const sources = walk(web).filter((f) => /\.(js|mjs)$/.test(f));
const imports = (file) => {
  const text = readFileSync(file, "utf8");
  const out = [];
  for (const m of text.matchAll(/(?:from|import)\s*\(?\s*["'](\.[^"']+)["']/g)) out.push(m[1]);
  for (const m of text.matchAll(/new Worker\(\s*["']([^"']+)["']/g)) out.push(m[1].startsWith(".") ? m[1] : `./${m[1]}`);
  for (const m of text.matchAll(/new URL\(\s*["'](\.[^"']+)["']\s*,\s*import\.meta\.url/g)) out.push(m[1]);
  return out;
};

test("every relative import and worker of the web app exists in web/ (outside the wasm packages)", () => {
  const found = new Set();
  for (const f of sources) {
    for (const rel of imports(f)) {
      const target = resolve(dirname(f), rel);
      if (built.has(relative(web, target).split("/")[0])) continue;
      assert.ok(existsSync(target), `${relative(root, f)} imports ${rel}, which is not in web/`);
      found.add(relative(web, target));
    }
  }
  // the pattern must see the imports that matter (a regex that finds nothing would pass)
  for (const must of ["gunzip.js", "gunzip-browser.js", "readhead.js", "archive-worker.js", "vendor/pako_inflate.mjs", "popsketch-worker.js", "indexsource.js"]) {
    assert.ok(found.has(must), `${must} is not found among the imports: ${[...found].join(", ")}`);
  }
});

test("build_site.sh copies every page, script and stylesheet and the vendored libraries by pattern, not by a list of names", () => {
  const sh = readFileSync(join(root, "scripts/build_site.sh"), "utf8");
  assert.match(sh, /web\/\*\.html/);
  assert.match(sh, /web\/\*\.js/);
  assert.match(sh, /web\/\*\.css/);
  assert.match(sh, /cp -r .*web\/vendor/);
  // and the pages and scripts it relies on are the ones with those extensions (no .mjs or .json needed at the top level of web/)
  const top = readdirSync(web).filter((n) => statSync(join(web, n)).isFile());
  for (const n of top) assert.match(n, /\.(html|js|css)$/, `web/${n} would not be copied`);
  // a vendored file that is imported must be in web/vendor
  for (const f of sources) for (const rel of imports(f)) if (rel.includes("vendor/")) assert.ok(existsSync(resolve(dirname(f), rel)), rel);
});
