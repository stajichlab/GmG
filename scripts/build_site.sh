#!/bin/bash
# usage: build_site.sh OUT_DIR [INDEX.gmgi ...]
# Assembles the static web app: web/*, the wasm package (pixi run wasm -> crates/gmg-wasm/pkg-web) and optional bundled indexes
# with index/manifest.json. Serve OUT_DIR with any static file server (the worker needs http(s), not file://).
set -euo pipefail
out=${1:?output directory}; shift || true
root=$(cd "$(dirname "$0")/.." && pwd)
rm -rf "$out"; mkdir -p "$out/pkg" "$out/index" "$out/examples"
cp "$root"/web/index.html "$root"/web/style.css "$root"/web/app.js "$root"/web/worker.js "$out/"
cp "$root"/examples/examples.tsv "$out/examples/"
cp "$root"/crates/gmg-wasm/pkg-web/gmg_wasm.js "$root"/crates/gmg-wasm/pkg-web/gmg_wasm_bg.wasm "$out/pkg/"
python3 - "$out" "$@" <<'PY'
import json, os, re, shutil, sys
out, files = sys.argv[1], sys.argv[2:]
man = []
for f in files:
    name = os.path.basename(f)
    shutil.copy(f, os.path.join(out, "index", name))
    mb = round(os.path.getsize(f) / 1e6)
    title = "Representatives (fast, fewer species calls)" if "reps" in name else "Full reference"
    v = re.search(r"(\d{4})\.(\d{2})(?:\.(\d+))?", name)  # YYYY.MM[.N]: the monthly freeze and the build of it
    if v:
        title += f", {v[1]}.{v[2]}" + (f" build {v[3]}" if v[3] else "")
    man.append({"file": name, "title": title, "size_mb": mb})
json.dump(man, open(os.path.join(out, "index", "manifest.json"), "w"), indent=1)
print("site:", out, [m["file"] for m in man])
PY
