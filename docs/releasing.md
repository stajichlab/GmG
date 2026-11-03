# Releasing GmG

Code releases are tags `vX.Y.Z`; reference indexes are separate releases (`index-YYYY.MM[.N]`, see [index-distribution.md](index-distribution.md)).

## A code release

1. A PR "Prepare vX.Y.Z": the version in `Cargo.toml`, `crates/gmg-py/Cargo.toml` and `crates/gmg-py/pyproject.toml`; the lockfiles (`cargo build`, and `cargo metadata --manifest-path crates/gmg-py/Cargo.toml`); `CHANGELOG.md`
   (a section `## vX.Y.Z (date)`, replacing `## Unreleased`). The three versions must be equal: `pypi.yml` refuses a tag that differs from `pyproject.toml`.
2. Wait for every check on the PR, then merge it.
3. From `main`: `git tag vX.Y.Z && git push origin vX.Y.Z`. This starts:
   * **Release** (`release.yml`): binaries (linux-64, osx-arm64), the wasm package, the image `ghcr.io/stajichlab/gmg:X.Y.Z` and `:latest`, and the GitHub release with its checksums.
     "Run workflow" by hand is a dry run that publishes nothing.
   * **PyPI** (`pypi.yml`): wheels for Linux x86_64, Linux aarch64 and macOS arm64 (smoke-tested in a clean venv), then the upload. The upload job waits in the `pypi` environment
     until a required reviewer approves it: Actions → the run → **Review deployments** → tick `pypi` → **Approve and deploy**.
4. Check: the GitHub release has the six assets; <https://pypi.org/project/guess-my-genome/> lists the version with three wheels; a clean `pip install guess-my-genome` imports `gmg` at the new version.
5. If a tagged version needs a fix, make a new patch version. Neither PyPI nor the tag history reuses a version number.

The Pages site (`pages.yml`) and the `:edge` image follow `main` and are not part of a release.

## PyPI: one-time setup (already done for `guess-my-genome`)

Trusted publishing: no API token is stored anywhere.

1. A PyPI account with two-factor authentication.
2. pypi.org → Account settings → Publishing → **Add a new pending publisher**: project `guess-my-genome`, owner `stajichlab`, repository `GmG`, workflow `pypi.yml`, environment `pypi`.
   The first upload creates the project. The workflow file name and the environment name must match exactly.
3. GitHub → Settings → Environments → `pypi` (add yourself as a required reviewer so every publish needs a click).
4. Optional rehearsal on TestPyPI (a separate account and publisher at test.pypi.org with environment `testpypi`): Actions → PyPI → Run workflow → tick "also publish to TestPyPI".

Wheels only: there is no source distribution, because the build reads `data/thresholds/*.json` from outside the crate and an sdist built by maturin does not contain them. Other platforms get "no matching distribution".
To add Intel Macs or Windows, add a row to the `wheels` matrix in `pypi.yml`.
The name `gmg` on PyPI belongs to an unrelated project, hence the distribution name `guess-my-genome` (the import name stays `gmg`).
