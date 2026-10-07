# Releasing kMate to PyPI and bioconda

> **STATUS (2026-10-06): 0.1.4.** The recipe in [`meta.yaml`](meta.yaml) pins the sha256
> of the 0.1.4 sdist uploaded to PyPI. History: 0.1.0 (June) could not count k-mers
> (`jellyfish` instead of `kmer-jellyfish`); 0.1.2 (2026-09-29) fixed that; 0.1.3 bundles
> the decomposition scripts, adds `decompose --gfa` and `bcftools`, and fixes
> `--haploidize` on haploid input; 0.1.4 caps Python below 3.13 (jellyfish bindings),
> adds `build-index --tmp-dir`, and reports low-evidence fits instead of returning NaN.

bioconda's `recipes/kmate/meta.yaml` already exists, so each release is a **version
bump** of that recipe. The recipe fetches the **PyPI sdist** (~170 KB), not the GitHub
auto-tarball, which archives the whole repo and would be rejected.

## 1. Build once, check, then pin the hash

Bump `version` in `pyproject.toml`, `src/kmate/__init__.py`, `conda/meta.yaml` and
`meta.yaml` here. Everything that goes into the sdist (`src/`, `README.md`,
`pyproject.toml`, `tests/test_*.py`) must be final before building; `conda/` is not in
the sdist, so the hash can be written into `meta.yaml` afterwards.

```bash
mamba activate kmate
python -m build --sdist --wheel --no-isolation -o dist/
python -m twine check dist/kmate-X.Y.Z*
sha256sum dist/kmate-X.Y.Z.tar.gz          # -> meta.yaml `sha256:`
```

> ⚠️ **Upload the artifact you hashed; do not rebuild it.** gzip embeds a timestamp, so
> a rebuild changes the hash and bioconda CI fails.

Before uploading, install the sdist into a pristine env holding only the recipe's run
deps and run the offline checks:

```bash
mamba create -p /tmp/kmate_rel -c conda-forge -c bioconda python numpy scipy pysam \
    kmer-jellyfish samtools bcftools pip                        # newest Python, no pins
/tmp/kmate_rel/bin/pip install --no-deps dist/kmate-X.Y.Z.tar.gz
PATH=/tmp/kmate_rel/bin:$PATH kmate selftest                  # must PASS
PATH=/tmp/kmate_rel/bin:$PATH python tests/test_haploidize.py
PATH=/tmp/kmate_rel/bin:$PATH python tests/test_cli_flags.py
PATH=/tmp/kmate_rel/bin:$PATH python tests/test_low_evidence.py
PATH=/tmp/kmate_rel/bin:$PATH python tests/test_jellyfish_counter.py
```

For a release that changes the panel builders, also run `tests/e2e/` (fresh install,
Chr1 panel compared with the production matrices; ~4 h).

## 2. Upload to PyPI, then confirm the served hash

```bash
python -m twine upload dist/kmate-X.Y.Z.tar.gz dist/kmate-X.Y.Z-py3-none-any.whl
curl -sL https://pypi.org/pypi/kmate/X.Y.Z/json \
  | python -c "import sys,json; d=json.load(sys.stdin); print([f['digests']['sha256'] for f in d['urls'] if f['packagetype']=='sdist'][0])"
```

The printed hash must equal `meta.yaml`'s.

## 3. Open the bioconda version-bump PR

Upload to PyPI first: the recipe fetches the sdist from PyPI, so CI fails on a 404 if
the release is not up yet.

```bash
# fork + clone bioconda-recipes (one-time)
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout master && git pull upstream master
git checkout -b kmate-X.Y.Z
cp <kMate checkout>/conda/bioconda/meta.yaml recipes/kmate/meta.yaml
git add recipes/kmate/meta.yaml
git commit -m "Update kmate to X.Y.Z"
git push -u origin kmate-X.Y.Z
gh pr create --repo bioconda/bioconda-recipes --base master --title "Update kmate to X.Y.Z"
```

bioconda's CI builds the recipe in a clean container and runs its `test:` block,
including `kmate selftest`. A reviewer merges it; then
`mamba install -c bioconda kmate=X.Y.Z` works. bioconda's autobump bot may open the PR
itself once the PyPI release is up; close one of the two if both appear.

## 4. After bioconda merges: check the bare install

The solver prefers the newest Python and will quietly fall back to an older kMate whose
recipe allows it. After the merge, the plain command must give the new version:

```bash
mamba create --dry-run -n x -c conda-forge -c bioconda kmate    # must list kmate X.Y.Z
```

(0.1.4 failed this: it capped Python below 3.13, so a bare install picked 0.1.2 on 3.14.)

## Notes

- `noarch: python` → one build, no per-platform matrix (kMate has no compiled code).
- Run deps all live on conda-forge/bioconda: numpy, scipy, pysam, kmer-jellyfish
  (not `jellyfish`, an unrelated Python library), samtools, bcftools.
- The selftest fixture and the vendored PanGenie scripts ship inside the sdist via
  `[tool.setuptools.package-data]`, so both work from a conda install.
