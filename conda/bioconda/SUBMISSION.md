# Submitting kMate to bioconda

The recipe in [`meta.yaml`](meta.yaml) is ready to submit. It fetches the **PyPI
sdist** (~122 KB), not the GitHub auto-tarball (which archives the whole ~185 MB
repo and would be rejected). Do the two steps below in order.

## Step 1 — publish the sdist to PyPI

`kmate` is available on PyPI (checked 2026-06-18). You need a PyPI account and an
API token (https://pypi.org/manage/account/token/).

```bash
mamba activate kmate
cd /global/scratch/users/tbellg/kmate

# build the clean sdist (already built once at dist/kmate-0.1.0.tar.gz)
python -m build --sdist          # needs `mamba install -c conda-forge python-build twine`

# optional but recommended: check metadata renders
twine check dist/kmate-0.1.0.tar.gz

# upload (will prompt for the token, or set TWINE_USERNAME=__token__ TWINE_PASSWORD=pypi-...)
twine upload dist/kmate-0.1.0.tar.gz
```

After upload, `pip install kmate` works, and the source URL in `meta.yaml`
resolves.

### Confirm the sha256

The `sha256` in `meta.yaml` (`ac20fee8…`) is for the sdist built here. **If you
rebuild the sdist, the hash changes** (gzip timestamps), so use the authoritative
hash from PyPI after uploading:

```bash
curl -sL https://pypi.org/pypi/kmate/0.1.0/json \
  | python -c "import sys,json; d=json.load(sys.stdin); print([f['digests']['sha256'] for f in d['urls'] if f['packagetype']=='sdist'][0])"
```

Paste that value into `meta.yaml`'s `source: sha256:` if it differs.

## Step 2 — open the bioconda PR

```bash
# fork + clone bioconda-recipes (one-time)
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout -b add-kmate

mkdir -p recipes/kmate
cp /global/scratch/users/tbellg/kmate/conda/bioconda/meta.yaml recipes/kmate/meta.yaml

git add recipes/kmate/meta.yaml
git commit -m "Add kmate"
git push -u origin add-kmate
gh pr create --repo bioconda/bioconda-recipes --base master \
  --title "Add kmate" --body "k-mer-based founder-mixture allele-frequency estimation for pool-seq."
```

Bioconda's CI then builds the recipe in a clean container and runs the `test:`
stage (`kmate --help`, `kmate --version`, `import kmate`). A reviewer merges it;
then `conda install -c bioconda kmate` works. Expect a round or two of CI
feedback — the recipe is a standard noarch-python one, so it should be light.

## Notes

- `noarch: python` → one build, no per-platform matrix (the single biggest
  difficulty reducer; kMate has no compiled extensions).
- All run deps already live on the channels: numpy/scipy (conda-forge),
  pysam/jellyfish/samtools (bioconda). Verified they co-resolve.
- The bundled `data/selftest/` fixture ships inside the sdist (via
  `[tool.setuptools.package-data]`), so `kmate selftest` works from a conda
  install too.
- Future releases: bump `version`, rebuild + re-upload the sdist, update the
  sha256, and open a version-bump PR (or let bioconda's autobump bot do it).
