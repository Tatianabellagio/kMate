# Submitting kMate to bioconda

The recipe in [`meta.yaml`](meta.yaml) is ready to submit. It fetches the **PyPI
sdist** (~122 KB), not the GitHub auto-tarball (which archives the whole ~185 MB
repo and would be rejected). Do the two steps below in order.

## Step 1 — publish the sdist to PyPI

The 0.1.1 sdist is **already built** at `dist/kmate-0.1.1.tar.gz`, and the `sha256`
in [`meta.yaml`](meta.yaml) is the hash of *that exact file*:

```
11ec193178131d2c51b9790430e0871d377b617cedc3bc203922b9617a499752
```

> ⚠️ **Upload that artifact, do not rebuild it.** Rebuilding produces a different
> hash (gzip embeds a timestamp), which would make `meta.yaml` wrong and fail
> bioconda CI. If you must rebuild, re-read the hash from PyPI afterwards using the
> command below and paste it into `meta.yaml`.

```bash
mamba activate kmate
cd <kMate checkout>

python -m twine check dist/kmate-0.1.1.tar.gz     # already PASSES
python -m twine upload dist/kmate-0.1.1.tar.gz    # needs a PyPI API token
```

After upload, `pip install kmate==0.1.1` works and the `url:` in `meta.yaml` resolves.

### Confirm the hash PyPI actually serves

```bash
curl -sL https://pypi.org/pypi/kmate/0.1.1/json \
  | python -c "import sys,json; d=json.load(sys.stdin); print([f['digests']['sha256'] for f in d['urls'] if f['packagetype']=='sdist'][0])"
```

It must equal the value above. If not, paste what PyPI reports into `meta.yaml`.

### What 0.1.1 fixes (say this in the PR)

The published 0.1.0 is broken in two independent ways, both verified against a clean
`mamba create -c conda-forge -c bioconda kmate`:

1. **It cannot count k-mers.** The recipe depended on `jellyfish`, which on
   conda-forge is a Python string-similarity library shipping no `jellyfish`
   binary. `kmate selftest` exited 1 with `jellyfish: command not found`. Fixed by
   depending on **`kmer-jellyfish`** (the real counter, 2.3.1).
2. **It is a pre-July-2026 snapshot.** It has no `--normalize`, `--unit` or
   `--emit-af-se`, and `per_founder` appears zero times in its source — so it would
   silently run the older EM in which k-mer-poor founders collapse toward zero.

The recipe's `test:` block now also runs `kmate selftest`, which exercises
jellyfish/samtools end-to-end. The old test ran only `--help`, `--version` and
`import kmate`, which is exactly why a package that could not count k-mers passed CI.

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
