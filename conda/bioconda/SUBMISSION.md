# Submitting kMate to bioconda

> **STATUS (2026-09-29): ready to submit.** kMate 0.1.2 is on PyPI and the recipe below
> is verified against it: the `sha256` matches what PyPI serves, and installing that sdist
> into a pristine conda env holding only the recipe's declared run deps gives
> `kmate selftest` PASS (R^2=0.996, RMSE=0.0032).
>
> bioconda currently serves **0.1.0**, which cannot count k-mers (wrong `jellyfish`
> dependency) and predates the per-founder normalization fix. This is an **update** to
> `recipes/kmate/meta.yaml`, not a new recipe.

The recipe in [`meta.yaml`](meta.yaml) is ready to submit. It fetches the **PyPI
sdist** (~122 KB), not the GitHub auto-tarball (which archives the whole ~185 MB
repo and would be rejected). Do the two steps below in order.

## Step 1 — publish the sdist to PyPI

0.1.2 is **already uploaded**. The `sha256` in [`meta.yaml`](meta.yaml) is the hash PyPI
serves, confirmed by download:

```
1a509f24e71026eb0cc99350bbfc245f98bf7d83f73650f9d543e6f4934f5d13
```

> ⚠️ **Upload that artifact, do not rebuild it.** Rebuilding produces a different
> hash (gzip embeds a timestamp), which would make `meta.yaml` wrong and fail
> bioconda CI. If you must rebuild, re-read the hash from PyPI afterwards using the
> command below and paste it into `meta.yaml`.

```bash
mamba activate kmate
cd <kMate checkout>

python -m twine check dist/kmate-0.1.2.tar.gz     # already PASSES
python -m twine upload dist/kmate-0.1.2.tar.gz    # needs a PyPI API token
```

After upload, `pip install kmate==0.1.2` works and the `url:` in `meta.yaml` resolves.

### Confirm the hash PyPI actually serves

```bash
curl -sL https://pypi.org/pypi/kmate/0.1.2/json \
  | python -c "import sys,json; d=json.load(sys.stdin); print([f['digests']['sha256'] for f in d['urls'] if f['packagetype']=='sdist'][0])"
```

It must equal the value above. If not, paste what PyPI reports into `meta.yaml`.

### What 0.1.2 fixes (say this in the PR)

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

## Step 2 — open the bioconda PR (this is an UPDATE, not a new recipe)

`recipes/kmate/meta.yaml` **already exists** in bioconda-recipes at version 0.1.0 —
that is the broken build. So this is a version bump to an existing recipe, not a new
submission.

```bash
# fork + clone bioconda-recipes (one-time)
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout master && git pull upstream master
git checkout -b kmate-0.1.2

# overwrite the existing recipe with ours
cp <kMate checkout>/conda/bioconda/meta.yaml recipes/kmate/meta.yaml

git add recipes/kmate/meta.yaml
git commit -m "Update kmate to 0.1.1"
git push -u origin kmate-0.1.2

gh pr create --repo bioconda/bioconda-recipes --base master \
  --title "Update kmate to 0.1.1" \
  --body "kMate 0.1.1.

Fixes two independent problems in the published 0.1.0, both reproduced against a
clean \\`mamba create -c conda-forge -c bioconda kmate\\`:

1. **The package cannot count k-mers.** The recipe depended on \\`jellyfish\\`, which on
   conda-forge is a Python string-similarity library and ships no \\`jellyfish\\`
   binary. kMate shells out to the k-mer counter, so \\`kmate selftest\\` exited 1 with
   \\`jellyfish: command not found\\`. Now depends on \\`kmer-jellyfish\\`.
2. **0.1.0 was built from a stale snapshot**, predating an estimator fix; it lacks
   \\`--normalize\\`, \\`--unit\\` and \\`--emit-af-se\\`.

The \\`test:\\` block now also runs \\`kmate selftest\\`, a bundled offline end-to-end
fixture that exercises jellyfish and samtools. The previous test ran only
\\`--help\\`/\\`--version\\`/\\`import\\`, which is why a package that could not count
k-mers passed CI. Verified locally: selftest passes from the 0.1.2 sdist in a clean
environment."
```

> **Order matters:** upload to PyPI (Step 1) *before* opening the PR. The recipe
> fetches the sdist from PyPI, so CI fails on a 404 if the release is not up yet.

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
