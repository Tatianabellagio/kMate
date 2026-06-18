# Conda packaging

`meta.yaml` is a **noarch: python** recipe (pure Python → one build, no platform
matrix). It builds from the working tree today; for bioconda it needs only a
versioned release tarball + sha256 (see the commented `url:`/`sha256:` lines).

## Build locally

```bash
mamba create -n conda-build -c conda-forge conda-build      # one-time
conda activate conda-build
conda build conda/ -c conda-forge -c bioconda --no-anaconda-upload
```

On a clean machine this builds the package, then runs the `test:` stage
(`kmate --help`, `kmate --version`, `import kmate`).

### On Savio / networked-home clusters

`conda build`'s sqlite channel-index and repodata caches use file `flock`, which
fails on NFS home dirs (`sqlite3.OperationalError: database is locked` /
`LockError`). Redirect every cache to node-local disk and give the solver room:

```bash
export CONDA_PKGS_DIRS=/tmp/$USER/pkgs CONDA_SOLVER=libmamba
conda build conda/ -c conda-forge -c bioconda \
    --croot /tmp/$USER/croot --no-anaconda-upload
```

Run it inside an allocation with enough memory (resolving bioconda under a tight
`--mem` cap can OOM the solver); request ~32 GB. The packaging itself is
validated independently of this heavyweight build:

- the pip **wheel** builds with all entry points + bundled `data/selftest/`
  (`pip wheel . --no-deps`);
- the run deps **co-resolve** on the channels
  (`mamba create --dry-run -c conda-forge -c bioconda python>=3.8 numpy scipy
  pysam jellyfish samtools` → clean);
- `meta.yaml` renders and `kmate --help` / `--version` / `import kmate` pass.

The authoritative full build is bioconda's CI (a clean container without the NFS
constraint).

## Bioconda submission

The submission recipe + step-by-step instructions live in
[`bioconda/`](bioconda/): [`bioconda/meta.yaml`](bioconda/meta.yaml) (fetches the
**PyPI sdist**, not the 185 MB GitHub auto-tarball) and
[`bioconda/SUBMISSION.md`](bioconda/SUBMISSION.md) (publish the sdist to PyPI →
open the bioconda-recipes PR). `v0.1.0` is already tagged.
