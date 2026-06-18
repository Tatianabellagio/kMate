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

## Bioconda submission (when ready)

1. Tag a release (`v0.1.0`) so GitHub serves a tarball; record its sha256.
2. Swap the `source: path: ../` in `meta.yaml` for the `url:` + `sha256:` lines.
3. Fork `bioconda/bioconda-recipes`, add `recipes/kmate/meta.yaml`, open a PR;
   their CI builds + tests it. Then `conda install -c bioconda kmate` works.
