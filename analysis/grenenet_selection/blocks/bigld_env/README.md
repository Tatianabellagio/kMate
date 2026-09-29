# BigLD env — reproducible install of `gpart::BigLD`

`gpart::BigLD` is the LD-block **fine-split** used by HapFM's block partition
(`recompute_blocks.py` → `HapFM/bin/BigLD.R` → `gpart::BigLD(..., CLQmode="density")`).
This folder makes that install **reproducible**, because the live install is a
hand-built R package (not a conda package), so `conda env export` won't capture it.

## Where it currently lives
Installed into the conda **`r_env`** R library
(`/global/home/users/tbellg/miniforge3/envs/r_env/lib/R/library/gpart`, R 4.4.1),
along with deps `optparse, igraph, data.table, Rcpp, GenomicRanges, IRanges`.
Persistent (NFS home). Day-to-day, just use `r_env`'s `Rscript` with
`export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH`.

## To recreate (e.g. new cluster / rebuilt env)
```bash
bash install_bigld.sh [R_PREFIX_BIN]   # defaults to the r_env bin
```

## Why this is needed (the gotchas)
- **conda/mamba HANGS on Savio's NFS home** (flock failures) — do NOT install R
  packages via conda here. Use R's own installers.
- **gpart is gone from current Bioconductor** (only ≤ Bioc 3.14 / R 4.1), so
  `BiocManager::install("gpart")` fails on R 4.4.1.
- The vendored `gpart_1.12.0_stripped.tar.gz` is the 1.12.0 source with two fixes:
  1. **human-genome deps stripped** — `Depends`/`Imports` reduced to what BigLD/CLQD
     actually use (igraph, Rcpp, data.table, GenomicRanges, IRanges, grid, base);
     the `Homo.sapiens` + `TxDb.Hsapiens.UCSC.hg38.knownGene` imports (only used by
     the unused GPART function) removed from DESCRIPTION + NAMESPACE; leftover blank
     lines removed (R CMD INSTALL rejects them).
  2. **`src/gpart_cpp.cpp`**: `#define PI M_PI` added (gcc ≥ 14 has no bare `PI`).
- Build needs system `libxml2` on the linker path: `export LIBRARY_PATH=/usr/lib64`.

## BigLD.R caller contract
`SNPinfo` must be exactly 3 columns in order `chrN, rsID, bp`; `geno` is a
data.frame with SNPs as columns. Full provenance: memory `gpart-bigld-install`.
