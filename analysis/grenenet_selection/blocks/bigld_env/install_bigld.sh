#!/bin/bash
# Reproducibly install gpart::BigLD (HapFM's LD-block fine-split engine) into a
# target R. We do NOT use conda for this (mamba hangs on Savio's NFS home FS).
#
# gpart was removed from current Bioconductor (only <=3.14 / R 4.1). We vendor the
# 1.12.0 source with two fixes already applied (see README.md):
#   - human-genome Depends/Imports stripped (BigLD/CLQD don't need them)
#   - src/gpart_cpp.cpp: `#define PI M_PI` added (gcc>=14 has no bare PI)
#
# Usage:  bash install_bigld.sh [R_PREFIX_BIN]
#   R_PREFIX_BIN defaults to the conda r_env bin (R 4.4.1).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RBIN="${1:-/global/home/users/tbellg/miniforge3/envs/r_env/bin}"
RSCRIPT="$RBIN/Rscript"; R="$RBIN/R"

echo "== target R: $($RSCRIPT -e 'cat(R.version.string)') =="

# build needs system libxml2 (igraph) + conda compilers on PATH
export PATH="$RBIN:$PATH"
export LIBRARY_PATH="/usr/lib64:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="/usr/lib64:${LD_LIBRARY_PATH:-}"

echo "== deps (R installers, not conda) =="
"$RSCRIPT" -e '
  options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=4)
  for (p in c("Rcpp","igraph","data.table","optparse","BiocManager"))
    if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  for (p in c("GenomicRanges","IRanges"))
    if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
  cat("deps OK\n")'

echo "== build + install vendored stripped gpart =="
TMP="$(mktemp -d)"; tar xzf "$HERE/gpart_1.12.0_stripped.tar.gz" -C "$TMP"
"$R" CMD INSTALL "$TMP/gpart"
rm -rf "$TMP"

echo "== verify =="
"$RSCRIPT" -e 'suppressMessages(library(gpart)); stopifnot(exists("BigLD")); cat("gpart::BigLD OK, v", as.character(packageVersion("gpart")), "\n")'
echo
echo "RUNTIME NOTE: before any Rscript that loads gpart, set:"
echo "  export LD_LIBRARY_PATH=/usr/lib64:\$LD_LIBRARY_PATH"
