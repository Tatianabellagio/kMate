#!/usr/bin/env Rscript
# Site-level LFMM across many axes, loading the (axis-independent) site-mean Δp matrix ONCE.
# For each axis reads <stem>_env_<axis>.csv, refits lfmm_ridge, writes <outdir>/scores_<axis>.csv
# with signed z-score + raw p + gif-calibrated p. Usage:
#   run_lfmm_scores_multiaxis.R <stem> <K> <outdir> <axis1,axis2,...>
suppressMessages(library(lfmm))
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; K <- as.integer(a[2]); outdir <- a[3]; axes <- strsplit(a[4], ",")[[1]]

dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)   # n_site n_rec
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)      # C-order [n_site x n_rec]
cat("Y (site x rec):", dim(Y), " K =", K, " axes:", length(axes), "\n")

for (ax in axes) {
  X <- as.matrix(read.csv(paste0(stem, "_env_", ax, ".csv"), header = TRUE))
  mod <- lfmm_ridge(Y = Y, X = X, K = K)
  pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = NULL)   # NO gif correction: raw p only
  z   <- as.numeric(pv$score)
  lam <- median(z^2, na.rm = TRUE) / qchisq(0.5, 1)              # inflation DIAGNOSTIC (reported, NOT applied)
  df  <- data.frame(score = z, praw = as.numeric(pv$pvalue))
  write.csv(df, paste0(outdir, "/scores_", ax, ".csv"), row.names = FALSE)
  cat(sprintf("  %-5s  lambda(diag, NOT applied)=%.3f  raw p<1e-5=%d\n", ax, lam, sum(pv$pvalue < 1e-5, na.rm = TRUE)))
}
cat("done\n")
