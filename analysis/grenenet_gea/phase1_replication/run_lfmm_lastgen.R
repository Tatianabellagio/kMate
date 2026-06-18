#!/usr/bin/env Rscript
# LFMM (lfmm_ridge + lfmm_test gif) on a phase-1-replication class matrix.
# Reads the binary Δp built by build_lfmm_input.py (C-order [n_pools x n_rec] float64)
# and the standardized env, runs K-factor ridge LFMM, writes calibrated p per record.
#
# Usage: run_lfmm_lastgen.R <stem> <K> <out_pval_csv>
#   <stem>: prefix used by build_lfmm_input.py -> <stem>_Y.f64, _dims.txt, _env.csv
suppressMessages({library(lfmm)})
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; K <- as.integer(a[2]); outp <- a[3]

dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)   # n_pools n_rec
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
# C-order [np x nr] -> R matrix filled by row
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)
X <- as.matrix(read.csv(paste0(stem, "_env.csv"), header = TRUE))
cat("Y (pools x rec):", dim(Y), " X:", dim(X), " K =", K, "\n")

mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
cat("GIF (lambda) =", round(pv$gif, 4),
    " | calibrated p<1e-5:", sum(pv$calibrated.pvalue < 1e-5, na.rm = TRUE), "\n")
write.csv(data.frame(pval = as.numeric(pv$calibrated.pvalue)), outp, row.names = FALSE)
cat("wrote", outp, "\n")
