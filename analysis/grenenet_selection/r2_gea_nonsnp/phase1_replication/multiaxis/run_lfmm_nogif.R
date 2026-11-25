#!/usr/bin/env Rscript
# The GEA's LFMM, for both arms: lfmm_ridge with K latent factors, then lfmm_test. p-values are RAW.
#
# calibrate = "gif" is passed only so lfmm_test returns lambda ($gif); $pvalue is the uncalibrated p and
# $calibrated.pvalue is not used. Lambda is recorded, never applied: dividing by it deflates only where lambda > 1,
# and on axes where K = 16 over-corrects (lambda < 1) it would inflate the statistics.
#
# Reads <stem>_Y.f64, <stem>_dims.txt, <stem>_env.csv (build_lfmm_input.py).
# Writes <out_csv>: pval (raw), B (effect of the z-scored axis on delta p), z (score; its sign is the direction),
# one row per column of Y, in order; and lambda to <gif_txt> unless it is "-".
# Environment: lfmm_env (R package lfmm).
# Usage: run_lfmm_nogif.R <stem> <K> <out_csv> <gif_txt|->

suppressMessages({library(lfmm)})
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; K <- as.integer(a[2]); outp <- a[3]; gifout <- a[4]

dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)   # n_pools n_rec
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)      # C-order [np x nr]
X <- as.matrix(read.csv(paste0(stem, "_env.csv"), header = TRUE))
cat("Y (pools x rec):", dim(Y), " X:", dim(X), " K =", K, "\n")

mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")   # $pvalue (raw), $gif, $B, $score
lam <- pv$gif
cat("GIF (lambda, NOT applied) =", round(lam, 4),
    " | RAW p<1e-5:", sum(pv$pvalue < 1e-5, na.rm = TRUE), "\n")
write.csv(data.frame(pval = as.numeric(pv$pvalue), B = as.numeric(pv$B), z = as.numeric(pv$score)),
          outp, row.names = FALSE)                                 # RAW p, no gif
if (!is.na(gifout) && gifout != "-") writeLines(as.character(lam), gifout)
cat("wrote", outp, "(raw p, B, z)\n")
