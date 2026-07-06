#!/usr/bin/env Rscript
# LFMM ridge K=16, NO gif calibration (multi-axis clq0.9 GEA).
# Writes the RAW (uncalibrated) lfmm_test p-value. Rationale: gif calibration
# divides the test statistics by lambda = median(z^2)/0.456, which is only a valid
# DEFLATION when lambda>1; when LFMM K=16 over-corrects on an axis (lambda<1, seen
# on some bioclim axes here) calibrate="gif" would INFLATE the statistics and
# manufacture false signal. Per project decision (2026-07-03) we do NOT apply gif
# on the multi-axis run; we save raw p and RECORD lambda per axis for transparency.
#
# Usage: run_lfmm_nogif.R <stem> <K> <out_pval_csv> <gif_out_txt>
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
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")  # gives $pvalue (raw) + $gif
lam <- pv$gif
cat("GIF (lambda, NOT applied) =", round(lam, 4),
    " | RAW p<1e-5:", sum(pv$pvalue < 1e-5, na.rm = TRUE), "\n")
write.csv(data.frame(pval = as.numeric(pv$pvalue)), outp, row.names = FALSE)  # RAW, no gif
if (!is.na(gifout)) writeLines(as.character(lam), gifout)
cat("wrote", outp, "(raw, uncalibrated)\n")
