#!/usr/bin/env Rscript
# LFMM ridge + test, writing BOTH the raw (uncalibrated) and GIF-calibrated p per
# record, plus the GIF. Lets us quantify how much GIF calibration costs in power.
# Usage: run_lfmm_both_p.R <stem> <env_csv> <K> <out_csv>
#   <stem>_Y.f64 (C-order [n x rec] float64), <stem>_dims.txt ("n rec")
suppressMessages(library(lfmm))
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; envf <- a[2]; K <- as.integer(a[3]); outc <- a[4]
dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE); np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little"); close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)
X <- as.matrix(read.csv(envf, header = TRUE))
mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
write.csv(data.frame(pval_raw = as.numeric(pv$pvalue),
                     pval_gif = as.numeric(pv$calibrated.pvalue)), outc, row.names = FALSE)
writeLines(as.character(pv$gif), sub("\\.csv$", ".gif.txt", outc))
cat(sprintf("n=%d rec=%d K=%d  GIF=%.3f  raw p<1e-5=%d  gif p<1e-5=%d -> %s\n",
            np, nr, K, pv$gif, sum(pv$pvalue < 1e-5, na.rm = TRUE),
            sum(pv$calibrated.pvalue < 1e-5, na.rm = TRUE), outc))
