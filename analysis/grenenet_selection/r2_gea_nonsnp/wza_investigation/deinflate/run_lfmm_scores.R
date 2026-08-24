#!/usr/bin/env Rscript
# Re-run lfmm_ridge + lfmm_test and SAVE the signed z-scores (pv$score), so we can
# compare deflation schemes (gif vs median+MAD vs Efron empirical null) downstream.
# Usage: run_lfmm_scores.R <stem> <K> <out_csv>
suppressMessages({library(lfmm)})
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; K <- as.integer(a[2]); outp <- a[3]

dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)      # C-order [np x nr]
X <- as.matrix(read.csv(paste0(stem, "_env.csv"), header = TRUE))
cat("Y:", dim(Y), " X:", dim(X), " K =", K, "\n")

mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = NULL)   # NO gif correction: raw p only
score <- as.numeric(pv$score)        # signed z-score per record
praw  <- as.numeric(pv$pvalue)       # raw two-sided p (uncalibrated)
lam <- median(score^2, na.rm = TRUE) / qchisq(0.5, 1)          # inflation DIAGNOSTIC (reported, NOT applied)
cat("lambda (diagnostic, NOT applied) =", round(lam, 4), "\n")
df <- data.frame(score = score, praw = praw)
write.csv(df, outp, row.names = FALSE)
cat("wrote", outp, " nrow =", nrow(df),
    " | score: median", round(median(score), 4), " sd", round(sd(score), 4), "\n")
