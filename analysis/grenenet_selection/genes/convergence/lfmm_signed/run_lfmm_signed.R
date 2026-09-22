#!/usr/bin/env Rscript
# The multi-axis GEA's LFMM, re-run ONLY to keep the sign it discarded.
# Identical to r2_gea_nonsnp/phase1_replication/multiaxis/run_lfmm_nogif.R -- same
# lfmm_ridge (K = 16), same lfmm_test, same raw (uncalibrated) p -- except that it also
# writes the effect size $B and the z-score $score, which carry the direction of the
# association. lfmm_ridge is deterministic, so the p-values must reproduce the saved ones;
# the calling job checks that before anything is used.
# Usage: run_lfmm_signed.R <stem> <K> <out_csv>
suppressMessages({library(lfmm)})
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; K <- as.integer(a[2]); outp <- a[3]
dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)
X <- as.matrix(read.csv(paste0(stem, "_env.csv"), header = TRUE))
mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
write.csv(data.frame(pval = as.numeric(pv$pvalue), B = as.numeric(pv$B),
                     z = as.numeric(pv$score)), outp, row.names = FALSE)
cat("wrote", outp, "\n")
