#!/usr/bin/env Rscript
# LFMM K-sweep for calibration: rerun lfmm_ridge + lfmm_test(gif) over a range of K
# on the SNP Δp matrix and record the genomic inflation factor (GIF) and p-value
# calibration per K. The well-calibrated K is the smallest K with GIF ~ 1 and a flat
# null p-value histogram. Loads Y once, sweeps K.
#
# Usage: run_lfmm_ksweep.R <stem> <out_csv> <comma_K_list>
#   <stem>: build_lfmm_input.py prefix -> <stem>_Y.f64, _dims.txt, _env.csv
suppressMessages(library(lfmm))
a <- commandArgs(trailingOnly = TRUE)
stem <- a[1]; outcsv <- a[2]; Ks <- as.integer(strsplit(a[3], ",")[[1]])

dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE)   # n_pools n_rec
np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little")
close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE)      # C-order [np x nr]
X <- as.matrix(read.csv(paste0(stem, "_env.csv"), header = TRUE))
cat("Y (pools x rec):", dim(Y), " X:", dim(X), " sweeping K =", paste(Ks, collapse = ","), "\n")

res <- data.frame()
for (K in Ks) {
  t0 <- Sys.time()
  mod <- lfmm_ridge(Y = Y, X = X, K = K)
  pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
  p   <- pv$calibrated.pvalue
  # null-flatness diagnostic: fraction of p in the top decile (should be ~0.10 if calibrated)
  frac_top <- mean(p > 0.9, na.rm = TRUE)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("K=%2d  GIF=%.4f  p<1e-5=%d  p<1e-3=%d  frac(p>0.9)=%.3f  [%ss]\n",
              K, pv$gif, sum(p < 1e-5, na.rm = TRUE), sum(p < 1e-3, na.rm = TRUE), frac_top, dt))
  res <- rbind(res, data.frame(K = K, gif = as.numeric(pv$gif),
                               n_p_lt_1e5 = sum(p < 1e-5, na.rm = TRUE),
                               n_p_lt_1e3 = sum(p < 1e-3, na.rm = TRUE),
                               frac_top_decile = frac_top, secs = dt))
  write.csv(res, outcsv, row.names = FALSE)   # incremental write so partial results survive
}
cat("wrote", outcsv, "\n")
