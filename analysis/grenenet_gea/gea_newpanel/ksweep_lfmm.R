#!/usr/bin/env Rscript
# K-sweep LFMM calibration: for one class, loop over axes x K, fit lfmm_ridge +
# lfmm_test, record the RAW GIF (= raw inflation lambda, pv$gif) and dump the RAW
# per-record p-values (float32) so downstream hit-calling needs no refit.
# Usage: ksweep_lfmm.R <cls>
suppressMessages(library(lfmm))
a <- commandArgs(trailingOnly = TRUE)
cls <- a[1]

ROOT <- "/global/scratch/users/tbellg/kmate/results/grenenet_gea"
LFMM <- file.path(ROOT, "gea_newpanel", "lfmm_site")
ENVD <- file.path(ROOT, "gea_newpanel", "env_site")
OUTD <- file.path(ROOT, "gea_newpanel", "ksweep_calib")
RAWD <- file.path(OUTD, "raw_p")

AXES <- c("bio5", "pc1", "bio12", "bio13", "bio16", "bio19")
KS   <- 1:6

stem <- file.path(LFMM, sprintf("lfmm_%s_site", cls))
dims <- scan(paste0(stem, "_dims.txt"), quiet = TRUE); np <- dims[1]; nr <- dims[2]
con <- file(paste0(stem, "_Y.f64"), "rb")
v <- readBin(con, what = "double", n = np * nr, size = 8, endian = "little"); close(con)
Y <- matrix(v, nrow = np, ncol = nr, byrow = TRUE); rm(v)
cat(sprintf("[%s] Y = %d x %d\n", cls, np, nr)); flush.console()

rows <- list()
for (axis in AXES) {
  X <- as.matrix(read.csv(file.path(ENVD, sprintf("env_site_%s.csv", axis)), header = TRUE))
  for (K in KS) {
    t0 <- Sys.time()
    mod <- lfmm_ridge(Y = Y, X = X, K = K)
    pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
    gif <- as.numeric(pv$gif)
    praw <- as.numeric(pv$pvalue)
    # dump raw p as float32
    fcon <- file(file.path(RAWD, sprintf("%s_%s_K%d.f32", cls, axis, K)), "wb")
    writeBin(as.numeric(praw), fcon, size = 4, endian = "little"); close(fcon)
    n1e5 <- sum(praw < 1e-5, na.rm = TRUE)
    n1e6 <- sum(praw < 1e-6, na.rm = TRUE)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rows[[length(rows) + 1]] <- data.frame(cls = cls, axis = axis, K = K,
      gif = gif, raw_p_lt_1e5 = n1e5, raw_p_lt_1e6 = n1e6, secs = round(dt, 1))
    cat(sprintf("  %s %s K=%d  GIF=%.4f  raw p<1e-5=%d  p<1e-6=%d  (%.0fs)\n",
                cls, axis, K, gif, n1e5, n1e6, dt)); flush.console()
  }
}
df <- do.call(rbind, rows)
write.csv(df, file.path(OUTD, sprintf("lambda_table_%s.csv", cls)), row.names = FALSE)
cat(sprintf("[%s] DONE -> lambda_table_%s.csv\n", cls, cls))
