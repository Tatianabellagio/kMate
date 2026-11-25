#!/usr/bin/env Rscript
# The GEA's LFMM: lfmm_ridge (K latent factors) + lfmm_test, RAW p (no genomic-control calibration).
# Used by both arms of the production GEA: the record arm (SNP, short indel; multiaxis/ma_lfmm_cr50.sbatch) and the
# event arm (SV / non-SNP Truvari events; genes/convergence/event_gea/run_event_lfmm_cr50.sbatch).
#
# Why raw p: calibrate = "gif" divides z^2 by lambda = median(z^2) / 0.456, which is only a deflation when lambda > 1;
# where K = 16 over-corrects (lambda < 1 on some axes) it would inflate the statistics. Project decision (2026-07-03):
# save raw p ($pvalue; $calibrated.pvalue is not used) and RECORD lambda.
#
# Writes <out_csv>: pval (raw), B (effect of the z-scored axis on delta p) and z (the score; its sign is the direction),
# one row per column of Y, in Y's column order; and lambda to <gif_txt> unless it is "-".
# B and z were added 2026-10-20 (GEA audit) so the event arm uses this script instead of a copy of it
# (lfmm_signed/run_lfmm_signed.R, now in archive_convergence/): on the production SV inputs (sv x bio1, sv x pc3) the two
# give identical pval, B and z, the pval column and lambda are identical to the previous version of this script, and
# both reproduce event_gea/cr50/lfmm_event_sv_gen9_{bio1,pc3}.csv.gz to <= 2.5e-12.
# Usage: run_lfmm_nogif.R <stem> <K> <out_csv> <gif_txt|->      (reads stem_Y.f64, stem_dims.txt, stem_env.csv)
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
