#!/usr/bin/env Rscript
# LFMM ridge climate-GEA on the gen-3 SV Δp matrix (mirrors the phase-1 method:
# lfmm_ridge + lfmm_test with GIF calibration). Reports the genomic inflation
# factor (GIF) so we can pick K where lambda ~ 1.
#
# Usage: run_lfmm.R <delta_p.csv> <env.csv> <K> <out_prefix>
#   delta_p.csv : loci x pools (header = pool names)   [Y, transposed here]
#   env.csv     : one standardized column (the bioclim var), pool order matches
#   K           : number of latent factors
#   out_prefix  : writes <prefix>.calibrated_pval.csv, .pval.csv, .beta.csv, .gif.txt
suppressMessages({library(lfmm); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
dp_file <- args[1]; env_file <- args[2]; K <- as.integer(args[3]); out <- args[4]

Y <- as.matrix(fread(dp_file, header = TRUE))   # [loci x pools]
Y <- t(Y)                                       # -> [pools x loci]
X <- as.matrix(read.csv(env_file, header = TRUE))
cat("Y (pools x loci):", dim(Y), " X (pools x 1):", dim(X), " K =", K, "\n")

mod <- lfmm_ridge(Y = Y, X = X, K = K)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")

write.csv(pv$calibrated.pvalue, paste0(out, ".calibrated_pval.csv"), row.names = FALSE)
write.csv(pv$pvalue,            paste0(out, ".pval.csv"),            row.names = FALSE)
write.csv(pv$B,                 paste0(out, ".beta.csv"),           row.names = FALSE)
writeLines(as.character(pv$gif), paste0(out, ".gif.txt"))
cat("K =", K, " GIF (lambda) =", round(pv$gif, 4),
    " | calibrated p<1e-5:", sum(pv$calibrated.pvalue < 1e-5, na.rm = TRUE), "\n")
