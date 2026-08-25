#!/usr/bin/env Rscript
# susie_rss on one block. args: z_csv (one per line)  R_txt (space-sep m x m)  n  out_csv
suppressMessages(library(susieR))
a <- commandArgs(trailingOnly = TRUE)
z <- scan(a[1], quiet = TRUE)
m <- length(z)
R <- matrix(scan(a[2], quiet = TRUE), nrow = m, byrow = TRUE)
n <- as.integer(a[3])
fit <- tryCatch(
  susie_rss(z = z, R = R, n = n, L = 10, estimate_residual_variance = FALSE),
  error = function(e) { message("susie error: ", conditionMessage(e)); NULL })
if (is.null(fit)) { write.csv(data.frame(idx = seq_len(m), pip = NA, cs = 0L), a[4], row.names = FALSE); quit(status = 0) }
pip <- fit$pip
cs <- integer(m)
if (!is.null(fit$sets$cs)) for (i in seq_along(fit$sets$cs)) cs[fit$sets$cs[[i]]] <- i
write.csv(data.frame(idx = seq_len(m), pip = pip, cs = cs), a[4], row.names = FALSE)
