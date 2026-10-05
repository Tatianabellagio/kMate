# Long-format 1001T expression for carrier-level tests: gene-level raw + d_log2_batch,
# isoform-level raw (TG.trans), and a random-gene control panel (gene-level d_log2_batch).
# Usage: Rscript export_expression_long.R <n_random> AGI ...
# v2: writes data/eqtl_v2/expr_{genes,isoforms,control}.csv (gitignored) -- a SEPARATE directory so the old list's data/eqtl exports are untouched
args <- commandArgs(trailingOnly = TRUE); nrand <- as.integer(args[1]); agis <- args[-1]
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(here, "..", "..", "..", "..", ".."))      # kmate root (this script is one directory deeper than the original)
eq <- file.path(root, "data", "eqtl")                                       # source: TG_data_20180606.Rdata
out <- file.path(root, "data", "eqtl_v2"); dir.create(out, showWarnings = FALSE)
e <- new.env(); load(file.path(eq, "TG_data_20180606.Rdata"), envir = e)
g <- get("TG.genes", e); tr <- get("TG.trans", e); m <- get("TG.meta", e)
keep <- m$batch_comb != "MU"; acc <- m$index[keep]
gn <- as.character(g$names); tn <- as.character(tr$names)
long <- function(mat, ids) data.frame(id = rep(ids, times = ncol(mat)),
                                      acc = rep(acc, each = length(ids)), value = as.vector(mat))
sel <- match(agis, gn); sel <- sel[!is.na(sel)]
G1 <- long(g$raw[sel, keep, drop = FALSE], gn[sel]); names(G1)[3] <- "raw"
G1$dlog <- as.vector(g$d_log2_batch[sel, keep, drop = FALSE])
write.csv(G1, file.path(out, "expr_genes.csv"), row.names = FALSE)
ti <- which(sub("\\.[0-9]+$", "", tn) %in% agis)
write.csv(long(tr$raw[ti, keep, drop = FALSE], tn[ti]), file.path(out, "expr_isoforms.csv"), row.names = FALSE)
set.seed(1)
expressed <- which(apply(g$raw[, keep], 1, median) > 5)
ctl <- sample(expressed, nrand)
write.csv(long(g$d_log2_batch[ctl, keep, drop = FALSE], gn[ctl]), file.path(out, "expr_control.csv"), row.names = FALSE)
cat("genes", length(sel), "| isoforms", length(ti), "| control genes", length(ctl), "\n")
