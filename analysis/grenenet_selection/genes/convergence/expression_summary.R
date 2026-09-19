# Is a gene expressed in 1001T rosettes at all? A truncation or splice change in a gene that
# is silent in every accession does nothing we could measure, and may mark a pseudogene.
# Usage: Rscript expression_summary.R AGI ...   -> results/expression_summary.csv
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
eq <- file.path(normalizePath(file.path(here, "..", "..", "..", "..")), "data", "eqtl")
e <- new.env(); load(file.path(eq, "TG_data_20180606.Rdata"), envir = e)
g <- get("TG.genes", e); meta <- get("TG.meta", e)
keep <- meta$batch_comb != "MU"
raw <- g$raw[, keep]
gene_med <- apply(raw, 1, median)
pct <- ecdf(gene_med)
nm <- as.character(g$names); sel <- match(args, nm)
out <- data.frame(gene = args,
                  median_raw = ifelse(is.na(sel), NA, gene_med[sel]),
                  genome_pctile = ifelse(is.na(sel), NA, round(100 * pct(gene_med[sel]), 1)),
                  frac_accessions_expressed = ifelse(is.na(sel), NA,
                      round(rowMeans(raw[sel[!is.na(sel)], , drop = FALSE] > 5)[match(args, nm[sel[!is.na(sel)]])], 3)))
cat("genome: median of per-gene median raw =", median(gene_med), "\n")
write.csv(out, file.path(here, "results", "expression_summary.csv"), row.names = FALSE)
print(out)
