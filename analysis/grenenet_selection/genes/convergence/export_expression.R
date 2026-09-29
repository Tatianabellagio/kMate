# Export 1001T expression (TG.genes$d_log2_batch) for a gene list, so the cis-eQTL can run
# in Python without reloading the 2 GB Rdata each time. Also writes TG.meta.
# Usage: Rscript export_expression.R AGI [AGI ...]
# Output (data/eqtl/, gitignored): cand_expression_d_log2_batch.csv, TG_meta.csv
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
proj <- normalizePath(file.path(here, "..", "..", "..", ".."))
eq <- file.path(proj, "data", "eqtl")
e <- new.env(); load(file.path(eq, "TG_data_20180606.Rdata"), envir = e)
meta <- get("TG.meta", e); genes <- get("TG.genes", e)
write.csv(meta, file.path(eq, "TG_meta.csv"), row.names = FALSE)
gn <- as.character(genes$names)
sel <- match(args, gn)
cat("found", sum(!is.na(sel)), "of", length(args), "; missing:", paste(args[is.na(sel)], collapse = ", "), "\n")
E <- genes$d_log2_batch[sel[!is.na(sel)], , drop = FALSE]
rownames(E) <- args[!is.na(sel)]; colnames(E) <- meta$index
write.csv(E, file.path(eq, "cand_expression_d_log2_batch.csv"))
