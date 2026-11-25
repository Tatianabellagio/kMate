# Do insertions and deletions have the same effect? -- the test that decides whether pooling them
# (as in _sv_burden_pooled_models.R) is legitimate.
#
# Pooling INS and DEL into one "SV burden" is only correct if their effects are equal; otherwise the
# pooled coefficient is a weighted average of two different things, with weights set by the relative
# number of events (2.4x more insertions than deletions), not by biology. This runs the contrast
# directly on the SPLIT models that are already fitted:
#
#   Model 1 (REML):  H0: beta_ins,c = beta_del,c        for c in {CDS, TE, promoter}
#   Model 2 (REML):  H0: gamma_ins,c = gamma_del,c      (the bio1c interactions), same contexts
#
# Wald test of (beta_ins - beta_del) using the fitted coefficients and their covariance:
#   se = sqrt(V_ii + V_jj - 2 V_ij),  z = (b_i - b_j) / se.
# The covariance term matters -- the two burdens are correlated across founders, so treating the
# coefficients as independent would overstate the difference.
#
# Reads the saved fits of _sv_burden_final_models.R (coef_*.csv, vcov_*.csv); fits nothing itself.
# Env: r_env (base R only). Output: results/sv_adaptive/final_models_pooled/class_contrast.csv
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
# FINAL_IN / FINAL_OUT default to the standard shared-variant fits, so existing invocations are unchanged.
# FINAL_IN=final_models_sv6_all FINAL_OUT=final_models_pooled_all runs the same contrast on the
# no-origin-filter split fits (see sv_burden_pooled_nofilter.ipynb).
FM <- file.path(here, "results", "sv_adaptive", Sys.getenv("FINAL_IN", "final_models"))
OUT <- file.path(here, "results", "sv_adaptive", Sys.getenv("FINAL_OUT", "final_models_pooled"))
dir.create(OUT, showWarnings = FALSE)

contrast <- function(model, a, b) {
  cf <- read.csv(file.path(FM, sprintf("coef_%s.csv", model)))
  V <- as.matrix(read.csv(file.path(FM, sprintf("vcov_%s.csv", model)), row.names = 1, check.names = FALSE))
  if (!(a %in% rownames(V) && b %in% rownames(V))) return(NULL)
  bi <- cf$estimate[match(a, cf$term)]; bj <- cf$estimate[match(b, cf$term)]
  se <- sqrt(V[a, a] + V[b, b] - 2 * V[a, b])
  zz <- (bi - bj) / se
  data.frame(model = model, term_ins = a, term_del = b, est_ins = bi, est_del = bj,
             diff = bi - bj, se_diff = se, z = zz, p = 2 * pnorm(-abs(zz)),
             r_coef = V[a, b] / sqrt(V[a, a] * V[b, b]))
}

res <- NULL
for (ctx in c("CDS", "TE", "promoter")) {
  res <- rbind(res, contrast("M1_reml", paste0("ins_", ctx), paste0("del_", ctx)))
  res <- rbind(res, contrast("M2_reml", paste0("ins_", ctx, ":bio1c"), paste0("del_", ctx, ":bio1c")))
}
res$p_bonferroni <- pmin(1, res$p * nrow(res))
write.csv(res, file.path(OUT, "class_contrast.csv"), row.names = FALSE)
print(format(res, digits = 3))
cat("\nA significant difference means the two classes are NOT exchangeable and the pooled model\n",
    "averages two different effects; a non-significant one means pooling loses little.\n")
