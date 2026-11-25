# Cross-check of `_sv_burden_onestage.R`: the same models A and B fitted with lme4qtl::relmatLmer
# (kinship passed as relmat) instead of the hand-built lme4 modular fit. Estimates should agree.
# Env: r_env (lme4, lme4qtl). Writes results/sv_adaptive/onestage/lme4qtl_coef_{A,B}.csv
suppressMessages({library(lme4); library(lme4qtl)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
D <- file.path(here, "results", "sv_adaptive", "onestage")
L <- read.csv(file.path(D, "long.csv"), colClasses = c(founder = "character"))
K <- as.matrix(read.csv(file.path(D, "kinship.csv"), row.names = 1, check.names = FALSE))
L$kin <- L$founder; L$fid <- factor(L$founder); L$garden <- factor(L$garden)
L$bio1c <- L$bio1 - mean(unique(L[, c("garden", "bio1")])$bio1)
TERMS  <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter", "intergenic"), paste, sep = "_"))
TESTED <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter"), paste, sep = "_"))
fixed <- paste(c(TERMS, "cactus"), collapse = " + ")
re_base <- "(1 | garden) + (1 | kin) + (1 | fid)"
re_slopes <- paste(sprintf("(0 + %s | garden)", TESTED), collapse = " + ")
inter <- paste(sprintf("%s:bio1c", TESTED), collapse = " + ")
coefs <- function(m) { s <- summary(m)$coefficients; data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], p = 2 * pnorm(-abs(s[, 3])), row.names = NULL) }
t0 <- Sys.time()
mA <- relmatLmer(as.formula(sprintf("y ~ %s + %s", fixed, re_base)), data = L, relmat = list(kin = K))
mB <- relmatLmer(as.formula(sprintf("y ~ %s + %s + %s + %s", fixed, inter, re_base, re_slopes)), data = L, relmat = list(kin = K))
write.csv(coefs(mA), file.path(D, "lme4qtl_coef_A.csv"), row.names = FALSE)
write.csv(coefs(mB), file.path(D, "lme4qtl_coef_B.csv"), row.names = FALSE)
cat(sprintf("lme4qtl A and B fitted in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
