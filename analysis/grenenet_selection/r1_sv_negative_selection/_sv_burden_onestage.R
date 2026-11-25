# One-stage founder burden model on all founder x garden observations (231 x 31 = 7,161 rows).
#
#   y_fs = burden_f . beta  +  (burden_f x bio1_s) . gamma  +  route_f
#          + garden_s  +  kin_f  +  founder_f  +  sum_j b_js burden_fj  +  e_fs
#
#   garden_s ~ N(0, s2_garden)          random garden intercept
#   kin_f    ~ N(0, s2_kin K)           founder genetic effect, K = SNP kinship (231 x 231)
#   founder_f ~ N(0, s2_founder)        the same founder measured in 31 gardens
#   b_js     ~ N(0, s2_j)               random garden slope of tested burden term j (independent)
#
# Kinship is put into lme4 through its modular interface (the approach lme4qtl wraps): the random
# effect (1 | kin) is fitted on a design Z L, where K = L L', so kin = L a with a ~ N(0, s2_kin I).
# Burden terms: SNPs, insertions, deletions x coding / promoter / intergenic, standardised (per SD);
# bio1 centred (per degree C). Inputs from `_sv_burden_onestage_export.py`.
#
# Models:
#   A  main effects only (+ garden, kin, founder)  -> must reproduce the two-stage mean-fitness fit
#   B  A + tested-term x bio1 + random garden slopes for the tested terms   (the direct climate test)
#   C  B without the random slopes                  (how much the slopes matter)
# Likelihood-ratio tests (ML) for each tested-term x bio1 interaction, dropping one at a time from B.
# Env: r_env (lme4). Writes results/sv_adaptive/onestage/{coef_A,coef_B,coef_C,lrt_B,varcomp}.csv
suppressMessages({library(lme4); library(Matrix)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
D <- file.path(here, "results", "sv_adaptive", "onestage")

L <- read.csv(file.path(D, "long.csv"), colClasses = c(founder = "character"))
K <- as.matrix(read.csv(file.path(D, "kinship.csv"), row.names = 1, check.names = FALSE))
L$kin <- factor(L$founder); L$fid <- factor(L$founder); L$garden <- factor(L$garden)
L$bio1c <- L$bio1 - mean(unique(L[, c("garden", "bio1")])$bio1)
K <- K[levels(L$kin), levels(L$kin)]
e <- eigen(K, symmetric = TRUE); Lk <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))   # K = Lk Lk'

TERMS  <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter", "intergenic"), paste, sep = "_"))
TESTED <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter"), paste, sep = "_"))
fixed  <- paste(c(TERMS, "cactus"), collapse = " + ")

fit_kin <- function(form, REML = TRUE) {
  ctl <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE)
  lf <- lFormula(as.formula(form), data = L, REML = REML, control = ctl)
  i <- which(names(lf$reTrms$cnms) == "kin")                       # the (1 | kin) block
  lf$reTrms$Ztlist[[i]] <- as(t(Lk) %*% lf$reTrms$Ztlist[[i]], "CsparseMatrix")
  lf$reTrms$Zt <- do.call(rbind, lf$reTrms$Ztlist)
  dfun <- do.call(mkLmerDevfun, lf)
  opt <- optimizeLmer(dfun)
  mkMerMod(environment(dfun), opt, lf$reTrms, fr = lf$fr)
}
coefs <- function(m) {
  s <- summary(m)$coefficients
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], z = s[, 3], p = 2 * pnorm(-abs(s[, 3])), row.names = NULL)
}
re_base   <- "(1 | garden) + (1 | kin) + (1 | fid)"
re_slopes <- paste(sprintf("(0 + %s | garden)", TESTED), collapse = " + ")
inter     <- paste(sprintf("%s:bio1c", TESTED), collapse = " + ")

t0 <- Sys.time()
fA <- sprintf("y ~ %s + %s", fixed, re_base)
fB <- sprintf("y ~ %s + %s + %s + %s", fixed, inter, re_base, re_slopes)
fC <- sprintf("y ~ %s + %s + %s", fixed, inter, re_base)
mA <- fit_kin(fA); mB <- fit_kin(fB); mC <- fit_kin(fC)
cat(sprintf("fitted A, B, C in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
for (nm in c("A", "B", "C")) write.csv(coefs(get(paste0("m", nm))), file.path(D, sprintf("coef_%s.csv", nm)), row.names = FALSE)
vc <- as.data.frame(VarCorr(mB)); write.csv(vc[, c("grp", "var1", "vcov")], file.path(D, "varcomp.csv"), row.names = FALSE)

# likelihood-ratio tests (ML) for each tested-term x bio1 interaction, one dropped at a time from B
mBml <- fit_kin(fB, REML = FALSE); lrt <- NULL
for (t in TESTED) {
  f0 <- sprintf("y ~ %s + %s + %s + %s", fixed, paste(sprintf("%s:bio1c", setdiff(TESTED, t)), collapse = " + "), re_base, re_slopes)
  m0 <- fit_kin(f0, REML = FALSE); chi <- 2 * (as.numeric(logLik(mBml)) - as.numeric(logLik(m0)))
  lrt <- rbind(lrt, data.frame(term = paste0(t, ":bio1c"), chisq = chi, p = pchisq(chi, 1, lower.tail = FALSE)))
}
write.csv(lrt, file.path(D, "lrt_B.csv"), row.names = FALSE)
cat(sprintf("done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

show <- function(m, lab, rows) { x <- coefs(m); x <- x[x$term %in% rows, ]; cat("\n", lab, "\n"); print(format(x, digits = 3), row.names = FALSE) }
show(mA, "A  main effects (per SD):", TESTED)
show(mB, "B  climate interactions (per SD per degree C), with random garden slopes:", paste0(TESTED, ":bio1c"))
show(mC, "C  climate interactions, WITHOUT random garden slopes:", paste0(TESTED, ":bio1c"))
cat("\nLRT (ML), each interaction dropped from B:\n"); print(format(lrt, digits = 3), row.names = FALSE)
cat("\nvariance components (model B):\n"); print(format(vc[, c("grp", "var1", "vcov")], digits = 3), row.names = FALSE)
