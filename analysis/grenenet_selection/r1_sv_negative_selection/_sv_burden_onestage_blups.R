# Model B of `_sv_burden_onestage.R` (burden x garden bio1 + random garden slopes), refitted with lme4qtl
# (validated identical to the hand-built fit) to save what the figures need:
#   vcov_B.csv    covariance matrix of the fixed effects (confidence bands of the fitted lines)
#   garden_B.csv  per garden and tested term: the garden's total effect  beta + gamma * bio1c + b_s
#                 (fixed part + the garden's random slope, its BLUP) and the random-slope conditional SD
# Env: r_env (lme4, lme4qtl). Writes results/sv_adaptive/onestage/{vcov_B,garden_B}.csv
suppressMessages({library(lme4); library(lme4qtl)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
D <- file.path(here, "results", "sv_adaptive", "onestage")
L <- read.csv(file.path(D, "long.csv"), colClasses = c(founder = "character"))
K <- as.matrix(read.csv(file.path(D, "kinship.csv"), row.names = 1, check.names = FALSE))
L$kin <- L$founder; L$fid <- factor(L$founder); L$garden <- factor(L$garden)
b1 <- unique(L[, c("garden", "bio1")]); L$bio1c <- L$bio1 - mean(b1$bio1)
TERMS  <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter", "intergenic"), paste, sep = "_"))
TESTED <- c(outer(c("snp", "ins", "del"), c("CDS", "promoter"), paste, sep = "_"))
f <- sprintf("y ~ %s + %s + (1 | garden) + (1 | kin) + (1 | fid) + %s", paste(c(TERMS, "cactus"), collapse = " + "),
             paste(sprintf("%s:bio1c", TESTED), collapse = " + "), paste(sprintf("(0 + %s | garden)", TESTED), collapse = " + "))
m <- relmatLmer(as.formula(f), data = L, relmat = list(kin = K))
V <- as.matrix(vcov(m)); write.csv(V, file.path(D, "vcov_B.csv"))
fe <- fixef(m); re <- ranef(m, condVar = TRUE)$garden; out <- NULL
b1$bio1c <- b1$bio1 - mean(b1$bio1)
for (t in TESTED) {
  # ranef(m)$garden is one data frame per slope term; find the column for t
  # all terms on the grouping factor 'garden' come back as columns of one data frame; its condVar attribute
  # is either a 3-d array (terms x terms x gardens) or a list with one array per term
  blk <- re; k <- match(t, names(re)); pv <- attr(re, "postVar"); sdv <- rep(NA_real_, nrow(re))
  if (is.list(pv) && length(pv) >= k) { a_ <- pv[[k]]; sdv <- sqrt(a_[1, 1, ]) } else if (is.array(pv) && dim(pv)[1] >= k) sdv <- sqrt(pv[k, k, ])
  g <- rownames(blk); x <- b1$bio1c[match(g, b1$garden)]
  out <- rbind(out, data.frame(term = t, garden = g, bio1 = b1$bio1[match(g, b1$garden)],
                               fixed_part = fe[t] + fe[paste0(t, ":bio1c")] * x, blup = blk[[t]], cond_sd = sdv,
                               total = fe[t] + fe[paste0(t, ":bio1c")] * x + blk[[t]]))
}
write.csv(out, file.path(D, "garden_B.csv"), row.names = FALSE)
cat("saved vcov_B.csv and garden_B.csv\n")
