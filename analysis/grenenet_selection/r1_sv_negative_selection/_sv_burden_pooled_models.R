# Sensitivity analysis: the founder burden models with insertions and deletions POOLED per context.
#
# Why. INS and DEL are defined relative to Col-0, not ancestral/derived (the SFS in Figure 1 is folded
# for exactly that reason), so a reviewer can reasonably ask whether splitting them is justified or
# whether one "SV burden" per context is the honest unit. This refits the final models with
# sv_<ctx> = (ins carried + del carried) / (ins called + del called), everything else identical to
# _sv_burden_final_models.R: same response, same covariates, same kinship, same random-effect
# structure, same REML/ML and LRT scheme.
#
# Read it as a sensitivity analysis, NOT a replacement: pooling is only the right unit if the two
# classes have the same effect. The companion script `_sv_burden_class_contrast.R` tests that
# directly (H0: beta_ins = beta_del within a context) on the split model; if the classes differ, the
# pooled coefficient is a weighted average of two different effects, with weights set by how many
# events of each class exist (2.4x more insertions than deletions), not by biology.
#
# Model 1p: y ~ snp_CDS + sv_CDS + sv_TE + sv_promoter + sv_intergenic + cactus
#               + (1 | garden) + (1 | kin) + (1 | fid)
#               + (0 + sv_CDS | garden) + (0 + sv_TE | garden) + (0 + sv_promoter | garden)
# Model 2p: Model 1p + {sv_CDS, sv_TE, sv_promoter}:bio1c
#
# Tasks (TASK env or SLURM array): 0 M1p_reml  1 M2p_reml  2 M2p_ml  3-5 M2p_ml_drop_<term>
# Burden rates are rebuilt from results/sv_adaptive/founder_burden_event.csv (the `shared_*_all`
# columns: shared variants, >= 2 carriers, called in >= 90% of founders) and z-scored across the 231
# founders -- verified to reproduce long.csv's split columns at r = 1.000000 before pooling.
# Env: r_env (lme4, lme4qtl). Output: results/sv_adaptive/final_models_pooled/.
suppressMessages({library(lme4); library(lme4qtl)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
# ONESTAGE_DIR / FINAL_OUT / BURDEN_ORIGIN all default to the standard shared-variant run, so existing
# invocations are unchanged. BURDEN_ORIGIN=all pools the no-origin-filter burdens instead, to mirror
# `_sv_burden_final_models.R`'s BURDEN_ORIGIN arm (see sv_burden_origin_sensitivity.ipynb).
IN <- file.path(here, "results", "sv_adaptive", Sys.getenv("ONESTAGE_DIR", "onestage"))
RAW <- file.path(here, "results", "sv_adaptive", "founder_burden_event.csv")
OUT <- file.path(here, "results", "sv_adaptive", Sys.getenv("FINAL_OUT", "final_models_pooled"))
ORIGIN <- Sys.getenv("BURDEN_ORIGIN", "shared")
stopifnot(ORIGIN %in% c("shared", "all"))
dir.create(OUT, showWarnings = FALSE)
task <- as.integer(Sys.getenv("TASK", Sys.getenv("SLURM_ARRAY_TASK_ID", "0")))

L <- read.csv(file.path(IN, "long.csv"), colClasses = c(founder = "character"))
K <- as.matrix(read.csv(file.path(IN, "kinship.csv"), row.names = 1, check.names = FALSE))
B <- read.csv(RAW, colClasses = c(founder = "character"))

# SUBSET: fit on one sequencing-technology panel only ("cactus" or "shortread"); "" = all 231 founders.
# Burden rates are z-scored WITHIN the subset and pooled_rate_scale.csv is written for that subset, so the
# per-carried-variant conversion stays internally consistent. `cactus` is dropped from the fixed effects
# because it is constant within a panel, and the kinship matrix is subset to the retained founders.
SUBSET <- Sys.getenv("SUBSET", "")
stopifnot(SUBSET %in% c("", "cactus", "shortread"))
if (nzchar(SUBSET)) {
  keep <- B$founder[if (SUBSET == "cactus") B$is_cactus == "True" | B$is_cactus == TRUE
                    else !(B$is_cactus == "True" | B$is_cactus == TRUE)]
  cat(sprintf("SUBSET=%s -> %d founders retained of %d\n", SUBSET, length(keep), nrow(B)))
  stopifnot(length(keep) > 0)
  B <- B[B$founder %in% keep, ]
  L <- L[L$founder %in% keep, ]
  K <- K[as.character(keep), as.character(keep), drop = FALSE]
}

# POOL_SET fixes which burden terms exist and which genomic contexts feed each one. Contexts named in no
# group contribute to NO burden -- that is deliberate for sv4 (reference TE) but was an oversight for UTR
# and intron, which are 20% of passing events and entered nothing until `all6` and `cds2` were added.
CTX6 <- c("CDS", "UTR", "intron", "TE", "promoter", "intergenic")
POOL_SET0 <- Sys.getenv("POOL_SET", "sv6")
stopifnot(POOL_SET0 %in% c("sv4", "sv6", "all6", "cds2"))
GRP <- switch(POOL_SET0,
  sv4  = list(CDS = "CDS", promoter = "promoter", intergenic = "intergenic"),
  sv6  = list(CDS = "CDS", TE = "TE", promoter = "promoter", intergenic = "intergenic"),
  all6 = list(CDS = "CDS", UTR = "UTR", intron = "intron", TE = "TE",
              promoter = "promoter", intergenic = "intergenic"),
  cds2 = list(CDS = "CDS", nonCDS = c("UTR", "intron", "TE", "promoter", "intergenic")))
CTX <- names(GRP)
z <- function(v) (v - mean(v)) / sd(v)
P <- data.frame(founder = B$founder)
cnt <- list()
# BURDEN_STRATIFY=freq: split the contexts named in STRATIFY_CTX into ALT-frequency bins (FREQ_BINS),
# reading the dedicated table from `_sv_burden_freq_bins.py`. Contexts not named stay pooled over all
# frequencies, so the model stays as close as possible to the unstratified one. Empty by default.
# NOTE the frequency table has no origin classes: it is `all` variants only, matching BURDEN_ORIGIN=all.
# NOTE `private` (MAC 1) is only interpretable within the cactus panel -- see the header of
# _sv_burden_freq_bins.py. Fitting it on all 231 founders measures assembly route, not biology.
STRATIFY <- Sys.getenv("BURDEN_STRATIFY", "")
stopifnot(STRATIFY %in% c("", "freq", "size"))
DEFBINS <- if (identical(STRATIFY, "size")) "small,medium,large" else "rare,common"
FBINS <- setdiff(strsplit(Sys.getenv("FREQ_BINS", DEFBINS), "[,: ]+")[[1]], "")
SCTX  <- setdiff(strsplit(Sys.getenv("STRATIFY_CTX", "CDS"), "[,: ]+")[[1]], "")
if (nzchar(STRATIFY)) {
  stopifnot(ORIGIN == "all", all(SCTX %in% CTX), POOL_SET0 %in% c("sv4", "sv6"))
  FR <- read.csv(file.path(here, "results", "sv_adaptive",
                           sprintf("founder_burden_event_%s.csv", STRATIFY)),
                 colClasses = c(founder = "character"))
  FR <- FR[match(B$founder, FR$founder), ]
  stopifnot(identical(FR$founder, B$founder))
  cat(sprintf("BURDEN_STRATIFY=%s: contexts %s split into bins %s\n", STRATIFY,
              paste(SCTX, collapse = "+"), paste(FBINS, collapse = "+")))
}
for (ctx in CTX) {
  if (nzchar(STRATIFY) && ctx %in% SCTX) {
    for (fb in FBINS) {
      carried <- FR[[paste0("ins_", ctx, "_", fb)]] + FR[[paste0("del_", ctx, "_", fb)]]
      called  <- FR[[paste0("called_ins_", ctx, "_", fb)]] + FR[[paste0("called_del_", ctx, "_", fb)]]
      stopifnot(all(called > 0))
      rate <- carried / called
      P[[paste0("sv_", ctx, "_", fb)]] <- z(rate)
      cnt[[paste0(ctx, "_", fb)]] <- c(sd_rate = sd(rate), med_called = median(called),
                                       med_carried = median(carried))
    }
    next
  }
  carried <- 0; called <- 0
  for (cc in GRP[[ctx]]) {                       # a group may pool several contexts (e.g. cds2's nonCDS)
    carried <- carried + B[[paste0("ins_", ORIGIN, "_", cc, "_all")]] + B[[paste0("del_", ORIGIN, "_", cc, "_all")]]
    called  <- called  + B[[paste0("called_ins_", ORIGIN, "_", cc, "_all")]] +
                         B[[paste0("called_del_", ORIGIN, "_", cc, "_all")]]
  }
  rate <- carried / called
  P[[paste0("sv_", ctx)]] <- z(rate)
  # kept for the per-variant conversion: SD of the raw rate and the median called count
  cnt[[ctx]] <- c(sd_rate = sd(rate), med_called = median(called), med_carried = median(carried))
}
# The per-variant scale must be tagged with everything that changes what the burdens ARE: a stratified run
# redefines the contexts (sv_CDS -> sv_CDS_rare etc.), so sharing one filename let each arm silently clobber
# the previous arm's scale. RE_DROP is excluded on purpose -- it changes only the random effects.
EXFIX <- Sys.getenv("EXTRA_FIXED", "")
SCALESUF <- paste0(if (nzchar(EXFIX)) paste0("_x_", gsub("[^A-Za-z0-9]+", "", EXFIX)) else "",
                   if (!(POOL_SET0 %in% c("sv4", "sv6"))) paste0("_", POOL_SET0) else "",
                   if (nzchar(SUBSET)) paste0("_", SUBSET) else "",
                   if (nzchar(STRATIFY)) paste0("_freq", gsub("[^A-Za-z0-9]+", "",
                                                paste(c(SCTX, FBINS), collapse = ""))) else "")
write.csv(do.call(rbind, lapply(names(cnt), function(x) data.frame(context = x, t(cnt[[x]])))),
          file.path(OUT, paste0("pooled_rate_scale", SCALESUF, ".csv")), row.names = FALSE)

# founder ORIGIN climate (B$bio1 is the accession's home bio1; L$bio1 is the GARDEN's) -- carried so that
# EXTRA_FIXED can add an origin x garden-climate term, the local-adaptation control.
P$obio1 <- B$bio1
# founder genotyping quality, so EXTRA_FIXED can add `callrate_sv`: within the cactus panel the rare-bin
# burden tracks it (r = -0.44, p = 4e-05) while private and common do not, so a stratified fit should control it
P$callrate_sv <- B$callrate_sv
L <- merge(L, P, by = "founder")
stopifnot(nrow(L) == length(unique(L$founder)) * length(unique(L$garden)))
cat(sprintf("fitting on %d founders x %d gardens = %d rows\n", length(unique(L$founder)), length(unique(L$garden)), nrow(L)))
L$kin <- L$founder; L$fid <- factor(L$founder); L$garden <- factor(L$garden)
b1 <- unique(L[, c("garden", "bio1")]); L$bio1c <- L$bio1 - mean(b1$bio1); b1$bio1c <- b1$bio1 - mean(b1$bio1)

# POOL_SET: sv6 (default) keeps the reference-TE context term; sv4 drops it entirely, so TE-context
# variants contribute to NO burden -- the pooled analogue of _sv_burden_final_models.R's SV_SET=sv4 arm.
POOL_SET <- POOL_SET0
L$obio1c <- L$obio1 - mean(unique(L[, c("founder", "obio1")])$obio1)
# expand any stratified context into its per-bin terms, keeping the original ordering
expand <- function(v) unlist(lapply(v, function(t) {
  ctx <- sub("^sv_", "", t)
  if (nzchar(STRATIFY) && ctx %in% SCTX) paste0("sv_", ctx, "_", FBINS) else t
}))
# which burdens get a garden-temperature interaction and a random garden slope
SV3 <- expand(switch(POOL_SET,
  sv4  = c("sv_CDS", "sv_promoter"),
  sv6  = c("sv_CDS", "sv_TE", "sv_promoter"),
  all6 = c("sv_CDS", "sv_TE", "sv_promoter"),
  cds2 = c("sv_CDS", "sv_nonCDS")))
ALLB <- expand(paste0("sv_", CTX))
fixed <- paste(c("snp_CDS", ALLB, "cactus"), collapse = " + ")
# `cactus` is constant within a technology panel -> unidentifiable, drop it
if (nzchar(SUBSET)) fixed <- sub(" \\+ cactus", "", fixed)
# RE_DROP: optional comma/colon-separated random-effect components to prune -- "fid" drops the independent
# founder intercept, an SV3 name drops that term's random garden slope. Needed for the single-technology
# panels, where some slope variances sit at the boundary (isSingular): a slope pinned at zero still costs a
# parameter and leaves that term's interaction SE anti-conservative, so the honest fix is to remove it.
# Empty by default, so every existing invocation is unchanged. Applied identically to the base and the
# drop-one-interaction refits, so the LRTs stay like-for-like.
# split on comma, colon, plus or space: sbatch --export treats commas as its own separator, so a
# multi-component RE_DROP must be passed colon-separated ("fid:sv_CDS").
RE_DROP <- setdiff(strsplit(Sys.getenv("RE_DROP", ""), "[,:+ ]+")[[1]], "")
stopifnot(all(RE_DROP %in% c("fid", SV3)))
RE_SV <- setdiff(SV3, RE_DROP)
re <- paste(c("(1 | garden)", "(1 | kin)", if (!("fid" %in% RE_DROP)) "(1 | fid)",
              sprintf("(0 + %s | garden)", RE_SV)), collapse = " + ")
# EXTRA_RE: append arbitrary random-effect terms. The motivating case is "(0 + obio1c | garden)", which
# lets LOCAL ADAPTATION have its own strength in every garden instead of being forced to scale linearly
# with garden temperature through the single obio1c:bio1c coefficient. Per-garden local-adaptation slopes
# correlate with garden temperature at only r = 0.39, so the linear term leaves most of that heterogeneity
# in the residual -- where a burden x temperature interaction can pick it up. Empty by default.
EXTRA_RE <- Sys.getenv("EXTRA_RE", "")
if (nzchar(EXTRA_RE)) re <- paste(re, EXTRA_RE, sep = " + ")
if (length(RE_DROP)) cat(sprintf("RE_DROP = %s  -> random effects: %s\n", paste(RE_DROP, collapse = ","), re))
# EXTRA_FIXED: optional extra fixed terms, e.g. "obio1c + obio1c:bio1c" to control founder origin climate
# and its interaction with garden climate (the local-adaptation confound). Empty by default.
EXTRA <- Sys.getenv("EXTRA_FIXED", "")
form <- function(inter) as.formula(paste("y ~", paste(c(fixed, if (length(inter)) sprintf("%s:bio1c", inter),
                                                        if (nzchar(EXTRA)) EXTRA),
                                                      collapse = " + "), "+", re))

SUF <- paste0(if (nzchar(EXTRA)) paste0("_x_", gsub("[^A-Za-z0-9]+", "", EXTRA)) else "",
              if (nzchar(EXTRA_RE)) paste0("_re", gsub("[^A-Za-z0-9]+", "", EXTRA_RE)) else "",
              if (!(POOL_SET %in% c("sv4", "sv6"))) paste0("_", POOL_SET) else "",
              if (nzchar(SUBSET)) paste0("_", SUBSET) else "",
              if (length(RE_DROP)) paste0("_re", gsub("[^A-Za-z0-9]+", "", paste(RE_DROP, collapse = ""))) else "",
              if (nzchar(STRATIFY)) paste0("_freq", gsub("[^A-Za-z0-9]+", "", paste(c(SCTX, FBINS), collapse = ""))) else "")

TASKS <- list(list(name = "M1p_reml", inter = character(0), REML = TRUE),
              list(name = "M2p_reml", inter = SV3, REML = TRUE),
              list(name = "M2p_ml", inter = SV3, REML = FALSE))
for (t in SV3) TASKS[[length(TASKS) + 1]] <- list(name = paste0("M2p_ml_drop_", t), inter = setdiff(SV3, t), REML = FALSE)
T <- TASKS[[task + 1]]
cat(sprintf("task %d: %s  (REML = %s)  origin = %s  set = %s  out = %s\nformula: %s\n",
            task, T$name, T$REML, ORIGIN, POOL_SET, OUT, deparse1(form(T$inter))))

t0 <- Sys.time()
m <- relmatLmer(form(T$inter), data = L, relmat = list(kin = K), REML = T$REML)
cat(sprintf("fitted in %.1f s  singular = %s\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), isSingular(m)))

s <- summary(m)$coefficients
write.csv(data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], z = s[, 3], p = 2 * pnorm(-abs(s[, 3]))),
          file.path(OUT, sprintf("coef_%s%s.csv", T$name, SUF)), row.names = FALSE)
vc <- as.data.frame(VarCorr(m))
write.csv(vc[, c("grp", "var1", "vcov")], file.path(OUT, sprintf("varcomp_%s%s.csv", T$name, SUF)), row.names = FALSE)
ll <- logLik(m)
write.csv(data.frame(model = paste0(T$name, SUF), REML = T$REML, logLik = as.numeric(ll), df = attr(ll, "df"),
                     singular = isSingular(m)), file.path(OUT, sprintf("loglik_%s%s.csv", T$name, SUF)), row.names = FALSE)

# --- self-documenting specification row -------------------------------------------------------------
# Every fit records what it was, so the model registry can never drift from the fits it describes.
# `_collect_model_registry.py` prefers these files and falls back to decoding the filename for older runs.
spec <- data.frame(
  model_id  = paste0(basename(OUT), "/", T$name, SUF),
  model     = T$name, suffix = SUF, estimator = if (T$REML) "REML" else "ML",
  formula   = deparse1(form(T$inter)),
  fixed     = fixed, random = re,
  unit = "pooled", origin = ORIGIN, te_term = POOL_SET,
  subset = if (nzchar(SUBSET)) SUBSET else "all_founders",
  extra_fixed = EXTRA, re_drop = paste(RE_DROP, collapse = ","),
  stratify = STRATIFY, stratify_ctx = paste(SCTX, collapse = ","),
  freq_bins = if (nzchar(STRATIFY)) paste(FBINS, collapse = ",") else "",
  onestage_dir = basename(IN),
  n_founders = length(unique(L$founder)), n_gardens = length(unique(L$garden)), n_rows = nrow(L),
  logLik = as.numeric(ll), df = attr(ll, "df"), singular = isSingular(m),
  fitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  git_commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
                        error = function(e) NA_character_),
  stringsAsFactors = FALSE)
write.csv(spec, file.path(OUT, sprintf("spec_%s%s.csv", T$name, SUF)), row.names = FALSE)
if (T$REML) {
  write.csv(as.matrix(vcov(m)), file.path(OUT, sprintf("vcov_%s%s.csv", T$name, SUF)))
  # per-garden total effect of each pooled SV term (fixed part + random slope), as in _sv_burden_final_models.R
  fe <- fixef(m); rg <- ranef(m, condVar = TRUE)$garden; pv <- attr(rg, "postVar"); out <- NULL
  # every random garden slope, not just the SV ones: EXTRA_RE can add others (e.g. obio1c), and their
  # per-garden values are exactly what a "how strong is local adaptation in this garden" figure needs
  for (t in setdiff(names(rg), "(Intercept)")) {
    k <- match(t, names(rg)); sdv <- rep(NA_real_, nrow(rg))
    if (is.list(pv) && length(pv) >= k) sdv <- sqrt(pv[[k]][1, 1, ]) else if (is.array(pv)) sdv <- sqrt(pv[k, k, ])
    x <- b1$bio1c[match(rownames(rg), b1$garden)]
    # lme4 may name the interaction either way round depending on model-frame order
    inm <- intersect(c(paste0(t, ":bio1c"), paste0("bio1c:", t)), names(fe))
    fixed_part <- fe[t] + (if (length(inm)) fe[inm[1]] * x else 0)
    out <- rbind(out, data.frame(term = t, garden = rownames(rg), bio1 = b1$bio1[match(rownames(rg), b1$garden)],
                                 fixed_part = fixed_part, slope_blup = rg[[t]], cond_sd = sdv, total = fixed_part + rg[[t]]))
  }
  write.csv(out, file.path(OUT, sprintf("garden_%s%s.csv", T$name, SUF)), row.names = FALSE)
  # per-observation partial residuals, so a per-garden "fitness against burden" figure can be drawn from
  # THIS fit instead of from raw per-garden regressions. The component-plus-residual for term t is
  #   resid + (that garden's TOTAL slope for t) * x_t
  # which has every other covariate -- kinship, the SNP burden, assembly route, origin climate, local
  # adaptation, the other burdens -- already removed, and whose within-garden slope is exactly the
  # model's conditional effect in that garden.
  rs <- data.frame(founder = as.character(L$founder), garden = as.character(L$garden),
                   y = L$y, fitted = fitted(m), resid = residuals(m))
  for (t in setdiff(names(rg), "(Intercept)")) if (t %in% names(L)) rs[[t]] <- L[[t]]
  write.csv(rs, file.path(OUT, sprintf("resid_%s%s.csv", T$name, SUF)), row.names = FALSE)
}
cat("done\n")
