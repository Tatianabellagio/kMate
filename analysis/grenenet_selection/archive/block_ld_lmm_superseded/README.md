# Retired: window-mode haploblock LD-LMM selection family

Archived 2026-07-07. This is a whole family of scripts testing **per-haploblock TEMPORAL
selection** (founding → gen3, default site 4) with a kinship/LD control — NOT a climate-GEA
test. `block_ld_lmm.py`/`_temporal.py` fit `s = mu*1 + g + e` where `g ~ N(0, sg^2 C_LD)` is
the founder-relationship (LD/haplotype-sharing) background and the residual `r_b = s_b - g_b`
is the claimed block-specific selection; `_linear.py` reruns it on the raw-frequency scale to
avoid logit-boundary artifacts; `_binom.py`/`_driftnull.py`/`_oneside.py`/`_sampvar.py` are
successive variance/null refinements (honest binomial sampling variance, a generative
selfing-projection null, one-sided adaptive-direction test).

**Why retired:** every script here reads from `results/grenenet_gea/hapfreq` — the
**window-mode** per-block founder-frequency store. `GLOBAL_MODE_DECISION.md` (2026-07-01)
decided evolved allele frequencies should be estimated in GLOBAL mode, not window/block mode,
for this system (~97% selfing + ~3 generations means window-mode recombination detection is
circular and underpowered — see that doc). `docs/RERUN_AFTER_FIX.md` (§ window-mode founder-GWAS
family) explicitly lists the whole `block_ld_lmm*` family, alongside `hapfreq/`,
`founder_gwas_multisite.py`, and `sv_adaptive/`, as **stale and not to be rerun** unless window
mode is revived — which it has not been; the window stores these scripts depend on have since
been deleted from disk (`archive/window_hapfreq_retired/`).

**Superseded by:** nothing 1:1 — the *question* (temporal, kinship-controlled, per-block
selection) doesn't have a live replacement on GLOBAL-mode AF; `BAYPASS_TEMPORAL_PLAN.md`
sketches the intended GLOBAL-mode successor (BayPass C2 contrast, Omega-conditioned) but is an
unexecuted plan, not a finished pipeline. If this thread is ever revived, it needs a window-mode
cohort rerun (per_founder default) first, per `docs/RERUN_AFTER_FIX.md`.
