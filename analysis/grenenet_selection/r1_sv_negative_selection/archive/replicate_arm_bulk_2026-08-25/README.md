# Replicate arm — bulk parallelism + PicMin (archived 2026-08-25)

Sections 1 and 2 of the former `parallelism_picmin.ipynb`. Archived under the
paper-streamlining pass: neither supports nor threatens the headline, so neither earns a
notebook. **The one-liners are the result** — keep them, don't re-run this.

## §1 — bulk parallelism: collapsed

> On post-Kf_w data the SV parallel-responder enrichment is small and frequency-limited:
> responder **1.03×**, repeatable ≥⅓ sites **1.19×**, strongly-repeatable ≥½ **1.30×**, and
> significant only in the **rare** MAF band (mid p=0.085, common p=0.194, both n.s.).
> Insertions carry what little there is (deletions ≈ SNP).

The pre-Kf_w headline — "escalates 1.2→1.7→3.6×, holds at all MAF" — **does not reproduce**.
That earlier number was frequency / founder-projection confounded.

## §2 — PicMin: near parity, and structurally unable to corroborate

> SVs are **not** meaningfully over-represented among repeated-adaptation loci across the 31
> site-lineages: fold **~1.07–1.12×** (q<0.1: 17.4% vs matched 15.8%; q<0.05: 14.9% vs
> ~13.3–13.9%). Nominally insertion-driven (ins ≈17% vs del ≈12% at q<0.05).

Two things worth carrying forward:

1. **PicMin is direction-agnostic.** Its per-locus input is drift-controlled parallelism
   `|z|` — an absolute value — so a variant purged hard everywhere and one sweeping up
   everywhere score identically. It measures *repeatability of response*, not direction, and
   therefore could not have corroborated a purging result even if it had come out strong.
   This is why its near-parity is **not** evidence against the negative-selection finding: a
   signal confined to a hot-site tail is exactly what a 31-lineage repeatability test dilutes.
2. **Near parity was never formally tested.** The post-Kf_w fold was not re-run through a
   significance test, so this is *near parity*, not *confirmed null*.

The absolute significant fraction is high (13–20%) because this founder-projection system has
pervasive parallel sorting; only the SV-vs-SNP **relative** contrast is interpretable.

## What was kept instead

§3 of the same notebook — the climate gradient of the parallelism excess — survived the Kf_w
fix and corroborates the headline via an independent statistic. It now lives on its own as
`notebooks/sv_parallelism_climate.ipynb`, built by
`../../_build_sv_parallelism_climate_nb.py`.

## Contents

    parallelism_picmin.ipynb        the original 3-section notebook (executed)
    scripts/_build_parallelism_picmin_nb.py   its builder (all 3 sections)
    scripts/_picmin.py              the PicMin implementation (Booker et al.)
    results/picmin.npz              PicMin output

`parallelism.npz` was **not** archived — it stays in `results/sv_adaptive/` because the live
§3 notebook reads it.
