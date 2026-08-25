# 06_sv_support_filter (archived 2026-08-25)

"How many SVs do we lose if we filter support at *x*?" — a pre-GEA decision notebook on the
kMate support stats (`n_called`, `info`, `se`).

**Archived because its recommendation is superseded.** It argued for `n_called ≥ 100`
(~78% of SVs kept). Production is stricter and set elsewhere: `lib.founder_panel_keep`
defaults to `min_mac=COMMON_MAC` (MAF≥0.05) **and `called_min=0.9`** — i.e. ≥208 of 231
founders, the same bar that defines the LD blocks and the SV landscape. Nothing downstream
reads the `n_called ≥ 100` cut this notebook proposed.

## The one thing worth carrying forward

> `info` is per-pool support (`h·V_called`, the projection denominator); `n_called` is
> panel-level. They correlate **0.95**. **If low support correlates with coverage — and
> coverage tracks site — then unmasked low-info cells can manufacture a climate
> association.**

That warning bears directly on the SV climate-gradient results, so it is worth knowing it
was considered. It is largely defused in practice by the production `called_min=0.9` floor
(high-support records only), but if any future analysis relaxes that floor, re-read this.

Also noted there and still open: the rare-SV identical-trajectory artifact (Chr3 clusters) is
*well supported* (n_called≈191), so a support cut does not address it — it needs a MAF filter
+ haploblock collapse.

## Contents

    06_sv_support_filter.ipynb
    scripts/_build_support_nb.py
