#!/usr/bin/env python
"""How many variants actually enter the burden model, and what the filters remove.

The burden terms are per-founder RATES (carried / called), so the natural question at a talk
-- "how many insertions and deletions are there in each context?" -- is not answerable from the
model output. This counts them, at the same unit and under the same filters as
`_sv_burden_build.py` (SV_UNIT=event), from `sv_event_catalog.csv`.

  a  the filter cascade: every truvari-collapsed SV event -> called in >= 90% of the 231
     founders -> >= 2 carrier founders (together: the catalog's `passes`) -> `shared`, i.e.
     carried by at least one cold-origin AND one warm-origin founder.
  b  what survives, by context x kind. Context is first-match CDS > UTR > intron > reference
     TE > promoter (<= 2 kb upstream) > intergenic. Filled bars are the three contexts that
     enter the reported model as terms; open bars are classified but not modelled (SV4).

No titles (house rule); panels carry an in-panel corner annotation. Writes the figure next to
the other main-figure plots and prints the table.

env: basic (plotting). Usage: python _sv_burden_variant_counts.py
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
G = f"{HERE}/results/sv_adaptive"
NB = os.path.join(os.path.dirname(HERE), "notebooks", "plots")
os.makedirs(NB, exist_ok=True)
CTXO = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
MODELLED = {"CDS", "promoter", "intergenic"}                      # SV4 terms
PAL = {"ins": {"CDS": "#2F4A1F", "UTR": "#375325", "intron": "#3E5B2C", "TE": "#466333",
               "promoter": "#4E6B3A", "intergenic": "#9DB58C"},
       "del": {"CDS": "#B04A00", "UTR": "#C85503", "intron": "#DF5F08", "TE": "#F26A10",
               "promoter": "#F07A2E", "intergenic": "#F8BE91"}}
KL = {"ins": "insertions", "del": "deletions"}

D = pd.read_csv(f"{G}/sv_event_catalog.csv")
call = D.n_called >= 0.9 * 231
mac = D.n_carriers >= 2
shared = D.origin == "shared"
STAGES = [("all SV events", np.ones(len(D), bool)),
          ("called in ≥90%\nof 231 founders", call),
          ("+ ≥2 carrier\nfounders", call & mac),
          ("+ shared\n(cold and warm)", call & mac & shared)]
assert (D.passes == (call & mac)).all(), "catalog `passes` is not call-rate & MAC"
B = D[call & mac & shared]

fig = plt.figure(figsize=(11.6, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.18], wspace=0.42)


def style(a):
    a.set_axisbelow(True); a.grid(color="0.9", lw=0.6, axis="x")
    for s in a.spines.values(): s.set_visible(False)


# a -- the filter cascade
a = fig.add_subplot(gs[0, 0]); style(a); h = 0.36
for j, k in enumerate(("ins", "del")):
    y = np.arange(len(STAGES)) + (h / 2 if k == "ins" else -h / 2)
    n = [int((m & (D.kind == k)).sum()) for _, m in STAGES]
    a.barh(y, n, height=h, color=PAL[k]["CDS"], label=KL[k])
    for yy, v in zip(y, n):
        a.annotate(f"{v:,}", xy=(v, yy), xytext=(4, 0), textcoords="offset points",
                   va="center", fontsize=7.2)
a.set_yticks(range(len(STAGES))); a.set_yticklabels([s for s, _ in STAGES], fontsize=8)
a.invert_yaxis(); a.set_xlabel("SV events", fontsize=8.5); a.set_xlim(0, 58000)
a.legend(frameon=False, fontsize=8, loc="lower right")
a.annotate("a", xy=(-0.42, 1.04), xycoords="axes fraction", fontsize=13)

# b -- what enters the model, by context
a = fig.add_subplot(gs[0, 1]); style(a)
for j, k in enumerate(("ins", "del")):
    y = np.arange(len(CTXO)) + (h / 2 if k == "ins" else -h / 2)
    n = [int(((B.context == c) & (B.kind == k)).sum()) for c in CTXO]
    for yy, c, v in zip(y, CTXO, n):
        inmod = c in MODELLED
        a.barh(yy, v, height=h, color=PAL[k][c] if inmod else "white",
               edgecolor=PAL[k][c], lw=1.1 if not inmod else 0)
        a.annotate(f"{v:,}", xy=(v, yy), xytext=(4, 0), textcoords="offset points",
                   va="center", fontsize=7.2, color="0.25" if not inmod else "black")
a.set_yticks(range(len(CTXO)))
a.set_yticklabels([c + ("" if c in MODELLED else "  (not a model term)") for c in CTXO], fontsize=8)
a.invert_yaxis(); a.set_xlabel("SV events entering the burden model", fontsize=8.5)
a.set_xlim(0, 3650)
a.annotate("b", xy=(-0.30, 1.04), xycoords="axes fraction", fontsize=13)
a.annotate(f"{len(B):,} events total\nfilled = modelled term", xy=(0.97, 0.05),
           xycoords="axes fraction", ha="right", va="bottom", fontsize=7.5, color="0.35")

for ext in ("png", "pdf"):
    fig.savefig(f"{NB}/sv_burden_variant_counts.{ext}", dpi=200, bbox_inches="tight")
print(f"[written] {NB}/sv_burden_variant_counts.png")

print("\nfilter cascade (SV events):")
for s, m in STAGES:
    print(f"  {s.replace(chr(10),' '):32s} {int(m.sum()):>7,}   "
          f"ins {int((m&(D.kind=='ins')).sum()):>6,}  del {int((m&(D.kind=='del')).sum()):>6,}")
T = (B.pivot_table(index="context", columns="kind", values="key", aggfunc="count")
       .reindex(CTXO).fillna(0).astype(int))
T["total"] = T.sum(axis=1); T["modelled"] = [c in MODELLED for c in T.index]
print("\nevents entering the burden model, by context:"); print(T.to_string())
print(f"\nTOTAL {int(T.total.sum()):,}   modelled contexts only "
      f"{int(T.loc[list(MODELLED),'total'].sum()):,}")
