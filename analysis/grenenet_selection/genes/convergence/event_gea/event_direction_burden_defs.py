#!/usr/bin/env python
"""Event-level GEA direction under the BURDEN model's own definitions.

The allele- and event-level GEA both show promoter SV deletions leaning slightly COLD relative
to intergenic ones, where the burden model (notebooks/sv_climate_maf.ipynb) finds promoter
deletions favoured WARM. Two definitional differences remain after the unit is matched:
  context  burden: CDS > UTR > intron > reference TE > promoter (<= 2 kb) > intergenic;
           the GEA classifier: gene body > promoter (<= 1 kb) > TE > intergenic
  filter   burden: `shared` variants only (carried by cold- AND warm-origin founders), called
           in >= 90% of founders with >= 2 carriers (`passes`)
sv_event_catalog.csv carries exactly those labels, row-aligned with the event units, so this
applies them to the event-level GEA's signed z unchanged. Same contrast and CI as
gea_direction.contrast_mean_z. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GEA = os.path.abspath(os.path.join(CONV, "..", ".."))
sys.path.insert(0, CONV)
import gea_region_by_axis as GR, gea_direction as GD           # noqa: E402

U = np.load(f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz")
sv_units = np.where(np.abs(U["dlen"]) > 50)[0]
C = pd.read_csv(f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_event_catalog.csv")
assert len(C) == len(sv_units)
E = pd.read_csv(f"{HERE}/sv_gen9_units.csv")
E["cat"] = E.unit.map(pd.Series(np.arange(len(sv_units)), index=sv_units))
c = C.iloc[E.cat.values].reset_index(drop=True)
E["context"], E["origin"], E["passes"] = c.context.values, c.origin.values, c.passes.values
E["kind"] = np.where(E.dlen < 0, "deletion", "insertion"); E["size"] = E.dlen.abs()
E["grp"] = E.context.where(E.context.isin(["promoter", "intergenic"]), "other")
E["stratum"] = (pd.cut(E["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                + pd.cut(E.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
rng = np.random.default_rng(1); rows = []
for ax in GD.TEMP:
    T = pd.read_csv(f"{HERE}/lfmm_event_sv_gen9_{ax}.csv.gz", usecols=["unit", "z"])
    E["z"] = E.unit.map(T.set_index("unit").z)
    for kind in ("deletion", "insertion"):
        for filt, m in [("all tested", E.kind == kind),
                        ("shared + passes (burden filter)", (E.kind == kind) & (E.origin == "shared") & E.passes.astype(bool))]:
            d = E[m]
            e, lo, hi, p = GD.contrast_mean_z(d, "promoter", "intergenic", rng)
            rows.append(dict(axis=ax, kind=kind, set=filt, n_promoter=int((d.grp == "promoter").sum()),
                             n_intergenic=int((d.grp == "intergenic").sum()), dz=e, lo=lo, hi=hi, p=p))
R = pd.DataFrame(rows)
R.to_csv(f"{HERE}/event_direction_burden_defs.csv", index=False)
pd.set_option("display.width", 200)
print("event-level GEA, burden context (promoter <= 2 kb, TE first); dz = mean z promoter - intergenic; z > 0 = warm")
print(R.round(3).to_string(index=False))
