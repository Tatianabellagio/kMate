#!/usr/bin/env python
"""Per dynld-UNIT founder vs evolved PC1-VE for ONE chrom, evolved BOTH ways
(global-mode AF vs block-based window-mode h) on identical gen9 records + 355 pools.
Per-chrom so each job loads only that chrom's columns (build_common_matrix is ~12GB).

Usage: founder_vs_evolved_dynld.py Chr1   (run in kmate env, >=48G)
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from recompute_blocks import build_common_matrix

GW = f"{lib.GEA}/gen9_window"
CM = f"{lib.GEA}/phase1_replication/class_matrices"
MD = f"{lib.GEA}/blocks_mcf90"
MAF, MINCF = 0.05, 0.9
CLS = ["snp", "sv", "smallindel"]

def pc1_ve(M):
    M = np.asarray(M, np.float64)
    if M.shape[1] < 2:
        return np.nan, np.nan
    cm = np.nanmean(M, 0); ind = np.where(np.isnan(M)); M = M.copy(); M[ind] = np.take(cm, ind[1])
    M = M - M.mean(0); sd = M.std(0); sd[sd == 0] = 1.0; M /= sd
    ev = np.linalg.svd(M, compute_uv=False) ** 2
    tot = ev.sum()
    return (np.nan, np.nan) if tot == 0 else (float(ev[0] / tot), float(tot ** 2 / (ev ** 2).sum()))

chrom = sys.argv[1]; chrlc = chrom.lower(); t0 = time.time()
# global gen9 AF for this chrom (load only chrom cols of each class)
Gp, gpos = [], []
for c in CLS:
    rc = pd.read_csv(f"{CM}/{c}_gen9.records.csv", usecols=["chrom", "pos"])
    cm = (rc.chrom == chrom).to_numpy()
    Gp.append(np.load(f"{CM}/{c}_gen9_af.npy", mmap_mode="r")[:, cm]); gpos.append(rc.pos.to_numpy()[cm])
Gc = np.hstack(Gp); gpos = np.concatenate(gpos)
# window gen9 AF for this chrom (same record set/order via records.csv)
wrec = pd.read_csv(f"{GW}/records.csv")
wm = (wrec.chrom == chrom).to_numpy()
Wc = np.asarray(np.load(f"{GW}/window_gen9_af.npy", mmap_mode="r")[:, wm]); wpos = wrec.pos.to_numpy()[wm]
og = np.argsort(gpos, kind="stable"); Gc = np.ascontiguousarray(Gc[:, og]); gpos = gpos[og]
ow = np.argsort(wpos, kind="stable"); Wc = np.ascontiguousarray(Wc[:, ow]); wpos = wpos[ow]
assert gpos.shape == wpos.shape and bool((gpos == wpos).all()), "global/window record mismatch"
print(f"{chrom}: G,W {Gc.shape} loaded ({time.time()-t0:.0f}s)", flush=True)

_, fgeno, fpos = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
fpos = np.asarray(fpos)
print(f"{chrom}: founder {fgeno.shape} ({time.time()-t0:.0f}s)", flush=True)

units = pd.read_csv(f"{MD}/{chrlc}_units_dynld_K500.tsv", sep="\t")
us, ue = units.start_pos.to_numpy(), units.end_pos.to_numpy()
nv, cov, pk = units.n_variants.to_numpy(), units.covered.to_numpy(), units.panel_kmers.to_numpy()
rows = []
for k in range(len(units)):
    s, e = int(us[k]), int(ue[k])
    flo, fhi = int(np.searchsorted(fpos, s)), int(np.searchsorted(fpos, e, side="right"))
    fve = pc1_ve(fgeno[:, flo:fhi])[0] if fhi - flo >= 2 else np.nan
    elo, ehi = int(np.searchsorted(gpos, s)), int(np.searchsorted(gpos, e, side="right"))
    nev = ehi - elo
    if nev >= 2:
        gve, ged = pc1_ve(Gc[:, elo:ehi]); wve, wed = pc1_ve(Wc[:, elo:ehi])
    else:
        gve = wve = ged = wed = np.nan
    rows.append((chrom, s, e, int(nv[k]), nev, bool(cov[k]), int(pk[k]), fve, gve, wve, ged, wed))

df = pd.DataFrame(rows, columns=["chrom", "start_pos", "end_pos", "n_variants", "n_evolved",
                                 "covered", "panel_kmers", "founder_ve", "evolved_ve_global",
                                 "evolved_ve_window", "eff_dim_global", "eff_dim_window"]).round(4)
out = f"{GW}/founder_vs_evolved_{chrlc}.csv"
df.to_csv(out, index=False)
print(f"{chrom}: wrote {len(df):,} units -> {out} ({time.time()-t0:.0f}s)", flush=True)
