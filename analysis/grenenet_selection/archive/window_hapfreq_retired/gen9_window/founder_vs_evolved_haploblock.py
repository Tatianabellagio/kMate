#!/usr/bin/env python
"""Founder vs evolved PC1-VE at the HAPLOBLOCK level (haplotype cluster), one chrom.
For each kept haplotype cluster (founding freq>=2/231) within a dynld unit, take its
SIGNATURE variants (alt in >=50% of cluster founders AND >=0.5 higher than outside) and
compute PC1-VE of (a) the founder genotypes and (b) the gen9 pool AF on those variants,
both global-mode and window-mode. The block-level analog is founder_vs_evolved_dynld.py.

Usage: founder_vs_evolved_haploblock.py Chr1   (kmate env, >=48G)
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
SIG_PRESENT, SIG_DIFF, MINFREQ = 0.5, 0.5, 2 / 231
CLS = ["snp", "sv", "smallindel"]

def pc1_ve(M):
    M = np.asarray(M, np.float64)
    if M.shape[1] < 2:
        return np.nan, np.nan
    cm = np.nanmean(M, 0); ind = np.where(np.isnan(M)); M = M.copy(); M[ind] = np.take(cm, ind[1])
    M = M - M.mean(0); sd = M.std(0); sd[sd == 0] = 1.0; M /= sd
    ev = np.linalg.svd(M, compute_uv=False) ** 2; tot = ev.sum()
    return (np.nan, np.nan) if tot == 0 else (float(ev[0] / tot), float(tot ** 2 / (ev ** 2).sum()))

chrom = sys.argv[1]; chrlc = chrom.lower(); t0 = time.time()
# evolved AF (global + window), same record set/order, per chrom
Gp, gpos = [], []
for c in CLS:
    rc = pd.read_csv(f"{CM}/{c}_gen9.records.csv", usecols=["chrom", "pos"])
    cm = (rc.chrom == chrom).to_numpy()
    Gp.append(np.load(f"{CM}/{c}_gen9_af.npy", mmap_mode="r")[:, cm]); gpos.append(rc.pos.to_numpy()[cm])
Gc = np.hstack(Gp); gpos = np.concatenate(gpos)
wrec = pd.read_csv(f"{GW}/records.csv"); wm = (wrec.chrom == chrom).to_numpy()
Wc = np.asarray(np.load(f"{GW}/window_gen9_af.npy", mmap_mode="r")[:, wm]); wpos = wrec.pos.to_numpy()[wm]
og = np.argsort(gpos, kind="stable"); Gc = np.ascontiguousarray(Gc[:, og]); gpos = gpos[og]
ow = np.argsort(wpos, kind="stable"); Wc = np.ascontiguousarray(Wc[:, ow]); wpos = wpos[ow]
assert bool((gpos == wpos).all())
print(f"{chrom}: evolved {Gc.shape} ({time.time()-t0:.0f}s)", flush=True)

_, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
positions = np.asarray(positions); geno = (raw >= 0.5).astype(np.int8)
lab = np.load(f"{GW}/clusters/{chrlc}_labels.npz"); labels = lab["labels"]; st = lab["start"]; en = lab["end"]; neff = lab["n_eff"]
print(f"{chrom}: founder {geno.shape} ({time.time()-t0:.0f}s)", flush=True)

def evolved_cols(sigpos):
    keep = []
    for p in sigpos:
        j = np.searchsorted(gpos, p)
        if j < len(gpos) and gpos[j] == p:
            keep.append(j)
    return np.array(keep, int)

rows = []
for ui in range(len(st)):
    lo = int(np.searchsorted(positions, st[ui])); hi = int(np.searchsorted(positions, en[ui], side="right"))
    sub = geno[:, lo:hi]; bpos = positions[lo:hi]
    if sub.shape[1] < 2:
        continue
    lab_u = labels[ui]
    for c in np.unique(lab_u):
        if c < 0:
            continue
        inm = lab_u == c; n_in = int(inm.sum())
        if n_in / 231.0 < MINFREQ:
            continue
        mean_in = sub[inm].mean(0); mean_out = sub[~inm].mean(0) if n_in < sub.shape[0] else np.zeros(sub.shape[1])
        sig = (mean_in >= SIG_PRESENT) & (mean_in - mean_out >= SIG_DIFF)
        if sig.sum() < 2:
            continue
        fve, fed = pc1_ve(sub[:, sig])
        ec = evolved_cols(bpos[sig])
        if len(ec) >= 2:
            gve, ged = pc1_ve(Gc[:, ec]); wve, wed = pc1_ve(Wc[:, ec])
        else:
            gve = wve = ged = wed = np.nan
        rows.append((chrom, ui, int(st[ui]), int(en[ui]), int(c), round(float(neff[ui]), 3),
                     n_in, int(sig.sum()), len(ec), fve, gve, wve, fed, ged, wed))

df = pd.DataFrame(rows, columns=["chrom", "unit_idx", "start", "end", "cluster", "n_eff",
                                 "n_founders", "n_sig", "n_sig_matched", "founder_ve",
                                 "evolved_ve_global", "evolved_ve_window", "eff_dim_founder",
                                 "eff_dim_global", "eff_dim_window"]).round(4)
out = f"{GW}/haploblock_founder_vs_evolved_{chrlc}.csv"
df.to_csv(out, index=False)
print(f"{chrom}: {len(df):,} haploblocks -> {out} ({time.time()-t0:.0f}s)", flush=True)
