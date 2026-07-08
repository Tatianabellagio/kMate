#!/usr/bin/env python
"""AUDIT of the passenger claim, in the GWAS's OWN currency.

The passenger scripts used a RAW founder-fitness gap (carrier vs non-carrier mean of the
site-averaged selection), not the kinship-corrected per-site EMMAX z the founder-GWAS actually
used to flag these blocks. That can mislead: raw gaps are confounded by relatedness, and the
site-averaged 'global gap' is ~0 for a heterogeneous (JOINT) hit even if causal.

This re-runs, for the SAME top-0.5% SV blocks, the IDENTICAL pipeline the GWAS used
(per-site QN selection -> LOCO-EMMAX kinship-corrected z -> per-site lambda GC -> Bolormaa
JOINT with the saved C) on:
  (i)  each block's LEAD hap-cluster  -> self-check: my JOINT must reproduce the npz chi2_joint
  (ii) the SVs sitting in those blocks -> the real question
plus founder-level r^2 between each SV's carrier vector and its block's lead hap-cluster.

Verdict logic: if SV JOINT is null-level (lambda~1, JOINT~S) AND r^2 low while the lead-hap
JOINT is high -> passenger CONFIRMED in the GWAS currency. If SV JOINT is comparable / r^2 high
-> the raw-gap passenger test misled. Env: kmate.
"""
import os, sys, glob
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype
from ecotype_selection_site import genome_h
import founder_gwas_multisite as FG          # build_trait_site, loco_emmax, lamgc

WIN = lib.OUT; SEED = lib.SEEDMIX
PANEL = "panel/arch3"; CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
MAC_MIN, MAC_GRM = 3, 12

raw = lib.multisite_gwas_raw("clq90_pc1")
sites, C, Cinv, S = raw["sites"], raw["C"], raw["Cinv"], raw["S"]
lam_site = np.array(raw["meta"]["per_site_lambda"])         # per-site GC used by the GWAS
chi2_npz = dict(zip(raw["unit"], raw["chi2_joint"]))        # block JOINT from the GWAS (max over its haps)

# ---- founder genotype + same poly/grm filters as the GWAS; align lead hap to Z rows ----
os.environ["MEMB_TAG"] = "clq90"
G, founders, reg = build_genotype("clq90"); nF = len(founders)
cnt = G.sum(0); blk = reg.block.to_numpy()
order = np.lexsort((-cnt, blk)); sbk = blk[order]
is_ref = np.zeros(len(cnt), bool); is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
grm_set = (~is_ref) & (cnt >= MAC_GRM) & (cnt <= nF - MAC_GRM)
Gp = G[:, poly].astype(np.float64); regp = reg[poly].reset_index(drop=True)
Gg = G[:, grm_set].astype(np.float64); chrom_g = reg.chrom.to_numpy()[grm_set]
regp_unit = (regp.chrom + ":" + regp.start.astype(str) + "-" + regp.end.astype(str)).to_numpy()
assert np.array_equal(regp_unit, raw["unit"]), "marker order drift vs GWAS npz"

# ---- per-site QN trait, exactly as the GWAS (cached window-mode h == GWAS's h source) ----
cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
H = cache["H"]; csamp = {s: i for i, s in enumerate(cache["samples"].astype(str))}
assert list(cache["founders"].astype(str)) == list(founders), "founder order mismatch"
gh = lambda s: (H[csamp[s]] if s in csamp else None)
seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
# spot-check the cache equals live WIN genome_h for a couple samples (trait-source audit)
pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(csamp)]
chk = [s for s in pt.sampleid.astype(str).unique()[:3]]
for s in chk:
    live = genome_h(s, WIN)
    if live is not None:
        assert np.allclose(H[csamp[s]], live, atol=1e-6), f"cache != WIN for {s}"
print(f"trait-source check OK (cache H == window-mode genome_h) on {len(chk)} samples")

yq_by_site = {}
for site, sd in pt.groupby("site"):
    if int(site) not in set(int(x) for x in sites):
        continue
    s_f, present = FG.build_trait_site(sd, p0, gh, nF)
    if s_f is None:
        continue
    yq_by_site[int(site)] = stats.norm.ppf((stats.rankdata(s_f) - 0.5) / nF)
order_sites = [int(s) for s in sites]
assert all(s in yq_by_site for s in order_sites), "missing a GWAS site in trait rebuild"

def joint_of(Gtest, chrom_test):
    """(n x nsite) GC-scaled z then JOINT via saved C, for test genotypes Gtest (nF x n)."""
    Zt = np.empty((Gtest.shape[1], len(order_sites)))
    for k, site in enumerate(order_sites):
        z, _ = FG.loco_emmax(yq_by_site[site], Gtest, np.asarray(chrom_test), Gg, chrom_g, nF)
        Zt[:, k] = z / np.sqrt(max(lam_site[k], 1e-9))
    Jt = np.einsum("mi,ij,mj->m", Zt, Cinv, Zt)
    return Jt, Zt

# ---- top-0.5% JOINT SV blocks (match the passenger test) ----
up = pd.read_csv(f"{FG.H}/multisite_founder_gwas_clq90_pc1.csv").groupby("unit", as_index=False).p_joint.min()
L = pd.read_csv(f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv").rename(columns={"block_id": "unit"})
up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
up = up[up.n_kept >= 2]
ntop = int(round(0.005 * len(up)))
topsv = up.nsmallest(ntop, "p_joint"); topsv = topsv[topsv.has_sv == 1].reset_index(drop=True)
units = list(topsv.unit)
print(f"top-0.5% JOINT SV blocks: {len(units)}")

# lead hap-cluster (GWAS unit) per block: argmax chi2 marker -> its Gp column
lead_col, lead_chrom, lead_unit = [], [], []
for u in units:
    rows = np.where(regp_unit == u)[0]
    lr = int(rows[np.argmax(raw["chi2_joint"][rows])])
    lead_col.append(lr); lead_chrom.append(regp.chrom.iloc[lr]); lead_unit.append(u)
Ghap = Gp[:, lead_col]                                      # (nF x nblock) lead hap membership
Jhap, _ = joint_of(Ghap, lead_chrom)
# self-check: reproduce npz chi2_joint for the lead haps
npz_J = np.array([chi2_npz[u] for u in lead_unit])
rc = np.corrcoef(Jhap, npz_J)[0, 1]
print(f"SELF-CHECK lead-hap JOINT vs npz chi2_joint: r={rc:.3f}, "
      f"median mine {np.median(Jhap):.1f} vs npz {np.median(npz_J):.1f} "
      f"({'PIPELINE OK' if rc > 0.9 else 'PIPELINE MISMATCH -- audit invalid'})")

# ---- SV genotypes in those blocks + r^2 to the block's lead hap ----
sv_g, sv_chrom, sv_unit, sv_pos, sv_r2 = [], [], [], [], []
unit_lead = dict(zip(lead_unit, lead_col))
for ch in CHROMS:
    cl = ch.lower(); meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64); dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
    vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
    vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
    na = np.asarray(vp.sum(0)).ravel(); ncl = np.asarray(vc.sum(0)).ravel()
    svmask = (dl > 50) & (na >= 12) & (na <= nF - 12) & (ncl / nF >= 0.9)
    blk = topsv[topsv.chrom == ch]
    for j in np.where(svmask)[0]:
        hit = blk[(blk.start_pos <= pos[j]) & (blk.end_pos >= pos[j])]
        if not len(hit):
            continue
        u = hit.iloc[0].unit; car = vp[:, j].toarray().ravel().astype(float)
        if car.sum() < 2 or car.sum() > nF - 2:
            continue
        sv_g.append(car); sv_chrom.append(ch); sv_unit.append(u); sv_pos.append(int(pos[j]))
        lh = Gp[:, unit_lead[u]]
        sv_r2.append(np.corrcoef(car, lh)[0, 1] ** 2)
SVG = np.array(sv_g).T                                      # (nF x n_sv)
Jsv, _ = joint_of(SVG, sv_chrom)
p_sv = stats.chi2.sf(Jsv, len(order_sites)); p_hap = stats.chi2.sf(Jhap, len(order_sites))

res = pd.DataFrame(dict(unit=sv_unit, sv_pos=sv_pos, r2_sv_leadhap=np.round(sv_r2, 2),
                        sv_JOINT=np.round(Jsv, 1), sv_p=p_sv.round(4)))
res["block_hap_JOINT"] = [round(float(np.array(Jhap)[lead_unit.index(u)]), 1) for u in sv_unit]
res = res.sort_values("sv_JOINT", ascending=False)
pd.set_option("display.width", 170, "display.max_rows", 60)
print(f"\n=== SV kinship-corrected JOINT (same pipeline as GWAS) vs their block's lead-hap JOINT ===")
print(res.to_string(index=False))
res.to_csv(f"{lib.GEA}/sv_adaptive/audit_sv_passenger_emmax.csv", index=False)

print(f"\nSUMMARY ({len(Jsv)} SVs in {len(units)} top-0.5% blocks):")
print(f"  SV JOINT: median {np.median(Jsv):.1f} (null chi2_{len(order_sites)} median {stats.chi2.ppf(0.5,len(order_sites)):.1f}), "
      f"lambda_GC {lib.lamgc(p_sv):.2f} | genome-wide FDR<.05: {int((lib.bh(p_sv)<0.05).sum())}")
print(f"  lead-hap JOINT: median {np.median(Jhap):.1f} (these are top-0.5% hits by construction)")
print(f"  founder r^2(SV, lead hap): median {np.median(sv_r2):.2f} | r2>=0.9 in {int((np.array(sv_r2)>=0.9).sum())}/{len(sv_r2)} | r2<0.5 in {int((np.array(sv_r2)<0.5).sum())}/{len(sv_r2)}")
share = Jsv / np.array([float(np.array(Jhap)[lead_unit.index(u)]) for u in sv_unit])
print(f"  SV JOINT / block-hap JOINT: median {np.median(share):.2f}  (1.0 = SV as strong as the selected haplotype; ~0 = passenger)")
