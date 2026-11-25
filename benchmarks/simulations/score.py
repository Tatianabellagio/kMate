#!/usr/bin/env python
"""Score every fit in runs.py against its pool's truth; write the tables the benchmark notebook draws
(analysis/grenenet_selection/notebooks/_build_kmate_benchmarks_nb.py) to config.RESULTS.

Allele frequencies: truth (pools/<pool>/truth_<panel>.tsv.gz) and estimate (fits/<fit>/fit.tsv) list the variant panel's
records in the same order, so they are paired row by row (checked against the panel; a key join is impossible because
multiallelic sites repeat (chrom, pos, ref_len, alt_len)). Scored records pass the project variant filter: called in at
least 50% of the founders and minor-allele count at least 2 among the called ones. Classes: SNP (ref and alt 1 bp),
short (longer allele under 50 bp), SV (50 bp or more). Founder frequencies: fits/<fit>/fit.h_per_chrom.npz against the
pool weights. Time: fits/<fit>/time.txt (/usr/bin/time -v).

hapFIRE (its pools from its own GrENE-Net SNP-only genomes, fits/<pool>.hapfire/) is scored on its own panel: SNPs with
minor-allele count at least 2 (fully imputed, so call rate 100%), truth = founder dosages x pool weights. Its pools are
drawn with the same founder list and seed as kMate's; checked.

Tables:
  af_accuracy.tsv, af_points_<scenario>.tsv.gz   allele-frequency scenarios: per class n, MAE, RMSE, R2, slope; every point
  pool_size_depth.tsv     per pool size x depth x seed: founder frequencies, allele frequencies (SNPs) and time, both tools
  founder_type.tsv        per founder: h, truth, long-read flag, per pool and arm (k-mer filter x normalization)
  units.tsv               estimation units: per founder set x mating x generations x seed x unit x class
  longread_grid.tsv       long-read founders: per pool size x depth x seed x arm (three kMate panels, hapFIRE)
  timing.tsv              speed: 50-founder pools, each tool alone on a whole node (wall-clock, CPU, peak memory)
Env: kmate. Compute node (~1 h).
"""
import os, re
from functools import lru_cache
import numpy as np, pandas as pd, scipy.sparse as sp
import config as C
from runs import figures, AF_SCENARIOS, Pool

KEY = ["chrom", "pos", "ref_len", "alt_len"]
LONG_READ = set(open(C.LONG_READ).read().split())


# ---- inputs ----------------------------------------------------------------------------------------------------------
@lru_cache(None)
def site_filter(panel, founders=None):
    """Per record of a variant panel (optionally only the given founders): passes the filter at call rate 50% / 90%,
    class, and the panel meta."""
    pre = C.VARIANT_PANELS[panel]["prefix"]
    meta = np.load(f"{pre}.meta.npz", allow_pickle=True)
    vp, vc = sp.load_npz(f"{pre}.var_pa.npz").tocsr(), sp.load_npz(f"{pre}.var_called.npz").tocsr()
    if founders:
        rows = np.isin(meta["founders"].astype(str), founders)
        vp, vc = vp[rows], vc[rows]
    F = vc.shape[0]
    n_call = np.asarray((vc > 0).sum(0)).ravel()
    n_alt = np.asarray(((vp > 0).multiply(vc > 0)).sum(0)).ravel()
    mac2 = np.minimum(n_alt, n_call - n_alt) >= 2
    L = np.maximum(meta["ref_len"], meta["alt_len"])
    cls = np.where((meta["ref_len"] == 1) & (meta["alt_len"] == 1), "SNP", np.where(L >= 50, "SV", "short"))
    return (n_call >= 0.5 * F) & mac2, (n_call >= 0.9 * F) & mac2, cls, meta


def af(fit):
    """(truth, estimate) per panel record, row-aligned; NaN where undefined."""
    T = pd.read_csv(f"{C.POOLS}/{fit.pool.id}/truth_{fit.variants}.tsv.gz", sep="\t", usecols=KEY + ["truth_af"])
    E = pd.read_csv(f"{C.FITS}/{fit.id}/fit.tsv", sep="\t", usecols=KEY + ["alt_freq"])
    meta = site_filter(fit.variants)[3]
    assert len(T) == len(E) == len(meta["pos"]) and (T[KEY].to_numpy() == E[KEY].to_numpy()).all() \
        and (T.pos.to_numpy() == meta["pos"]).all(), f"{fit.id}: truth, estimate and panel not row-aligned"
    return T.truth_af.to_numpy(), E.alt_freq.to_numpy()


def pool_weights(pool):
    w = pd.read_csv(f"{C.POOLS}/{pool.id}/pool_weights.tsv", sep="\t")
    return dict(zip(w.founder.astype(str), w.weight))


def h(fit):
    z = np.load(f"{C.FITS}/{fit.id}/fit.h_per_chrom.npz", allow_pickle=True)
    return z["founders"].astype(str), z[C.CHROM].astype(float)


def run_time(path):
    """wall-clock s, CPU s, peak memory MB from /usr/bin/time -v."""
    t = open(path).read()
    wall = [float(x) for x in re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", t).group(1).split(":")]
    cpu = sum(float(re.search(rf"{k} time \(seconds\): ([\d.]+)", t).group(1)) for k in ("User", "System"))
    rss = float(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1)) / 1024
    return sum(v * 60 ** i for i, v in enumerate(reversed(wall))), cpu, rss


# ---- metrics ---------------------------------------------------------------------------------------------------------
def af_metrics(est, truth):
    e = est - truth
    return dict(n=len(e), MAE=np.abs(e).mean(), RMSE=np.sqrt((e ** 2).mean()),
                R2=1 - (e ** 2).sum() / ((truth - truth.mean()) ** 2).sum(), slope=np.polyfit(truth, est, 1)[0])


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def h_metrics(est, truth, prefix):
    """R2, RMSE, RMSE / sd(truth) (NaN when every founder has the same truth, e.g. the 231-founder pool) and founders
    recovered (true pool members among the top-K estimates, K = pool members)."""
    ss = ((truth - truth.mean()) ** 2).sum()
    flat = ss < (1e-6 * truth.mean()) ** 2 * len(truth)
    r = rmse(est, truth)
    k = int((truth > 0).sum())
    return {f"{prefix}_h_R2": np.nan if flat else 1 - ((est - truth) ** 2).sum() / ss, f"{prefix}_h_RMSE": r,
            f"{prefix}_h_RMSE_norm": np.nan if flat else r / truth.std(),
            f"{prefix}_n_found": int(len(set(np.argsort(-est)[:k]) & set(np.flatnonzero(truth > 0))))}


# ---- hapFIRE ---------------------------------------------------------------------------------------------------------
@lru_cache(None)
def greneNet(founders=None):
    """GrENE-Net SNPs (optionally of the given founders only) with minor-allele count >= 2: sites, dosages, founder ids."""
    gt = np.load(f"{C.B}/vs_hapfire/work/greneNet_gt_matrix.npz", allow_pickle=True)
    fo = gt["founders"].astype(str); cols = np.isin(fo, founders) if founders else np.ones(len(fo), bool)
    dos = gt["dosage"][:, cols].astype(np.float64); s = dos.sum(1)
    keep = np.minimum(s, dos.shape[1] - s) >= 2
    sites = pd.DataFrame({"chrom": gt["chrom"].astype(str), "pos": gt["pos"]})[keep].reset_index(drop=True)
    assert not sites.duplicated().any(), "GrENE-Net sites repeat (chrom, pos): the hapFIRE join would pair wrong records"
    return sites, dos[keep], fo[cols]


def hapfire(fit, founders=None):
    """hapFIRE's founder frequencies with their truth, SNP allele frequencies with their truth, and time."""
    w = pool_weights(fit.pool); d = f"{C.FITS}/{fit.id}"
    e = pd.read_csv(f"{d}/fit_ecotype_frequency.txt", sep="\t", header=None, names=["founder", "freq"], dtype={"founder": str})
    sites, dos, fo = greneNet(founders)
    truth = sites.assign(truth_af=dos @ np.array([w.get(f, 0.0) for f in fo]))
    sf = pd.read_csv(f"{d}/fit_snp_frequency.txt", sep="\t", header=None, names=["chrom", "pos", "freq"],
                     dtype={"chrom": str}).merge(truth, on=["chrom", "pos"], validate="one_to_one")
    assert len(sf) > 0.9 * len(truth), f"{fit.id}: hapFIRE estimates only {len(sf)} of {len(truth)} SNPs"
    ta, tb = run_time(f"{d}/align_time.txt"), run_time(f"{d}/hapfire_time.txt")
    return dict(h=e.freq.to_numpy(), h_truth=np.array([w.get(f, 0.0) for f in e.founder]), weights=w,
                af_est=sf.freq.to_numpy(), af_truth=sf.truth_af.to_numpy(), align=ta, hapfire=tb)


def same_founders(w, v, what):
    assert {k: round(x, 9) for k, x in w.items() if x} == {k: round(x, 9) for k, x in v.items() if x}, \
        f"hapFIRE pool differs from kMate's: {what}"


# ---- tables ----------------------------------------------------------------------------------------------------------
def af_scenarios():
    rows = []
    for name, pool in AF_SCENARIOS.items():
        fit = next(f for f in figures()["af_scenarios"] if f.pool == pool)
        truth, est = af(fit)
        keep, _, cls, _ = site_filter("p231")
        m = keep & np.isfinite(truth) & np.isfinite(est)        # truth undefined where no pool founder is called
        pd.DataFrame(dict(cls=cls[m], truth_af=truth[m], est_af=est[m])).to_csv(
            f"{C.RESULTS}/af_points_{name}.tsv.gz", sep="\t", index=False, float_format="%.5g")
        rows += [dict(scenario=name, cls=c, **af_metrics(est[m & (cls == c)], truth[m & (cls == c)])) for c in ("SNP", "SV", "short")]
    return pd.DataFrame(rows)


def pool_size_depth():
    keep, keep90, cls, _ = site_filter("p231")
    snp = keep & (cls == "SNP"); snp90 = keep90 & (cls == "SNP")
    fits = figures()["pool_size_depth"]
    hf_fit = {(f.pool.n, f.pool.depth, f.pool.seed): f for f in fits if f.tool == "hapfire"}
    rows = []
    for fit in (f for f in fits if f.tool == "kmate"):
        p = fit.pool; w = pool_weights(p)
        fo, hk = h(fit); t = np.array([w.get(f, 0.0) for f in fo])
        wall, cpu, rss = run_time(f"{C.FITS}/{fit.id}/time.txt")
        truth, est = af(fit); ok = np.isfinite(truth) & np.isfinite(est)
        r = dict(N=p.n, depth=p.depth, seed=p.seed, n_truth_founders=int((t > 0).sum()), **h_metrics(hk, t, "kmate"),
                 kmate_elapsed_s=wall, kmate_cpu_s=cpu, kmate_max_rss_mb=rss)
        for suf, m in (("", ok & snp), ("_cr90", ok & snp90)):
            a = af_metrics(est[m], truth[m])
            r.update({f"kmate_af_R2{suf}": a["R2"], f"kmate_af_RMSE{suf}": a["RMSE"], f"kmate_n_snps{suf}": a["n"]})
        hf = hapfire(hf_fit[p.n, p.depth, p.seed]); same_founders(w, hf["weights"], fit.id)
        a = af_metrics(hf["af_est"], hf["af_truth"])
        r.update(h_metrics(hf["h"], hf["h_truth"], "hapfire"), hapfire_af_R2=a["R2"], hapfire_af_RMSE=a["RMSE"], hapfire_n_snps=a["n"],
                 hapfire_af_RMSE_cr90=a["RMSE"],                                      # fully imputed: call rate 100%
                 hapfire_align_elapsed_s=hf["align"][0], hapfire_align_cpu_s=hf["align"][1],
                 hapfire_proper_elapsed_s=hf["hapfire"][0], hapfire_proper_cpu_s=hf["hapfire"][1],
                 hapfire_total_elapsed_s=hf["align"][0] + hf["hapfire"][0], hapfire_total_cpu_s=hf["align"][1] + hf["hapfire"][1],
                 hapfire_max_rss_mb=max(hf["align"][2], hf["hapfire"][2]))
        rows.append(r)
    return pd.DataFrame(rows)


def founder_type():
    rows = []
    for fit in figures()["founder_type"]:
        p = fit.pool; w = pool_weights(p)
        fo, hk = h(fit)
        arm = fit.kmers.removeprefix("p231_") + ("_global" if fit.normalize == "global" else "")   # unfilt, filt, *_global
        regime = f"n{p.n}_{p.founder_type}heavy" if p.founder_type else f"n{p.n}_g{p.gen}"
        rows.append(pd.DataFrame(dict(regime=regime, arm=arm, founder=fo, h=hk, truth=[w.get(f, 0.0) for f in fo],
                                      long_read=[f in LONG_READ for f in fo])))
    return pd.concat(rows)


def units():
    rows = []
    for fit in figures()["units_p231"] + figures()["units_p78"]:
        p = fit.pool
        keep, _, cls, _ = site_filter(fit.variants)
        truth, est = af(fit)
        for c in ("SNP", "short", "SV"):
            s = keep & (cls == c) & np.isfinite(truth)
            e, t = est[s], truth[s]; has = np.isfinite(e)                   # windows under 200 k-mers abstain (NaN)
            a = af_metrics(e[has], t[has])
            rows.append(dict(founders=p.founders, mating="selfing 97%" if p.selfing else "outcrossing", generations=p.gen,
                             seed=p.seed, unit=fit.unit, window_kb=int(fit.unit[2:-1]) if fit.unit.startswith("bp") else np.nan,
                             cls=c, n_pass=int(s.sum()), frac_estimated=has.mean(), MAE=a["MAE"], RMSE=a["RMSE"], R2=a["R2"]))
    return pd.DataFrame(rows)


def longread_grid():
    f78 = tuple(sorted(LONG_READ))
    rows = []
    for fit in figures()["longread_grid"]:
        p = fit.pool; w = pool_weights(p); assert set(f for f, v in w.items() if v) <= set(f78)
        if fit.tool == "hapfire":
            hf = hapfire(fit, f78); same_founders(w, pool_weights(Pool("p78", p.n, seed=p.seed, depth=p.depth)), fit.id)
            r = rmse(hf["af_est"], hf["af_truth"])
            rows.append(dict(N=p.n, depth=p.depth, seed=p.seed, arm="hapfire", h_rmse78=rmse(hf["h"], hf["h_truth"]),
                             h_misassigned=0.5 * np.abs(hf["h"] - hf["h_truth"]).sum(), h_on_short=0.0,
                             h_true_founder=hf["h"][hf["h_truth"] > 0][0] if p.n == 1 else np.nan,
                             af_rmse=r, af_n=len(hf["af_est"]), af_rmse_cr90=r, af_n_cr90=len(hf["af_est"]),
                             elapsed_s=hf["align"][0] + hf["hapfire"][0], cpu_s=hf["align"][1] + hf["hapfire"][1]))
            continue
        fo, hk = h(fit); t = np.array([w.get(f, 0.0) for f in fo]); in78 = np.isin(fo, f78)
        keep, keep90, cls, _ = site_filter(fit.variants, f78)                 # the filter among the 78 founders
        truth, est = af(fit); ok = np.isfinite(truth) & np.isfinite(est); snp = cls == "SNP"
        wall, cpu, _ = run_time(f"{C.FITS}/{fit.id}/time.txt")
        rows.append(dict(N=p.n, depth=p.depth, seed=p.seed, arm=fit.kmers, h_rmse78=rmse(hk[in78], t[in78]),
                         h_misassigned=0.5 * np.abs(hk - t).sum(), h_on_short=hk[~in78].sum(),
                         h_true_founder=hk[t > 0][0] if p.n == 1 else np.nan,
                         af_rmse=rmse(est[ok & keep & snp], truth[ok & keep & snp]), af_n=int((ok & keep & snp).sum()),
                         af_rmse_cr90=rmse(est[ok & keep90 & snp], truth[ok & keep90 & snp]), af_n_cr90=int((ok & keep90 & snp).sum()),
                         af_rmse_all=rmse(est[ok & keep], truth[ok & keep]), af_n_all=int((ok & keep).sum()),
                         elapsed_s=wall, cpu_s=cpu))
    return pd.DataFrame(rows)


def timing():
    rows = {}
    for fit in figures()["timing"]:
        p = fit.pool; r = rows.setdefault((p.n, p.depth, p.seed), dict(N=p.n, depth=p.depth, seed=p.seed))
        d = f"{C.FITS}/{fit.id}"
        if fit.tool == "kmate":
            r["kmate_elapsed_s"], r["kmate_cpu_s"], r["kmate_max_rss_mb"] = run_time(f"{d}/time.txt")
        else:
            a, b = run_time(f"{d}/align_time.txt"), run_time(f"{d}/hapfire_time.txt")
            r.update(hapfire_align_elapsed_s=a[0], hapfire_proper_elapsed_s=b[0], hapfire_total_elapsed_s=a[0] + b[0],
                     hapfire_total_cpu_s=a[1] + b[1], hapfire_max_rss_mb=max(a[2], b[2]))
    return pd.DataFrame(rows.values())


if __name__ == "__main__":
    os.makedirs(C.RESULTS, exist_ok=True)
    for name, make in (("af_accuracy", af_scenarios), ("pool_size_depth", pool_size_depth), ("founder_type", founder_type),
                       ("units", units), ("longread_grid", longread_grid), ("timing", timing)):
        T = make()
        T.to_csv(f"{C.RESULTS}/{name}.tsv", sep="\t", index=False)
        print(f"wrote {name}.tsv ({len(T):,} rows)", flush=True)
