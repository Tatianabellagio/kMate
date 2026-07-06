#!/usr/bin/env python
"""SuSiE-RSS per block + driver-vs-passenger readout.

For each fine-map-input block npz (site-level z + founder LD R), run susie_rss
(n=31 sites) and read out, per block, whether an SV is the causal/lead variant
beyond the best SNP:

  * top-PIP variant class (is the SV the single strongest signal?)
  * max SV PIP vs max SNP PIP
  * credible sets: does any 95% CS contain an SV? (SV in a CS = a driver/independent
    signal the SNPs can't fully explain; SV absent from every CS = passenger)

Aggregates across the SV-containing selected blocks -> the headline count:
how many blocks have an SV as a driver vs passenger.

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY run_susie.py
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA_DIR = os.path.dirname(HERE)
sys.path.insert(0, GEA_DIR)
import lib

RSCRIPT = "/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript"
SUSIE_R = f"{HERE}/susie_one.R"
DP = f"{lib.GEA}/driver_passenger"


def run_block(npz_path: str, n: int, tmpd: str, env: dict):
    d = dict(np.load(npz_path, allow_pickle=True))
    block = str(d["block"]); z = d["z"].astype(float); R = d["R"].astype(float)
    m = len(z)
    zf = f"{tmpd}/{block}.z.txt"; Rf = f"{tmpd}/{block}.R.txt"; of = f"{tmpd}/{block}.pip.csv"
    np.savetxt(zf, z)
    np.savetxt(Rf, R)
    subprocess.run([RSCRIPT, SUSIE_R, zf, Rf, str(n), of], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    pipdf = pd.read_csv(of)
    pip = pipdf["pip"].to_numpy(); cs = pipdf["cs"].to_numpy()
    is_sv = d["is_sv"].astype(bool); cls = d["cls"].astype(str)
    pos = d["pos"]; zz = d["z"]
    top = int(np.nanargmax(pip)) if np.isfinite(pip).any() else 0
    sv_pip = pip[is_sv]; snp_pip = pip[cls == "snp"]
    n_cs = int(cs.max())
    # credible sets containing an SV
    cs_with_sv = sorted({int(c) for c in cs[is_sv] if c > 0})
    # per-CS composition (size + class of its top-PIP member)
    cs_info = []
    for c in range(1, n_cs + 1):
        mask = cs == c
        if not mask.any():
            continue
        j = np.where(mask)[0]
        lead = j[np.argmax(pip[j])]
        cs_info.append(dict(cs=c, size=int(mask.sum()), has_sv=bool(is_sv[mask].any()),
                            lead_cls=cls[lead], lead_pip=float(pip[lead])))
    return dict(
        block=block, m=m, n_sv=int(is_sv.sum()),
        top_cls=cls[top], top_pip=float(pip[top]), top_is_sv=bool(is_sv[top]),
        max_sv_pip=float(np.nanmax(sv_pip)) if is_sv.any() else np.nan,
        max_snp_pip=float(np.nanmax(snp_pip)) if (cls == "snp").any() else np.nan,
        n_cs=n_cs, n_cs_with_sv=len(cs_with_sv),
        sv_max_absz=float(np.abs(zz[is_sv]).max()) if is_sv.any() else np.nan,
        max_absz=float(np.abs(zz).max()),
        cs_detail=";".join(f"CS{ci['cs']}[n={ci['size']},lead={ci['lead_cls']}"
                           f"{'+SV' if ci['has_sv'] else ''},pip={ci['lead_pip']:.2f}]"
                           for ci in cs_info),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=f"{DP}/finemap_inputs")
    ap.add_argument("--n", type=int, default=31, help="effective sample size (sites)")
    ap.add_argument("--out", default=f"{DP}/finemap_susie_results.csv")
    args = ap.parse_args()

    env = dict(os.environ); env["LD_LIBRARY_PATH"] = "/usr/lib64:" + env.get("LD_LIBRARY_PATH", "")
    npzs = sorted(glob.glob(f"{args.indir}/*.npz"))
    print(f"[SuSiE-RSS] {len(npzs)} blocks, n={args.n}", flush=True)
    rows = []
    with tempfile.TemporaryDirectory(dir=args.indir) as tmpd:
        for f in npzs:
            try:
                r = run_block(f, args.n, tmpd, env)
            except subprocess.CalledProcessError as e:
                print(f"  ! {os.path.basename(f)} R FAILED: {e.stderr.decode()[:200]}", flush=True)
                continue
            rows.append(r)
            drv = "SV-DRIVER" if (r["top_is_sv"] or r["n_cs_with_sv"] > 0) else "passenger"
            print(f"  {r['block']:12s} m={r['m']:>4} sv={r['n_sv']:>2} | top={r['top_cls']:>10}"
                  f"(pip{r['top_pip']:.2f}) SVpip={r['max_sv_pip']:.2f} SNPpip={r['max_snp_pip']:.2f}"
                  f" | CS={r['n_cs']}(SV in {r['n_cs_with_sv']}) -> {drv}", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(args.out, index=False)

    print("\n" + "=" * 70)
    print("DRIVER-VS-PASSENGER SUMMARY (SV-containing selected blocks)")
    print("=" * 70)
    nb = len(res)
    driver = (res["top_is_sv"] | (res["n_cs_with_sv"] > 0))
    print(f"  blocks fine-mapped              : {nb}")
    print(f"  any credible set at all         : {int((res['n_cs']>0).sum())}/{nb}")
    print(f"  SV is top-PIP variant           : {int(res['top_is_sv'].sum())}/{nb}")
    print(f"  SV in >=1 95% credible set      : {int((res['n_cs_with_sv']>0).sum())}/{nb}")
    print(f"  => SV-DRIVER (either)           : {int(driver.sum())}/{nb}")
    print(f"  => SV passenger                 : {int((~driver).sum())}/{nb}")
    print(f"  SV max-PIP >= best SNP max-PIP  : {int((res['max_sv_pip']>=res['max_snp_pip']).sum())}/{nb}")
    print(f"  SV max-PIP: median={res['max_sv_pip'].median():.3f} max={res['max_sv_pip'].max():.3f}")
    print("=" * 70)
    print(f"  -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
