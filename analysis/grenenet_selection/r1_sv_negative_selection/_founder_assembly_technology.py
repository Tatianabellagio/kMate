#!/usr/bin/env python
"""Does long-read SEQUENCING TECHNOLOGY or assembly contiguity explain founder insertion content?

This was the one technical check flagged as untestable in earlier passes -- the metadata is
not in the kMate repo. It IS available in the assembly release:

    /global/scratch/users/tbellg/pang/long_read_seq_ara/ASSEMBLIES_Best_version_of_dataset.csv

carrying, per Assembly_ID: Primary_Sequencing_Technology (CLR / HiFi / ONT / ONT_R10.4 /
ONT_HiFi / RS II / Sanger), Assembler, N50_contigs, Largest_contigs, Gaps_Scaffolds,
Sum_Scaffolds, and publication fields.

Why it matters. HiFi, CLR and ONT differ substantially in how well they resolve insertions
relative to a reference -- HiFi is accurate in repeats, CLR is error-prone, ONT is long but
noisier. If the accessions sequenced with the better technology happened to come from colder
places, "cold-origin ecotypes carry more inserted sequence" would be a platform artifact.

Tests (cactus founders only -- PanGenie founders have no assembly):
  T1  technology vs ecotype origin climate      (Kruskal-Wallis across platforms)
  T2  technology vs inserted kb                 (Kruskal-Wallis)
  T3  contiguity (N50, gaps) vs inserted kb and vs origin climate
  T4  does kb ~ origin climate hold WITHIN each technology group?
  T5  partial correlation of kb ~ origin / kb ~ gamma with technology (dummy-coded),
      N50 and assembly size all held fixed

Env: kmate. Writes results/sv_adaptive/founder_assembly_technology.csv + printed report.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
ASM_DIR = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only"
TECH_CSV = "/global/scratch/users/tbellg/pang/long_read_seq_ara/ASSEMBLIES_Best_version_of_dataset.csv"


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) if len(np.unique(Z[:, j])) > 2 else Z[:, j]
                          for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None); by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    GRp = np.load(f"{OUT}/founder_graph_representation.npz", allow_pickle=True)
    keep = C["keep"].astype(bool)
    fid = C["founders"].astype("U6")[keep]
    home1 = C["home_bio1"][keep]; g1 = C["gamma_bio1"][keep]
    KB = GRp["kb_ins"]; asm_mb = GRp["asm_mb"]

    have_asm = {os.path.basename(f).split(".")[0] for f in glob.glob(f"{ASM_DIR}/*.chr.fa")}
    T = pd.read_csv(TECH_CSV)
    T["Assembly_ID"] = T["Assembly_ID"].astype(str)
    T["Accession_ID"] = T["Accession_ID"].astype(str)
    T = T[T.Assembly_ID.isin(have_asm)].copy()
    # normalise a casing typo in the source table
    T["tech"] = T["Primary_Sequencing_Technology"].astype(str).str.replace("ONT_HiFI", "ONT_HiFi")
    print(f"[meta] {len(T)} of {len(have_asm)} graph assemblies matched in the metadata table")

    acc2 = {}
    for _, r0 in T.iterrows():
        acc2.setdefault(r0["Accession_ID"], r0)
    tech = np.array([acc2[f]["tech"] if f in acc2 else "" for f in fid], dtype=object)
    n50 = np.array([float(acc2[f]["N50_contigs"]) if f in acc2 else np.nan for f in fid])
    gaps = np.array([float(acc2[f]["Gaps_Scaffolds"]) if f in acc2 else np.nan for f in fid])
    asmb = np.array([str(acc2[f]["Assembler"]) if f in acc2 else "" for f in fid], dtype=object)
    ok = np.array([t != "" for t in tech]) & np.isfinite(KB) & np.isfinite(home1)
    print(f"[match] {int(ok.sum())} founders with technology metadata\n")

    rows = []
    print("=" * 88)
    print("T1/T2  Technology groups")
    print(f"{'technology':<12}{'n':>4}{'origin bio1':>14}{'kb inserted':>14}{'N50 (Mb)':>12}{'asm (Mb)':>11}")
    groups = []
    for t in sorted(set(tech[ok])):
        m = ok & (tech == t)
        if m.sum() < 2:
            continue
        groups.append((t, m))
        print(f"{t:<12}{int(m.sum()):>4}{np.median(home1[m]):>14.2f}{np.median(KB[m]):>14.1f}"
              f"{np.nanmedian(n50[m])/1e6:>12.2f}{np.nanmedian(asm_mb[m]):>11.1f}")
    if len(groups) >= 2:
        kw_h = stats.kruskal(*[home1[m] for _, m in groups])
        kw_k = stats.kruskal(*[KB[m] for _, m in groups])
        kw_g = stats.kruskal(*[g1[m] for _, m in groups])
        print(f"\n  Kruskal-Wallis across platforms:")
        print(f"    ecotype origin bio1 : H={kw_h.statistic:.2f} p={kw_h.pvalue:.4f}")
        print(f"    kb inserted         : H={kw_k.statistic:.2f} p={kw_k.pvalue:.4f}")
        print(f"    gamma_bio1          : H={kw_g.statistic:.2f} p={kw_g.pvalue:.4f}")
        rows += [dict(test="KW", var=v, H=s.statistic, p=s.pvalue) for v, s in
                 (("origin_bio1", kw_h), ("kb_ins", kw_k), ("gamma", kw_g))]
    print(f"\n  assemblers: {pd.Series(asmb[ok]).value_counts().head(6).to_dict()}")

    print("\n" + "=" * 88)
    print("T3  Contiguity")
    for nm, v in (("N50_contigs", n50), ("Gaps_Scaffolds", gaps), ("assembly size Mb", asm_mb)):
        m = ok & np.isfinite(v)
        r_k, p_k = stats.spearmanr(v[m], KB[m])
        r_h, p_h = stats.spearmanr(v[m], home1[m])
        r_g, p_g = stats.spearmanr(v[m], g1[m])
        print(f"  {nm:<18} vs kb inserted {r_k:+.3f} (p={p_k:.3f}) | vs origin bio1 {r_h:+.3f} "
              f"(p={p_h:.3f}) | vs gamma {r_g:+.3f} (p={p_g:.3f})   [n={int(m.sum())}]")
        rows.append(dict(test="T3", var=nm, rho_kb=r_k, p_kb=p_k, rho_home=r_h, p_home=p_h))

    print("\n" + "=" * 88)
    print("T4  Does kb ~ origin climate hold WITHIN each technology group?")
    for t, m in groups:
        if m.sum() < 6:
            print(f"  {t:<12} n={int(m.sum())} -- too few")
            continue
        r_h, p_h = stats.spearmanr(KB[m], home1[m])
        r_g, p_g = stats.spearmanr(KB[m], g1[m])
        print(f"  {t:<12} n={int(m.sum()):>3}  kb~origin rho={r_h:+.3f} (p={p_h:.4f})   "
              f"kb~gamma rho={r_g:+.3f} (p={p_g:.4f})")
        rows.append(dict(test="T4", var=t, n=int(m.sum()), rho_home=r_h, p_home=p_h,
                         rho_gamma=r_g, p_gamma=p_g))

    print("\n" + "=" * 88)
    print("T5  Everything held fixed at once (technology dummies + N50 + assembly size)")
    tl = sorted(set(tech[ok]))
    D = np.column_stack([(tech[ok] == t).astype(float) for t in tl[:-1]])   # drop one level
    Z = np.column_stack([D, n50[ok], asm_mb[ok]])
    Z = Z[:, np.isfinite(Z).all(0)] if np.isfinite(Z).all() else np.column_stack(
        [D, np.nan_to_num(n50[ok], nan=np.nanmedian(n50[ok])),
         np.nan_to_num(asm_mb[ok], nan=np.nanmedian(asm_mb[ok]))])
    r0h, p0h = stats.spearmanr(KB[ok], home1[ok])
    r0g, p0g = stats.spearmanr(KB[ok], g1[ok])
    rh, ph = partial_spearman(KB[ok], home1[ok], Z)
    rg, pg = partial_spearman(KB[ok], g1[ok], Z)
    print(f"  kb ~ origin bio1 : raw {r0h:+.3f} (p={p0h:.1e})  ->  controlled {rh:+.3f} (p={ph:.1e})")
    print(f"  kb ~ gamma_bio1  : raw {r0g:+.3f} (p={p0g:.1e})  ->  controlled {rg:+.3f} (p={pg:.1e})")
    rows.append(dict(test="T5", var="kb_origin", rho_raw=r0h, rho_ctrl=rh, p_ctrl=ph))
    rows.append(dict(test="T5", var="kb_gamma", rho_raw=r0g, rho_ctrl=rg, p_ctrl=pg))

    np.savez_compressed(f"{OUT}/founder_assembly_technology.npz",
                        founders=fid, tech=tech.astype("U12"), n50=n50, gaps=gaps,
                        asm_mb=asm_mb, kb_ins=KB, ok=ok)
    pd.DataFrame(rows).to_csv(f"{OUT}/founder_assembly_technology.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_assembly_technology.{{npz,csv}}")


if __name__ == "__main__":
    main()
