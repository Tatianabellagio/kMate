#!/usr/bin/env python
"""Does SV content track CLIMATIC distance to the Col-0 backbone, or just the climate gradient?

The panel is TAIR10/Col-0-referenced, so every SV measure is defined against Col-0. §2.3 of
the notebook showed SV content is not explained by GENETIC divergence from Col-0 (SNP ALT
count is not climate-structured, rho=-0.063 p=0.34). This asks the complementary question the
user raised: is SV content related to how climatically UNLIKE Col-0's home a founder's home is?

Two hypotheses that must be told apart, because they make different predictions:

  DIRECTIONAL   SV content varies monotonically with home climate (colder home -> more
                inserted sequence). This is what §2.2/§2.4 found (bp_ins_frac vs home bio1
                rho = -0.500).
  DISTANCE      SV content rises with |home climate - Col-0 home climate| in either
                direction. This is the shape a reference-related artifact would take: the
                further a genome's environment from the reference's, the more of it looks
                "inserted".

THE TRAP. Col-0's recorded bio1 (13.1 C) sits in the WARM half of the panel (founder median
~9.6 C), not at an extreme. So |delta bio1| is largely "how cold is this founder", i.e. a
partly monotone re-encoding of the directional gradient. A raw correlation of SV content with
distance therefore proves nothing on its own. Both are fitted jointly, and each is partialled
out of the other.

CAVEAT ON COL-0'S CLIMATE RECORD -- read before using any number here. The 1001G table places
ecotype 6909 in the USA at 38.3 N, -92.3 W (Columbia, Missouri) with bio1 13.1 C, bio12 1023
mm. That is where the Laibach/Redei lineage was propagated and named, NOT where Col-0 was
collected; the accession's origin is conventionally placed near Landsberg an der Warthe /
Gorzow Wielkopolski (Poland-Germany, ~52.7 N 15.2 E), a markedly colder and drier-summer
climate. Distances computed from the recorded value are therefore distances to central
Missouri. A European-origin sensitivity arm is included, anchored on the panel founders
nearest that putative origin.

Env: kmate. Reads founder_sv_content.npz, founder_climate_confound.npz, the 1001g bioclim
and geo tables. Writes results/sv_adaptive/founder_col0_climate_distance.csv.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
KEY = "/global/scratch/users/tbellg/gea_grene-net/key_files"
BIOS = [f"bio{i}" for i in range(1, 20)]
COL0 = "6909"
ORIGIN_LAT, ORIGIN_LON, ORIGIN_R = 52.75, 15.25, 400.0   # putative European origin, km radius

MEAS = [("n_ins", "ins, count"), ("bp_ins", "ins, bp"),
        ("n_ins_frac", "ins, count / load"), ("bp_ins_frac", "ins, bp / load"),
        ("n_del", "del, count"), ("bp_del", "del, bp"),
        ("n_del_frac", "del, count / load"), ("bp_del_frac", "del, bp / load"),
        ("n_sv", "SV, count"), ("bp_sv", "SV, bp"),
        ("n_sv_frac", "SV, count / load"), ("bp_sv_frac", "SV, bp / load")]


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None); by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    S = np.load(f"{OUT}/founder_sv_content.npz", allow_pickle=True)
    keep = C["keep"].astype(bool)
    fid = C["founders"].astype("U6")[keep]
    g1 = C["gamma_bio1"][keep]
    M = {k: S[k][keep].astype(float) for k, _ in MEAS}

    bio = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
    bio["ecotypeid"] = bio["ecotypeid"].astype(str)
    bm = bio.set_index("ecotypeid")
    H = np.column_stack([[float(bm.loc[f, b]) for f in fid] for b in BIOS])
    home1 = H[:, 0]; home18 = H[:, 17]

    geo = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv")
    geo["ecotype_id"] = geo["ecotype_id"].astype(str)
    gm = geo.set_index("ecotype_id")
    lat = np.array([float(gm.loc[f, "Latitude_corrected"]) for f in fid])
    lon = np.array([float(gm.loc[f, "Longitude_corrected"]) for f in fid])

    # standardize bioclim on the FOUNDER panel, then place Col-0 in the same space
    mu, sd = H.mean(0), H.std(0)
    Hz = (H - mu) / sd
    c0 = (np.array([float(bm.loc[COL0, b]) for b in BIOS]) - mu) / sd
    c0_bio1 = float(bm.loc[COL0, "bio1"]); c0_bio18 = float(bm.loc[COL0, "bio18"])
    g0 = gm.loc[COL0]
    print(f"[Col-0 record] {g0['country']} lat={g0['Latitude_corrected']} lon={g0['Longitude_corrected']}"
          f"  bio1={c0_bio1:.2f}C bio18={c0_bio18:.0f}mm")
    print(f"[panel]  founder bio1: min {home1.min():.1f} median {np.median(home1):.1f} max {home1.max():.1f}")
    print(f"  -> Col-0 sits at percentile {100*(home1 < c0_bio1).mean():.0f} of the founder bio1 range "
          f"(NOT an extreme; |delta bio1| is therefore partly a re-encoding of 'how cold')\n")

    # --- European-origin sensitivity anchor -----------------------------------
    R = 6371.0
    d_km = R * np.arccos(np.clip(
        np.sin(np.deg2rad(ORIGIN_LAT)) * np.sin(np.deg2rad(lat)) +
        np.cos(np.deg2rad(ORIGIN_LAT)) * np.cos(np.deg2rad(lat)) *
        np.cos(np.deg2rad(lon - ORIGIN_LON)), -1, 1))
    near = d_km <= ORIGIN_R
    print(f"[EU-origin anchor] {int(near.sum())} founders within {ORIGIN_R:.0f} km of "
          f"({ORIGIN_LAT}, {ORIGIN_LON}); their median bio1 = {np.median(home1[near]):.2f}C")
    eu = np.median(Hz[near], axis=0)
    eu_bio1 = float(np.median(home1[near]))

    # --- distance variables ---------------------------------------------------
    D = {
        "signed d bio1 (home - Col0)":      home1 - c0_bio1,
        "|d bio1| to Col0":                 np.abs(home1 - c0_bio1),
        "|d bio18| to Col0":                np.abs(home18 - c0_bio18),
        "19-axis climate dist to Col0":     np.linalg.norm(Hz - c0, axis=1),
        "|d bio1| to EU-origin anchor":     np.abs(home1 - eu_bio1),
        "19-axis dist to EU-origin anchor": np.linalg.norm(Hz - eu, axis=1),
    }
    print("How much is each distance just the directional gradient?")
    for k, v in D.items():
        r = stats.spearmanr(v, home1).correlation
        print(f"  {k:<36} vs home bio1: rho={r:+.3f}")

    rows = []
    print("\n" + "=" * 116)
    print("SV measures vs climatic distance to Col-0   (and each partialled against the directional gradient)")
    for dname, dv in D.items():
        print(f"\n--- {dname} ---")
        print(f"{'measure':<22}{'raw rho':>10}{'p':>11}{'| home bio1':>14}{'p':>11}")
        for k, lab in MEAS:
            r, p = stats.spearmanr(M[k], dv)
            rp, pp = partial_spearman(M[k], dv, home1[:, None])
            star = "***" if pp < 1e-3 else "**" if pp < 0.01 else "*" if pp < 0.05 else " "
            print(f"{lab:<22}{r:>+10.3f}{p:>11.1e}{rp:>+14.3f}{pp:>11.1e} {star}")
            rows.append(dict(distance=dname, measure=k, label=lab, rho=r, p=p,
                             partial_home_bio1=rp, p_partial=pp))

    # --- head-to-head: directional vs distance, on the headline measure -------
    print("\n" + "=" * 116)
    print("HEAD-TO-HEAD on bp_ins_frac: is it direction or distance?")
    v = M["bp_ins_frac"]
    r_dir, p_dir = stats.spearmanr(v, home1)
    r_dist, p_dist = stats.spearmanr(v, D["19-axis climate dist to Col0"])
    rp_dir, pp_dir = partial_spearman(v, home1, D["19-axis climate dist to Col0"][:, None])
    rp_dist, pp_dist = partial_spearman(v, D["19-axis climate dist to Col0"], home1[:, None])
    print(f"  directional (home bio1)        raw {r_dir:+.3f} (p={p_dir:.1e})  |  "
          f"controlling distance {rp_dir:+.3f} (p={pp_dir:.1e})")
    print(f"  distance to Col-0 (19 axes)    raw {r_dist:+.3f} (p={p_dist:.1e})  |  "
          f"controlling direction {rp_dist:+.3f} (p={pp_dist:.1e})")
    print("  -> whichever SURVIVES the other is the real variable.")

    print("\nSame head-to-head for gamma_bio1 (the founder climate response):")
    rg_dir, pg_dir = partial_spearman(g1, home1, D["19-axis climate dist to Col0"][:, None])
    rg_dist, pg_dist = partial_spearman(g1, D["19-axis climate dist to Col0"], home1[:, None])
    print(f"  gamma vs home bio1 | distance   {rg_dir:+.3f} (p={pg_dir:.1e})")
    print(f"  gamma vs distance  | home bio1  {rg_dist:+.3f} (p={pg_dist:.1e})")

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_col0_climate_distance.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_col0_climate_distance.csv")


if __name__ == "__main__":
    main()
