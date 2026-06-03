#!/usr/bin/env python
"""Panel + kMate-input statistics for the paper.

Production panel (arch3, 231 founders):
  VCF   : panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz
  V_pa  : panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz
  K_pa  : data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr{N}.{kmer_pa,meta}.npz

Computes, per chrom and genome-wide:
  - VCF size + record count
  - matrix file sizes + dimensions
  - record class: SNP / small indel (<=50bp) / SV (>50bp); INS/DEL/MNP; SV size bins
  - allele-count / frequency spectrum (carriers among called founders)
  - missingness (per-record call rate; overall genotype missing rate)
  - K_pa: n k-mers, n bubbles, k-mers/bubble, founders/k-mer (density), AC range
Writes a markdown report to results/panel_stats/PANEL_STATS.md
"""
import os, gc, glob
import numpy as np

PROJ = "/global/scratch/users/tbellg/kmate"
NF = 231                      # founders
SV_BP = 50                    # SV threshold (length, bp)
CHROMS = [1, 2, 3, 4, 5]
OUTDIR = f"{PROJ}/results/panel_stats"
os.makedirs(OUTDIR, exist_ok=True)


def mb(path):
    return os.path.getsize(path) / 1e6 if os.path.exists(path) else float("nan")


def colcounts(npz_path, ncols):
    """Per-column nnz of a CSR (rows=founders): bincount of the column-index array."""
    z = np.load(npz_path, allow_pickle=False)
    idx = z["indices"]
    c = np.bincount(idx, minlength=ncols).astype(np.int64)
    del idx, z
    gc.collect()
    return c


agg = {}        # genome-wide accumulators
rows_md = []

for c in CHROMS:
    cl = f"chr{c}"; Cl = f"Chr{c}"
    vcf = f"{PROJ}/panel/arch3/{cl}/merged_231_{cl}_final.vcf.gz"
    vp  = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.var_pa.npz"
    vc  = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.var_called.npz"
    vm  = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz"
    kp  = f"{PROJ}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{Cl}.kmer_pa.npz"
    km  = f"{PROJ}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{Cl}.meta.npz"

    m = np.load(vm, allow_pickle=True)
    pos = m["pos"]; rl = m["ref_len"].astype(np.int64); al = m["alt_len"].astype(np.int64)
    N = len(rl)
    Lmax = np.maximum(rl, al)
    Ldiff = np.abs(al - rl)

    snp = (rl == 1) & (al == 1)
    sv  = Ldiff > SV_BP                          # SV by length change >50bp
    small = (~snp) & (~sv)                        # everything else: small indels/MNPs <=50bp
    ins = al > rl; dele = rl > al; mnp = (rl == al) & (rl > 1)

    # carriers (AC) and called per record
    car = colcounts(vp, N)                        # founders carrying ALT
    cal = colcounts(vc, N)                         # founders with a called genotype
    af = np.where(cal > 0, car / cal, np.nan)      # ALT freq among called
    miss_frac = 1.0 - cal / NF                      # per-record missing fraction

    # K_pa
    kmeta = np.load(km, allow_pickle=True)
    nk = len(kmeta["kmer_index"]); nb = len(kmeta["bubble_start"])
    kpz = np.load(kp, allow_pickle=False)
    k_nnz = len(kpz["indices"])                    # total founder-kmer presences
    kfc = colcounts(kp, nk)                        # founders per k-mer
    del kpz; gc.collect()

    d = dict(
        N=N, vcf_mb=mb(vcf),
        vp_mb=mb(vp)+mb(vc)+mb(vm), kp_mb=mb(kp)+mb(km),
        n_snp=int(snp.sum()), n_small=int(small.sum()), n_sv=int(sv.sum()),
        n_ins=int(ins.sum()), n_del=int(dele.sum()), n_mnp=int(mnp.sum()),
        sv_50_100=int(((Ldiff>50)&(Ldiff<=100)).sum()),
        sv_100_1k=int(((Ldiff>100)&(Ldiff<=1000)).sum()),
        sv_1k_10k=int(((Ldiff>1000)&(Ldiff<=10000)).sum()),
        sv_10k=int((Ldiff>10000).sum()),
        sv_max=int(Ldiff.max()),
        ac1=int((car==1).sum()), ac_le5=int((car<=5).sum()),
        called_entries=int(cal.sum()),
        miss_mean=float(miss_frac.mean()), miss_med=float(np.median(miss_frac)),
        n_pos=int(np.unique(pos).size),
        nk=nk, nb=nb, k_nnz=k_nnz,
        kfc_min=int(kfc.min()), kfc_max=int(kfc.max()), kfc_mean=float(kfc.mean()),
    )
    for k, v in d.items():
        agg[k] = agg.get(k, 0) + v if isinstance(v, int) else agg.get(k, [])
    # store per-chrom for weighted means later
    agg.setdefault("_per", []).append((c, d, car, cal, af, miss_frac, kfc))

    print(f"Chr{c}: N={N:,} SNP={d['n_snp']:,} small={d['n_small']:,} SV={d['n_sv']:,} "
          f"| miss={d['miss_mean']*100:.2f}% | kmers={nk:,} bubbles={nb:,}")
    del m, pos, rl, al, Lmax, Ldiff, car, cal, af, miss_frac, kfc, kmeta
    gc.collect()

# ---- genome-wide rollups ----
per = agg["_per"]
N_tot = sum(d["N"] for _, d, *_ in per)
car_all = np.concatenate([p[2] for p in per])
cal_all = np.concatenate([p[3] for p in per])
af_all  = np.concatenate([p[4] for p in per])
miss_all = np.concatenate([p[5] for p in per])
kfc_all = np.concatenate([p[6] for p in per])

def s(key): return sum(d[key] for _, d, *_ in per)

genome = dict(
    N=N_tot,
    vcf_mb=s("vcf_mb"), vp_mb=s("vp_mb"), kp_mb=s("kp_mb"),
    n_snp=s("n_snp"), n_small=s("n_small"), n_sv=s("n_sv"),
    n_ins=s("n_ins"), n_del=s("n_del"), n_mnp=s("n_mnp"),
    sv_50_100=s("sv_50_100"), sv_100_1k=s("sv_100_1k"),
    sv_1k_10k=s("sv_1k_10k"), sv_10k=s("sv_10k"),
    sv_max=max(d["sv_max"] for _, d, *_ in per),
    ac1=s("ac1"), ac_le5=s("ac_le5"),
    overall_miss=1.0 - cal_all.sum() / (NF * N_tot),
    miss_mean=float(miss_all.mean()), miss_med=float(np.median(miss_all)),
    nk=s("nk"), nb=s("nb"), k_nnz=s("k_nnz"),
    kfc_mean=float(kfc_all.mean()), kfc_min=int(kfc_all.min()), kfc_max=int(kfc_all.max()),
    af_med=float(np.nanmedian(af_all)),
)

# ---- write markdown ----
def pct(n, d): return f"{100*n/d:.1f}%"
L = []
L.append("# Production panel statistics (arch3, 231 founders) — for the paper\n")
L.append(f"Generated by `scripts/panel_stats_for_paper.py`. SV cutoff = length change > {SV_BP} bp.\n")
L.append(f"Panel: **{NF} founders** (78 cactus + 153 PanGenie assemblies), 5 chromosomes (TAIR10).\n")

L.append("\n## 1. File sizes & record counts\n")
L.append("| Chrom | VCF (MB) | records | V_pa set (MB) | K_pa set (MB) | k-mers | bubbles |")
L.append("|---|--:|--:|--:|--:|--:|--:|")
for c, d, *_ in per:
    L.append(f"| Chr{c} | {d['vcf_mb']:.0f} | {d['N']:,} | {d['vp_mb']:.0f} | {d['kp_mb']:.0f} | {d['nk']:,} | {d['nb']:,} |")
L.append(f"| **Total** | **{genome['vcf_mb']:.0f}** | **{genome['N']:,}** | "
         f"**{genome['vp_mb']:.0f}** | **{genome['kp_mb']:.0f}** | **{genome['nk']:,}** | **{genome['nb']:,}** |")
L.append("\nV_pa set = var_pa + var_called + meta. K_pa set = kmer_pa + meta. "
         "Matrix dims: V_pa = 231 × records; K_pa = 231 × k-mers (both int8 CSR).\n")

L.append("\n## 2. Variant classes (V_pa records)\n")
L.append("| Chrom | SNP | small indel ≤50bp | SV >50bp | (INS / DEL / MNP) |")
L.append("|---|--:|--:|--:|--:|")
for c, d, *_ in per:
    L.append(f"| Chr{c} | {d['n_snp']:,} | {d['n_small']:,} | {d['n_sv']:,} | "
             f"{d['n_ins']:,} / {d['n_del']:,} / {d['n_mnp']:,} |")
g = genome
L.append(f"| **Total** | **{g['n_snp']:,}** ({pct(g['n_snp'],g['N'])}) | "
         f"**{g['n_small']:,}** ({pct(g['n_small'],g['N'])}) | "
         f"**{g['n_sv']:,}** ({pct(g['n_sv'],g['N'])}) | "
         f"{g['n_ins']:,} / {g['n_del']:,} / {g['n_mnp']:,} |")

L.append("\n## 3. SV size spectrum (length change, bp)\n")
L.append("| 50–100 | 100–1k | 1k–10k | >10k | largest |")
L.append("|--:|--:|--:|--:|--:|")
L.append(f"| {g['sv_50_100']:,} | {g['sv_100_1k']:,} | {g['sv_1k_10k']:,} | {g['sv_10k']:,} | {g['sv_max']:,} bp |")

L.append("\n## 4. Allele-frequency / count spectrum (among called founders)\n")
L.append(f"- Singletons (AC=1): **{g['ac1']:,}** ({pct(g['ac1'],g['N'])})")
L.append(f"- AC ≤ 5: **{g['ac_le5']:,}** ({pct(g['ac_le5'],g['N'])})")
L.append(f"- Median ALT frequency among called founders: **{g['af_med']:.3f}**")

L.append("\n## 5. Missingness (V_pa called-mask)\n")
L.append(f"- Overall genotype missing rate: **{g['overall_miss']*100:.2f}%** "
         f"(1 − called entries / ({NF} × {g['N']:,}))")
L.append(f"- Per-record missing fraction: mean **{g['miss_mean']*100:.2f}%**, median **{g['miss_med']*100:.2f}%**")
L.append("| Chrom | mean missing | median missing |")
L.append("|---|--:|--:|")
for c, d, *_ in per:
    L.append(f"| Chr{c} | {d['miss_mean']*100:.2f}% | {d['miss_med']*100:.2f}% |")

L.append("\n## 6. K_pa (k-mer panel, k=31, filt2inv)\n")
L.append(f"- Total k-mers: **{g['nk']:,}**  | bubbles: **{g['nb']:,}**  | "
         f"mean k-mers/bubble: **{g['nk']/g['nb']:.1f}**")
L.append(f"- Founder–k-mer presences (matrix nnz): **{g['k_nnz']:,}**")
L.append(f"- Mean founders carrying a k-mer: **{g['kfc_mean']:.1f}** / {NF}  "
         f"(density {100*g['k_nnz']/(NF*g['nk']):.1f}%)")
L.append(f"- Founders-per-k-mer range: **{g['kfc_min']}–{g['kfc_max']}** "
         f"(filt2inv removes AC=1 and AC={NF} → min ≥2, max ≤{NF-1}) ✓")

L.append("\n## 7. Decomposition / density\n")
L.append(f"- Distinct positions: **{g['nk'] and sum(d['n_pos'] for _,d,*_ in per):,}** vs "
         f"{g['N']:,} records "
         f"(symbolic-ID decomposition expands co-located multiallelics)")

txt = "\n".join(L) + "\n"
open(f"{OUTDIR}/PANEL_STATS.md", "w").write(txt)
print("\n" + txt)
print(f"\nWrote {OUTDIR}/PANEL_STATS.md")
