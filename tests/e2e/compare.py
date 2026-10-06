"""Compare the Chr1 panel built end-to-end from kMate 0.1.3 with the production
(arch3) matrices. Records are matched on (pos, ref, alt) -- full alleles, never
pos alone (multi-allelic positions) -- and founders by name."""
import gzip, numpy as np, scipy.sparse as sp, pandas as pd

PROD = "/global/home/users/tbellg/kmate"
P_IDX = f"{PROD}/panel/pangenie_index/pang_135_haploid"
P_VAR = f"{PROD}/panel/arch3/chr1/var_pa_231_arch3_chr1"
P_KMR = f"{PROD}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"

def hdr(s): print(f"\n=== {s} ===", flush=True)

# graph samples are assembly IDs, production founders are 1001G ecotype IDs; this
# table maps 73 of the 78 production assembly founders independently of either panel.
MAP = pd.read_csv(f"{PROD}/panel/pangenie_genotyping/data/trueloo_founders.tsv", sep="\t", dtype=str)
ECO2ASM = dict(zip(MAP["ecotype"], MAP["assembly_id"]))

def founder_rows(ours, prod):
    """Row indices (ours, prod) for founders present in both, via ECO2ASM."""
    oi = {f: i for i, f in enumerate(ours)}
    pairs = [(oi[ECO2ASM[e]], j) for j, e in enumerate(prod) if e in ECO2ASM and ECO2ASM[e] in oi]
    return [a for a, _ in pairs], [b for _, b in pairs]

# ---- 1. k-mer index (all chromosomes) -------------------------------------------
hdr("1. k-mer index vs production pang_135_haploid")
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    a = gzip.open(f"index/ours_{c}_kmers.tsv.gz", "rt").read().splitlines()
    b = gzip.open(f"{P_IDX}/ours_{c}_kmers.tsv.gz", "rt").read().splitlines()
    same = sum(x == y for x, y in zip(a, b))
    print(f"{c}: ours {len(a):,} lines, production {len(b):,}; identical lines {same:,}"
          + ("  IDENTICAL" if a == b else ""))

# ---- 2. var_pa / var_called ------------------------------------------------------
hdr("2. var_pa / var_called vs production arch3 (shared founders, shared records)")
mo = np.load("var_pa/var_pa_Chr1.meta.npz", allow_pickle=True)
mp = np.load(f"{P_VAR}.meta.npz", allow_pickle=True)
Vo = sp.load_npz("var_pa/var_pa_Chr1.var_pa.npz").tocsc(); Co = sp.load_npz("var_pa/var_pa_Chr1.var_called.npz").tocsc()
Vp = sp.load_npz(f"{P_VAR}.var_pa.npz").tocsc();          Cp = sp.load_npz(f"{P_VAR}.var_called.npz").tocsc()
fo, fp = list(mo["founders"].astype(str)), list(mp["founders"].astype(str))
rfo, rfp = founder_rows(fo, fp)
print(f"founders: ours {len(fo)}, production {len(fp)}, mapped pairs {len(rfo)}")
ko = list(zip(mo["pos"].tolist(), mo["ref"].tolist(), mo["alt"].tolist()))
kp = list(zip(mp["pos"].tolist(), mp["ref"].tolist(), mp["alt"].tolist()))
dup_o, dup_p = len(ko) - len(set(ko)), len(kp) - len(set(kp))
print(f"records: ours {len(ko):,} (duplicate keys {dup_o:,}), production {len(kp):,} (duplicate keys {dup_p:,})")
ip = {k: i for i, k in enumerate(kp)}
pairs = [(i, ip[k]) for i, k in enumerate(ko) if k in ip]
print(f"shared records {len(pairs):,}; ours-only {len(ko)-len(pairs):,}; production-only {len(kp)-len(pairs):,}")
oi = np.array([a for a, _ in pairs]); pi = np.array([b for _, b in pairs])
for name, A, B in [("var_pa", Vo, Vp), ("var_called", Co, Cp)]:
    a = (A[:, oi][rfo, :].toarray() > 0); b = (B[:, pi][rfp, :].toarray() > 0)
    mism = int((a != b).sum())
    print(f"{name}: {a.size:,} cells compared, mismatches {mism:,} ({100*mism/a.size:.4f}%)"
          f"; mismatching records {int((a != b).any(0).sum()):,}")

# ---- 3. kmer_pa ------------------------------------------------------------------
hdr("3. kmer_pa vs production filt2inv (shared founders, shared k-mers)")
qo = np.load("kmer_pa/kmer_pa_Chr1.meta.npz", allow_pickle=True)
qp = np.load(f"{P_KMR}.meta.npz", allow_pickle=True)
Ko = sp.load_npz("kmer_pa/kmer_pa_Chr1.kmer_pa.npz").tocsc(); Kp = sp.load_npz(f"{P_KMR}.kmer_pa.npz").tocsc()
go, gp = list(qo["founders"].astype(str)), list(qp["founders"].astype(str))
so = pd.Index(qo["kmer_index"]); spi = pd.Index(qp["kmer_index"])
common = so.intersection(spi)
ro, rp = founder_rows(go, gp)
print(f"k-mers: ours {len(so):,}, production {len(spi):,}, shared {len(common):,}; founder pairs {len(ro)}")
ci_o = so.get_indexer(common); ci_p = spi.get_indexer(common)
mism = 0; cells = 0
for s in range(0, len(common), 2_000_000):
    a = Ko[:, ci_o[s:s+2_000_000]][ro, :].toarray() > 0
    b = Kp[:, ci_p[s:s+2_000_000]][rp, :].toarray() > 0
    mism += int((a != b).sum()); cells += a.size
print(f"kmer_pa: {cells:,} cells compared, mismatches {mism:,} ({100*mism/max(cells,1):.4f}%)")

# ---- 4. allele frequencies from the same pool ------------------------------------
hdr("4. SEEDMIX S1 Chr1 alt_freq: new 135-assembly panel vs production 231 panel")
eo = pd.read_csv("S1_new_Chr1.tsv", sep="\t"); ep = pd.read_csv("S1_prod_Chr1.tsv", sep="\t")
assert len(eo) == len(ko) and len(ep) == len(kp), "TSV rows must follow var_pa meta order"
x = eo["alt_freq"].values[oi]; y = ep["alt_freq"].values[pi]
ok = np.isfinite(x) & np.isfinite(y)
rl = mo["ref_len"][oi]; al = mo["alt_len"][oi]
cls = np.where(np.maximum(rl, al) >= 50, "SV", np.where((rl == 1) & (al == 1), "SNP", "indel"))
print(f"finite in both: {ok.sum():,} of {len(oi):,} shared records")
for c in ["ALL", "SNP", "indel", "SV"]:
    m = ok if c == "ALL" else ok & (cls == c)
    if m.sum() < 2: continue
    d = x[m] - y[m]
    print(f"{c:5s} n={m.sum():>9,}  r={np.corrcoef(x[m], y[m])[0,1]:.4f}  MAE={np.abs(d).mean():.4f}"
          f"  |diff|>0.1: {100*np.mean(np.abs(d) > 0.1):.2f}%")
