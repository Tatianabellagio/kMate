#!/usr/bin/env python
"""Prepare gggenomes input tables for a 'pangenome-style' CAM5 haplotype figure.

Bins (rows), top -> bottom:
  3 isoform rows  (AT2G27030.1/.2/.3)  -> gene model drawn as exon/CDS boxes
  227 founder rows (ordered cold->hot by home-climate bio1) -> ALT variants as
       SNP/indel ticks.

All bins share the same genomic window so columns line up. Mirrors the JBrowse
MAF sketch but as a clean, fully-controlled gggenomes figure.

Emits (to this dir):
  gg_seqs.tsv   seq_id, bin_id, length, kind, bio1   (3 isoforms + 227 founders)
  gg_model.tsv  seq_id(isoform), start, end, type     (exon/CDS/UTR boxes)
  gg_feats.tsv  seq_id(founder), start, end, cls, kendall_p, nlp  (ALT ticks)
"""
import numpy as np, pandas as pd, pysam

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/cam5_replication/pang_cam5"
WIN_LO, WIN_HI = 11531800, 11534333          # region (gene span + small 5' pad)
LEN = WIN_HI - WIN_LO
ISOFORMS = ["AT2G27030.1", "AT2G27030.2", "AT2G27030.3"]

# ---- founder climate (bio1 at origin) ----
bio1 = pd.read_csv("/tmp/ecotype_bio1.tsv", sep="\t").set_index("ecotype")["bio1"].to_dict()

# ---- genotypes + per-variant GEA stats ----
vcf = pysam.VariantFile(f"{HERE}/cam5_region.vcf.gz")
samp = [int(s) for s in vcf.header.samples]
ann = pd.read_csv(f"{HERE}/cam5_variants_gen9.tsv", sep="\t").drop_duplicates("pos", keep="first").set_index("pos")

feat_rows = []
for r in vcf.fetch("Chr2", WIN_LO, WIN_HI):
    cls = "snp" if (len(r.ref) == 1 and len(r.alts[0]) == 1) else "indel"
    kp = float(ann.loc[r.pos, "kendall_p"]) if r.pos in ann.index else np.nan
    nlp = -np.log10(kp) if (kp == kp and kp > 0) else 0.0
    x = r.pos - WIN_LO
    for s in vcf.header.samples:
        if r.samples[s].get("GT", (None,))[0] == 1:        # ALT carrier -> tick
            feat_rows.append((str(s), x, x + max(len(r.ref), 1), cls, kp, nlp))
feats = pd.DataFrame(feat_rows, columns=["seq_id", "start", "end", "cls", "kendall_p", "nlp"])
feats.to_csv(f"{HERE}/gg_feats.tsv", sep="\t", index=False)

# ---- 3' climate haplotype: carriers + zoom matrix (window matches tube-map) ----
HW_LO, HW_HI, MINN, LEAD = 11533880, 11534160, 3, 11533967
hrec, Gh = [], []
for r in vcf.fetch("Chr2", HW_LO, HW_HI):
    Gh.append([(r.samples[s].get("GT", (None,))[0]
                if r.samples[s].get("GT", (None,))[0] is not None else -1) for s in vcf.header.samples])
    hrec.append((r.pos, "snp" if (len(r.ref) == 1 and len(r.alts[0]) == 1) else "indel"))
Gh = np.array(Gh).T
H, cnt = np.unique(Gh, axis=0, return_counts=True)
keep = cnt >= MINN; H, cnt = H[keep], cnt[keep]
clim = H[int(np.argmax((H == 1).sum(1)))]                  # most-alt class = climate hap
carriers = {samp[i] for i in range(len(samp)) if np.array_equal(Gh[i], clim)}

# zoom matrix (long): founder x site genotype over the 21 hap-window sites
mrows = []
for j, (pos, cls) in enumerate(hrec):
    for i, e in enumerate(samp):
        g = Gh[i, j]
        mrows.append((str(e), j, pos, cls, {0: "ref", 1: "alt", -1: "miss"}[g],
                      int(pos == LEAD)))
mat = pd.DataFrame(mrows, columns=["seq_id", "site_idx", "pos", "cls", "geno", "is_lead"])
mat.to_csv(f"{HERE}/gg_hapmatrix.tsv", sep="\t", index=False)

# ---- seqs: 3 isoform rows, then founders cold->hot (+ carrier flag, + hap-cluster order) ----
order = sorted(samp, key=lambda e: (bio1.get(e, 99), e))
alt_in_win = {samp[i]: int((Gh[i] == 1).sum()) for i in range(len(samp))}
seq_rows = [(iso, iso, LEN, "isoform", np.nan, 0, -1) for iso in ISOFORMS]
seq_rows += [(str(e), str(e), LEN, "founder", bio1.get(e, np.nan),
              int(e in carriers), alt_in_win[e]) for e in order]
seqs = pd.DataFrame(seq_rows, columns=["seq_id", "bin_id", "length", "kind", "bio1",
                                       "hap_carrier", "alt_in_win"])
seqs.to_csv(f"{HERE}/gg_seqs.tsv", sep="\t", index=False)

# ---- gene model: exon/CDS/UTR boxes per isoform ----
feat = pd.read_csv(f"{HERE}/cam5_features.tsv", sep="\t",
                   names=["type", "start", "end", "strand", "attr"])
def iso(a):
    for kv in str(a).split(";"):
        if kv.startswith("Parent="): return kv.split("=")[1].split(",")[0]
    return None
g = feat[feat.type.isin(["CDS", "five_prime_UTR", "three_prime_UTR"])].copy()
g["seq_id"] = g.attr.map(iso)
g = g[g.seq_id.isin(ISOFORMS)].copy()
g["start"] = g.start - WIN_LO
g["end"] = g.end - WIN_LO
g["type"] = g.type.replace({"five_prime_UTR": "UTR", "three_prime_UTR": "UTR"})
g[["seq_id", "start", "end", "type"]].to_csv(f"{HERE}/gg_model.tsv", sep="\t", index=False)

print(f"founders with bio1: {sum(e in bio1 for e in samp)}/{len(samp)}  (range "
      f"{min(bio1.values()):.1f}..{max(bio1.values()):.1f} C)")
print(f"feats (ALT ticks): {len(feats)}  snp={sum(feats.cls=='snp')} indel={sum(feats.cls=='indel')}")
# ---- collapse founders to UNIQUE HAPLOTYPES (by drawn ALT pattern over the gene window) ----
wsites, Gw = [], []
for r in vcf.fetch("Chr2", WIN_LO, WIN_HI):
    rl, al = len(r.ref), len(r.alts[0]); p0 = r.pos - WIN_LO
    if rl == 1 and al == 1:                       # SNP
        vtype, size, s0, e0 = "snp", 0, p0, p0 + 1
    elif rl == al:                                # equal-length multi-base substitution
        vtype, size, s0, e0 = "mnp", 0, p0, p0 + rl
    elif rl > al:                                 # deletion: bar spans the deleted bases
        size = rl - al; vtype, s0, e0 = "del", p0 + al, p0 + al + size
    else:                                         # insertion: point, sized by inserted bp
        size = al - rl; vtype, s0, e0 = "ins", p0 + rl, p0 + rl
    wsites.append((s0, e0, vtype, size))
    Gw.append([1 if r.samples[s].get("GT", (None,))[0] == 1 else 0 for s in vcf.header.samples])
Gw = np.array(Gw).T                                    # founders x sites, 1=ALT drawn
pat, inv = np.unique(Gw, axis=0, return_inverse=True)  # unique ALT patterns
nhap = pat.shape[0]
b1 = np.array([bio1.get(e, np.nan) for e in samp])

hseq, hfeat = [], []
for k in range(nhap):
    members = np.where(inv == k)[0]
    n_eco = len(members)
    mb = np.nanmean(b1[members]) if np.any(~np.isnan(b1[members])) else np.nan
    hid = f"hap{k:03d}"
    hseq.append((hid, hid, LEN, "hap", mb, n_eco))
    for j in np.where(pat[k] == 1)[0]:
        s0, e0, vtype, size = wsites[j]
        hfeat.append((hid, s0, e0, vtype, size))
hs = pd.DataFrame(hseq, columns=["seq_id", "bin_id", "length", "kind", "mean_bio1", "n_eco"])
# order: isoforms first, then haps cold->hot by mean climate (nan last)
hs = hs.sort_values("mean_bio1", kind="mergesort", na_position="last")
iso = pd.DataFrame([(i, i, LEN, "isoform", np.nan, 0) for i in ISOFORMS], columns=hs.columns)
pd.concat([iso, hs]).to_csv(f"{HERE}/gg_hapseqs.tsv", sep="\t", index=False)
pd.DataFrame(hfeat, columns=["seq_id", "start", "end", "vtype", "size"]).to_csv(
    f"{HERE}/gg_hapfeats.tsv", sep="\t", index=False)
pd.DataFrame({"name": ["bio1_mean"], "x": [round(float(np.nanmean(b1)), 3)]}).to_csv(
    f"{HERE}/gg_hapmeta.tsv", sep="\t", index=False)

print(f"gene-model boxes: {len(g)}  isoforms: {sorted(g.seq_id.unique())}")
print(f"climate-hap carriers: n={len(carriers)}  | zoom matrix sites: {len(hrec)} "
      f"({sum(c=='snp' for _,c in hrec)} SNP, {sum(c=='indel' for _,c in hrec)} indel)")
print(f"unique haplotypes: {nhap} from {len(samp)} founders "
      f"(largest group n={int(np.bincount(inv).max())})")
print("wrote gg_seqs.tsv, gg_feats.tsv, gg_model.tsv, gg_hapmatrix.tsv, gg_hapseqs.tsv, gg_hapfeats.tsv")
