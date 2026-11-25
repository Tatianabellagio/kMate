#!/usr/bin/env python
"""Build + execute notebooks/gea_hit_landscape.ipynb: what the non-SNP GEA hits are and where
they land -- SV hits first, small indels as the comparison -- against everything the GEA tested.

Reads the GEA's own outputs, via results/gea_hit_landscape.csv.gz (gea_hit_landscape_data.py), plus
(section 1b) the SV repeat classification in r1_sv_negative_selection/results/sv_adaptive/graffite/. No chart titles (project convention). Env: basic.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
GS = os.path.abspath(os.path.join(HERE, "..", ".."))
NBDIR = f"{GS}/notebooks"
OUT = f"{NBDIR}/gea_hit_landscape.ipynb"

C = []
md = lambda s: C.append(new_markdown_cell(s))     # noqa: E731
co = lambda s: C.append(new_code_cell(s))         # noqa: E731

md(r"""# Where the non-SNP GEA hits land — SV hits, and small indels for comparison

**What this is.** A description of the GEA's own hits: what kind of variant they are, and where
in the genome they sit, each set against **everything the GEA tested** (the right background:
only variants common enough to test can become hits).

**Rebuilt 2026-09-30 on the re-run GEA.** Two things changed, and both move the numbers here:
a **founder call-rate floor of ≥ 0.50** on top of MAF ≥ 0.05 (a variant called in only the 78
long-read founders is a different estimand, not a noisier one — see `notebooks/flc_case.ipynb`),
and **SVs are now tested as truvari events**, not as arch3 allele records. The previous version
of this notebook, with neither change, is at `notebooks/old_lfmm/gea_hit_landscape.ipynb`.

**Hits.** The GEA's own call: raw LFMM p-values (gen 9, 22 climate axes = bio1–19 + pc1–3), per
class × axis MAF > 0.05 and call rate ≥ 0.50, Bonferroni at 0.05 / units tested; a unit is a hit
if it clears that bar on **any** axis in any of the three class scans (SV, small indel, and the
two pooled). Nothing is re-estimated here.

**Classes and types.** *SV* = the unit is an SV event (|Δlen| > 50 bp); *small indel* otherwise.
*Deletion* = REF longer than ALT, *insertion* = ALT longer, both relative to Col-0 (not to the
ancestral state). Equal-length MNPs are neither and are set aside, counted in section 0.

**Regions.** The same classifier the candidate pool uses
(`genes/dissection/screen_sig_blocks.classify`), run on every tested record
(`annotate_gea_universe.py`, on the event's representative record — the region call depends on
the REF span only). First match wins, gene body before promoter:
CDS > UTR / non-coding exon > intron > **promoter (≤ 1 kb upstream of a TSS, strand-aware, outside
any gene body)** > reference TE > intergenic (proximal) > gene desert.

**Unit.** One row per **GEA unit**: a truvari event for SVs, a record for short indels and MNPs
(which are one unit per record either way). That is the unit the burden and purging work already
used (`sv_purging_fixed.ipynb`), and the unit the SV literature uses. Units carry an integer
`unit` id, so the scans are matched on that id and never on coordinates — 2.14% of arch3
positions carry more than one record, and a coordinate join there pairs the wrong allele.
Section 9 shows what the record-level arm would have given instead, in both directions.

**Enrichment.** For a region *R*: the odds that a variant in *R* is a hit, over the odds for a
variant outside *R*, Mantel–Haenszel over size-bin × MAF-bin strata (hits skew common and
towards particular sizes, so an unadjusted share would mislead), with 95% CIs from resampling
the GEA's LD blocks — hits cluster in blocks, so records are not independent.""")

co(r"""import os, sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
CONV = os.path.abspath("../genes/convergence"); sys.path.insert(0, CONV)
import plot_theme as TH
import gea_region_by_axis as GR          # or_with_ci (MH + block bootstrap), size / MAF bins
TH.apply()
PLOTS = os.path.abspath("plots/gea_hit_landscape"); os.makedirs(PLOTS, exist_ok=True)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 400)

K4 = ["chrom", "pos", "ref_len", "alt_len"]
U = pd.read_csv(f"{CONV}/results/gea_hit_landscape.csv.gz")
n_raw = len(U)
# split records -> one; prefer the record the GEA called (its min p can be the higher of the two:
# the call may come from the SV scan, the min p here is from the pooled scan)
U = U.sort_values(["hit", "min_p"], ascending=[False, True]).drop_duplicates(K4).reset_index(drop=True)
REGION = {"1_CDS": "CDS", "2_UTR": "UTR / ncRNA exon", "4_exon_noncoding": "UTR / ncRNA exon",
          "5_intron": "intron", "3_promoter": "promoter (≤1 kb)", "6_TE": "TE",
          "7_proximal_intergenic": "intergenic", "8_gene_desert": "gene desert"}
RORDER = ["CDS", "UTR / ncRNA exon", "intron", "promoter (≤1 kb)", "TE", "intergenic", "gene desert"]
U["region"] = U.tier_1kb.map(REGION)
U["stratum"] = (pd.cut(U["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                + pd.cut(U.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
# SV colour convention of notebooks/sv_burden_final_models.ipynb: insertions forest green,
# deletions orange, and NO blue or red for variant classes -- those are the garden-temperature
# scale, used for direction in section 7. (sv_main_fig1.ipynb swaps green and orange; the
# final-models convention is the one followed here.)
KIND = {"deletion": "#F07A2E", "insertion": "#4E6B3A"}
REF = TH.MUTED                                           # "all tested" is context, drawn grey
rng = np.random.default_rng(1)
print(f"{n_raw:,} tested records -> {len(U):,} after collapsing split records; "
      f"{int(U.hit.sum()):,} GEA hits (the candidate pool has 2,938)")""")

md(r"""## 0 · The numbers

Tested records and GEA hits by class and type; *per 10k* is the hit rate, *blocks* the number of
distinct LD blocks the hits fall in (hits cluster, so blocks are the conservative count).""")

co(r"""T = U.groupby(["cls", "kind"]).agg(tested=("hit", "size"), hits=("hit", "sum"))
T["hits per 10k tested"] = (1e4 * T.hits / T.tested).round(1)
T["hit blocks"] = U[U.hit].groupby(["cls", "kind"]).block.nunique()
T["share of class hits"] = (T.hits / T.groupby(level=0).hits.transform("sum")).round(3)
print(T.fillna(0).astype({"hit blocks": int}).to_string())
sv = U[(U.cls == "sv") & U.kind.isin(["deletion", "insertion"])].copy()
si = U[(U.cls == "smallindel") & U.kind.isin(["deletion", "insertion"])].copy()
print(f"\nSV hits: {int(sv.hit.sum())} ({int(sv[sv.hit].kind.eq('deletion').sum())} deletions, "
      f"{int(sv[sv.hit].kind.eq('insertion').sum())} insertions) in {sv[sv.hit].block.nunique()} blocks")""")

md(r"""## 1 · What the SV hits are — size and frequency

**A–B** size of SV deletions and insertions: GEA hits (bars) against every SV of that type the
GEA tested (grey line), each as a share of its own set so the shapes compare. **C** minor allele
frequency: hits by type against all tested SVs. The enrichment tests below stratify on
frequency and size anyway, so a difference in either cannot masquerade as a regional one.""")

co(r"""fig, axs = plt.subplots(1, 3, figsize=(14, 3.7), gridspec_kw=dict(wspace=0.32))
bins = np.linspace(np.log10(50), np.log10(sv["size"].max()) + 0.01, 30)
tick = [50, 100, 300, 1000, 3000, 10000, 30000]
for ax, k, L in zip(axs[:2], ["deletion", "insertion"], "AB"):
    t = sv[sv.kind == k]; h = t[t.hit]
    ax.hist(np.log10(t["size"]), bins=bins, weights=np.full(len(t), 1 / len(t)), histtype="step",
            color=REF, lw=1.3, label=f"all tested ({len(t):,})")
    ax.hist(np.log10(h["size"]), bins=bins, weights=np.full(len(h), 1 / len(h)), color=KIND[k],
            alpha=0.85, rwidth=0.88, label=f"GEA hits ({len(h):,})")
    ax.set_xticks(np.log10(tick)); ax.set_xticklabels([f"{v:,}" for v in tick])
    ax.set_xlim(bins[0], min(bins[-1], np.log10(40000)))
    ax.set_xlabel(f"{k} size (bp)"); ax.set_ylabel("share of variants")
    TH.grid_only(ax, "y"); TH.panel(ax, L, f"SV {k}s", y=1.12); ax.legend(loc="upper right")
ax = axs[2]; mb = np.linspace(0.05, 0.5, 19)
ax.hist(sv.MAF, bins=mb, weights=np.full(len(sv), 1 / len(sv)), histtype="step", color=REF, lw=1.3,
        label=f"all tested SVs ({len(sv):,})")
for k in KIND:
    h = sv[sv.hit & (sv.kind == k)]
    ax.hist(h.MAF, bins=mb, weights=np.full(len(h), 1 / len(h)), histtype="step", color=KIND[k],
            lw=2, label=f"hits, {k}s ({len(h)})")
ax.set_xlabel("minor allele frequency"); ax.set_ylabel("share of variants")
TH.grid_only(ax, "y"); TH.panel(ax, "C", "frequency", y=1.12); ax.legend(loc="upper right")
plt.savefig(f"{PLOTS}/sv_hits_size_maf.png", dpi=170, bbox_inches="tight"); plt.show()
for k in KIND:
    t = sv[sv.kind == k]
    print(f"{k:9s} median size: hits {t[t.hit]['size'].median():,.0f} bp, tested {t['size'].median():,.0f} bp | "
          f"median MAF: hits {t[t.hit].MAF.median():.3f}, tested {t.MAF.median():.3f}")""")

md(r"""## 1b · What the SV hits carry — TE content

Each SV unit's repeat class, from the classification of every SV record with GraffiTE's masking protocol
(`r1_sv_negative_selection/TE_CLASSIFICATION_PIPELINE.md`). A deletion's class describes the reference
sequence it removes, an insertion's the sequence it adds. A unit is **at least half TE** if its event is
`te_full` (≥ 80% of the allele masked as TE) or `te_most` (≥ 50%). The event takes the class of its
first record. Units are matched to events through integer ids (unit → representative af_store record →
panel record → Truvari event), never coordinates.

**A–B** share at least half TE among all tested SVs of that type and among the GEA hits. The odds ratio
is the same Mantel–Haenszel / block-bootstrap statistic as section 2, over finer strata here:
size (50–100–300–1k–3k–10k–∞ bp) × MAF × **replacement status**. TE content rises steeply with size,
and replacement events (HPRC `COMPLEX`, or a swapped base) are much less often whole TEs at the same
size, so either could fake an enrichment.""")

co(r"""R1 = os.path.abspath("../r1_sv_negative_selection/results/sv_adaptive")
PANEL = os.path.abspath("../../../panel/arch3")
AFI = []                                   # af_store non-SNP index -> (chrom, panel record)
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    m = np.load(f"{PANEL}/{c.lower()}/var_pa_231_arch3_{c.lower()}.meta.npz", allow_pickle=True)
    rl, al = m["ref_len"].astype(int), m["alt_len"].astype(int)
    AFI.append(pd.DataFrame(dict(chrom=c, rec=np.where(~((rl == 1) & (al == 1)))[0], rl=rl[~((rl == 1) & (al == 1))])))
AFI = pd.concat(AFI, ignore_index=True)
EU = np.load(f"{R1}/nonsnp_event_units.npz")
SVT = sv.copy()
rep = AFI.iloc[EU["rep"][SVT.unit.to_numpy()]]
assert (rep.rl.to_numpy() == SVT.ref_len.to_numpy()).all(), "unit -> record mapping is off"
RCK = pd.read_csv(f"{R1}/graffite/sv_repeat_class_records.csv", usecols=["key", "event_key"]).set_index("key")
EVC = pd.read_csv(f"{R1}/graffite/sv_repeat_class_events.csv").set_index("event_key")
ek = (rep.chrom + "|" + rep.rec.astype(str)).map(RCK.event_key).to_numpy()
for col in ("cargo", "any_replacement", "top_superfamily"):
    SVT[col] = EVC[col].reindex(ek).to_numpy()
assert SVT.cargo.notna().all(), "SV units without a repeat class"
SVT["te"] = SVT.cargo.isin(["te_full", "te_most"])
TSZ = [50, 100, 300, 1000, 3000, 10000, np.inf]
rng_te = np.random.default_rng(11)          # own stream: leaves the shared `rng` (sections 2-3) untouched
rows = []
fig, axs = plt.subplots(1, 2, figsize=(10, 3.7), gridspec_kw=dict(wspace=0.38))
for ax, k, L in zip(axs, ["deletion", "insertion"], "AB"):
    d = SVT[SVT.kind == k].copy()
    d["stratum"] = (pd.cut(d["size"], TSZ, labels=False).astype(str) + "|"
                    + pd.cut(d.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str) + "|"
                    + d.any_replacement.astype(str))
    o, lo, hi, p = GR.or_with_ci(d, "hit", "te", rng_te)
    h = d[d.hit]
    ax.bar([0, 1], [100 * d.te.mean(), 100 * h.te.mean()], color=[REF, KIND[k]], width=0.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"all tested\n({len(d):,})", f"GEA hits\n({len(h)})"])
    ax.set_ylabel("% at least half TE"); ax.set_ylim(0, 145 * max(d.te.mean(), h.te.mean()))
    ax.annotate(f"OR {o:.2f} [{lo:.2f}, {hi:.2f}], p = {p:.2g}", xy=(0.03, 0.97), xycoords="axes fraction",
                va="top", fontsize=8)
    TH.grid_only(ax, "y"); TH.panel(ax, L, f"SV {k}s", y=1.12)
    rows.append(dict(kind=k, tested=len(d), hits=len(h), pct_te_tested=100 * d.te.mean(),
                     pct_te_hits=100 * h.te.mean(), OR=o, ci_lo=lo, ci_hi=hi, p=p,
                     med_size_hits=h["size"].median(), med_size_tested=d["size"].median(),
                     pct_repl_hits=100 * h.any_replacement.mean(), pct_repl_tested=100 * d.any_replacement.mean()))
plt.savefig(f"{PLOTS}/sv_hits_te_content.png", dpi=170, bbox_inches="tight"); plt.show()
print("TE CONTENT OF THE SV HITS  (at least half TE = te_full + te_most)\n")
print(pd.DataFrame(rows).round(3).to_string(index=False))
for k in KIND:
    h = SVT[SVT.hit & (SVT.kind == k) & SVT.te]
    print(f"\n{k} hits at least half TE: {len(h)}; superfamilies {h.top_superfamily.value_counts().head(5).to_dict()}")""")

md(r"""**Reading it.** **Insertion hits are TE-enriched; deletion hits are not.** 53.6% of the 151
insertion hits are at least half TE against 44.3% of all tested insertions (OR 1.67, 95% CI 1.16–2.36,
block-bootstrap p = 0.008, within size × MAF × replacement strata). For the 131 deletion hits it is
42.0% against 38.0% (OR 1.16, CI 0.76–1.68, p = 0.45). Hits and background have about the same
median size and the same replacement share, so neither explains the insertion result. RC/Helitron is the
commonest superfamily among the TE hits of both types, but no family dominates
(`notebooks/gea_te_cargo.ipynb` tests families one by one).

Two limits. "Insertion" and "deletion" are relative to Col-0, not ancestral state, so an insertion hit can
be an older element lost from Col-0's lineage. And "tested" is already TE-poorer than the panel as a
whole, because the MAF and call-rate floors remove rare, TE-rich insertions. The enrichment is relative to
what the GEA could test.""")

md(r"""## 2 · Where the SV hits land

**A–B** share of SV hits in each region (bars), deletions and insertions separately, with the
share among **all tested SVs of the same type** as a grey tick — a bar past its tick means that
region holds more of the hits than of the tested variants. *n* in each label is the hit count.
**C** the same comparison as an enrichment: Mantel–Haenszel odds ratio over size × MAF strata,
95% block-bootstrap CI, log scale; filled = CI excludes 1.""")

co(r"""def region_stats(D):
    rows = []
    for k in KIND:
        d = D[D.kind == k]
        for r in RORDER:
            e = d.assign(_x=d.region.eq(r))
            if e._x.sum() == 0:
                continue
            o, lo, hi, p = GR.or_with_ci(e, "hit", "_x", rng)
            rows.append(dict(kind=k, region=r, hits=int((e.hit & e._x).sum()), tested=int(e._x.sum()),
                             hit_share=e[e.hit]._x.mean(), tested_share=e._x.mean(),
                             OR=o, lo=lo, hi=hi, p=p))
    return pd.DataFrame(rows)


XL = (0.25, 4.0)                      # odds-ratio axis; a CI past it is drawn to the edge


def region_figure(D, S, cls_label, fname):
    fig = plt.figure(figsize=(14, 4.3)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.95)
    y = np.arange(len(RORDER))
    for i, (k, L) in enumerate(zip(KIND, "AB")):
        ax = fig.add_subplot(gs[0, i]); s = S[S.kind == k].set_index("region").reindex(RORDER)
        ax.barh(y, s.hit_share, color=KIND[k], height=0.62, alpha=0.9)
        ax.scatter(s.tested_share, y, marker="|", s=260, color=REF, lw=2, zorder=3,
                   label="share among all tested")
        ax.set_yticks(y); ax.set_yticklabels([f"{r}   n={int(n) if n == n else 0}" for r, n in zip(RORDER, s.hits)])
        ax.invert_yaxis(); ax.set_xlabel(f"share of {cls_label} {k} hits")
        ax.set_xlim(0, 1.12 * np.nanmax([s.hit_share.max(), s.tested_share.max()]))
        TH.grid_only(ax, "x"); TH.panel(ax, L, f"{cls_label} {k}s", y=1.13)
        if i == 0:
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0))
    ax = fig.add_subplot(gs[0, 2])
    for j, k in enumerate(KIND):
        s = S[S.kind == k].set_index("region").reindex(RORDER); off = -0.17 if j == 0 else 0.17
        ok = s.OR.notna() & (s.OR > 0)
        yy = y[ok.values] + off; ss = s[ok].copy()
        ss["lo"] = ss.lo.clip(lower=XL[0]); ss["hi"] = ss.hi.clip(upper=XL[1])   # CIs run to the edge
        sig = (ss.lo > 1) | (ss.hi < 1)
        ax.errorbar(ss.OR, yy, xerr=[ss.OR - ss.lo, ss.hi - ss.OR], fmt="none", ecolor=KIND[k], lw=1.4)
        ax.scatter(ss.OR[sig], yy[sig.values], color=KIND[k], s=46, zorder=3, label=f"{k}s")
        ax.scatter(ss.OR[~sig], yy[~sig.values], facecolor="white", edgecolor=KIND[k], s=46, lw=1.4, zorder=3)
    ax.axvline(1, color=TH.MUTED, lw=1)
    ax.set_xscale("log"); ax.set_xlim(*XL); ax.minorticks_off()
    ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.set_yticks(y); ax.set_yticklabels(RORDER); ax.invert_yaxis()
    ax.set_xlabel("odds ratio, hit vs not")
    TH.grid_only(ax, "x"); TH.panel(ax, "C", "enrichment", y=1.13)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
    plt.savefig(f"{PLOTS}/{fname}", dpi=170, bbox_inches="tight"); plt.show()


S_sv = region_stats(sv)
region_figure(sv, S_sv, "SV", "sv_hits_by_region.png")
print(S_sv.assign(hit_share=S_sv.hit_share.round(3), tested_share=S_sv.tested_share.round(3),
                  OR=S_sv.OR.round(2), lo=S_sv.lo.round(2), hi=S_sv.hi.round(2), p=S_sv.p.round(3)).to_string(index=False))""")

md(r"""## 2b · Why promoters are already common among the tested variants

The grey ticks in section 2 are the share among **tested** variants, not the genome's share —
and among tested SV deletions promoters are about a quarter, far above the ~16% of the genome's
bases they cover. Two checks on where that comes from.

**A — is the promoter call right?** Every variant the classifier calls promoter is checked against
an *independent* signed distance to the nearest TSS (a different gene loader): a promoter call
should sit 0–1 kb upstream of a TSS, and a non-genic variant *not* called promoter should not.

**B — what would chance give?** A variant is an interval, not a base: a 2 kb deletion dropped at
random overlaps a gene body — which the classifier ranks first — far more often than a single base,
so a per-base genome share is the wrong yardstick for SVs. The baseline here is random positions
with the **same size distribution** as the tested variants of each class × type, run through the
**same classifier** (`gea_region_baseline.py`). **Left:** tested share ÷ that random share, per
region. **Right:** the same among variants *outside* gene bodies only, which removes the gene-body
effect and asks whether promoters hold more than their share of the non-genic variants.
Uniform placement includes sequence where nothing can be called (centromeric repeats, gaps),
which mainly moves the TE comparison.""")

co(r"""CA = pd.read_csv(f"{CONV}/results/gea_region_checkA.csv")
print("A  share of variants sitting 0-1 kb upstream of the nearest TSS (independent gene loader):")
print(CA.assign(share_0_1kb_upstream=CA.share_0_1kb_upstream.round(3)).to_string(index=False))
B = pd.read_csv(f"{CONV}/results/gea_region_baseline.csv")
SER = [("sv", "deletion"), ("sv", "insertion"), ("smallindel", "deletion"), ("smallindel", "insertion")]
MK = {"sv": "o", "smallindel": "s"}; LAB = {"sv": "SV", "smallindel": "small-indel"}
fig, axs = plt.subplots(1, 2, figsize=(12.5, 3.8), gridspec_kw=dict(wspace=0.45))
for ax, col, regs, L, txt in [(axs[0], "ratio", ["gene body", "promoter", "TE", "intergenic"], "A", "all tested variants"),
                              (axs[1], "ratio_outside_genes", ["promoter", "TE", "intergenic"], "B", "outside gene bodies")]:
    y = np.arange(len(regs))
    for j, (c, k) in enumerate(SER):
        d = B[(B.cls == c) & (B.kind == k)].set_index("region").reindex(regs)
        ax.scatter(np.log2(d[col]), y + (j - 1.5) * 0.14, marker=MK[c], s=52, color=KIND[k],
                   edgecolor="white", linewidth=0.6, zorder=3, label=f"{LAB[c]} {k}s")
    ax.axvline(0, color=TH.MUTED, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(regs); ax.invert_yaxis()
    ax.set_xticks([-2, -1, 0, 1, 2]); ax.set_xticklabels(["¼×", "½×", "1×", "2×", "4×"]); ax.set_xlim(-2, 2)
    ax.set_xlabel("tested share ÷ size-matched random share")
    TH.grid_only(ax, "x"); TH.panel(ax, L, txt, y=1.13)
axs[1].legend(loc="center right")          # the 2x-4x band is empty in both panels
plt.savefig(f"{PLOTS}/tested_vs_size_matched_random.png", dpi=170, bbox_inches="tight"); plt.show()
print(B.assign(**{c: B[c].round(3) for c in ["random_share", "tested_share", "ratio", "ratio_outside_genes"]}).to_string(index=False))
gb = B[B.region == "gene body"]; po = B[B.region == "promoter"]
print(f"\ngene bodies hold {gb.ratio.min():.2f}-{gb.ratio.max():.2f}x their size-matched share; "
      f"outside genes, promoters hold {po.ratio_outside_genes.min():.2f}-{po.ratio_outside_genes.max():.2f}x theirs")""")

md(r"""**Reading it.** The classifier's promoter calls check out (A). The large "promoter excess"
among tested variants is mostly a **gene-body deficit**: common indels and SVs are strongly
depleted from genes — purifying selection keeps variants in genes from reaching the 5% frequency
the GEA needs — so every region outside genes is inflated against the genome. Outside gene
bodies (B), promoters hold only modestly more than their share: a little for SVs, more for small
indels, where AT-rich poly-A/T runs typical of Arabidopsis promoters are a plausible mutational
source (not tested here). None of this enters the climate test in section 2, which compares hits
against these same tested variants and so holds all of it fixed.""")

md(r"""## 3 · The same for small indels

Small indels are ~12× more numerous among the hits, so this is where the estimates are tight.
Same layout and the same statistic as section 2.""")

co(r"""S_si = region_stats(si)
region_figure(si, S_si, "small-indel", "smallindel_hits_by_region.png")
print(S_si.assign(hit_share=S_si.hit_share.round(3), tested_share=S_si.tested_share.round(3),
                  OR=S_si.OR.round(2), lo=S_si.lo.round(2), hi=S_si.hi.round(2), p=S_si.p.round(3)).to_string(index=False))
pd.concat([S_sv.assign(cls="sv"), S_si.assign(cls="smallindel")]).to_csv(
    f"{CONV}/results/gea_hits_region_enrichment.csv", index=False)""")

md(r"""## 4 · Across the climate axes

For every axis with at least 20 hits (SVs and small indels together, deletions and insertions
together — split further, most cells are empty): each region's share of **that axis's** hits
divided by its share among all tested variants, on a log₂ scale. Green = more hits than the
region's size predicts, purple = fewer, grey = as expected. Cells where fewer than 3 hits are
expected are left blank: a ratio there is noise. The count after each axis name is its number
of hits. The axes are correlated with each other, so rows are not independent.""")

co(r"""A = U[U.kind.isin(["deletion", "insertion"])]
exp_share = A.region.value_counts(normalize=True).reindex(RORDER)
rows, counts = [], {}
for ax_ in GR.AXES:
    h = A[A.hit_axes.fillna("").str.split(",").map(lambda s: ax_ in s)]
    if len(h) < 20:
        continue
    counts[ax_] = len(h)
    obs = h.region.value_counts().reindex(RORDER).fillna(0); exp = exp_share * len(h)
    rows.append(np.where(exp >= 3, np.log2((obs + 0.5) / (exp + 0.5)), np.nan))
M = np.array(rows); names = list(counts)
cmap = LinearSegmentedColormap.from_list("enr", ["#6A4C93", "#E6E6E6", "#1B5E20"])
fig, ax = plt.subplots(figsize=(8.8, 0.34 * len(names) + 1.6))
im = ax.imshow(M, cmap=cmap, vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(len(RORDER))); ax.set_xticklabels(RORDER, rotation=30, ha="right")
ax.set_yticks(range(len(names))); ax.set_yticklabels([f"{a}  ({counts[a]})" for a in names])
ax.grid(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("hits ÷ expected (log₂)", color=TH.LABEL)
cb.outline.set_visible(False)
plt.savefig(f"{PLOTS}/hits_by_region_by_axis.png", dpi=170, bbox_inches="tight"); plt.show()
print("regions consistently above expectation (log2 > 0 on most axes):")
print(pd.Series((M > 0).mean(0), index=RORDER).round(2).to_string())""")

md(r"""## 5 · Distance to the nearest transcription start site

Signed distance from each variant to the nearest TSS, in the gene's orientation (negative =
upstream). Only the ±3 kb window is drawn, but each bar is a share of the set's **total**, so the
curves also show how much of each set lies near a TSS at all (variants farther away are in no
bar). SV hits by type against all tested SVs.""")

co(r"""fig, axs = plt.subplots(1, 2, figsize=(13, 3.6), gridspec_kw=dict(wspace=0.25), sharey=True)
b = np.arange(-3000, 3001, 200)
for ax, k, L in zip(axs, ["deletion", "insertion"], "AB"):
    t = sv[sv.kind == k]; h = t[t.hit]
    tw, hw = t[t.dist_tss.abs() <= 3000], h[h.dist_tss.abs() <= 3000]
    ax.hist(tw.dist_tss, bins=b, weights=np.full(len(tw), 1 / len(t)), histtype="step",
            color=REF, lw=1.3, label=f"all tested ({len(t):,})")
    ax.hist(hw.dist_tss, bins=b, weights=np.full(len(hw), 1 / len(h)), histtype="step",
            color=KIND[k], lw=2.2, label=f"GEA hits ({len(h)})")
    ax.axvline(0, color=TH.MUTED, lw=1)
    ax.set_xlabel("distance to nearest TSS (bp; negative = upstream)")
    TH.grid_only(ax, "y"); TH.panel(ax, L, f"SV {k}s", y=1.12); ax.legend(loc="upper right")
axs[0].set_ylabel("share of all variants in the set")
plt.savefig(f"{PLOTS}/sv_hits_tss_distance.png", dpi=170, bbox_inches="tight"); plt.show()
for k in KIND:
    t = sv[sv.kind == k]; near = lambda d: (d.dist_tss.abs() <= 1000).mean()
    print(f"{k:9s} within 1 kb of a TSS: hits {near(t[t.hit]):.3f}, tested {near(t):.3f}")""")

md(r"""## 6 · The SV hits, by region

Every SV hit with its region and the gene the region call is about (`gene` from the candidate
pool's classification — for a promoter, the gene whose promoter it is; for intergenic / TE /
gene desert, no host gene). Sorted by region, then by p. Also written to
`genes/convergence/results/gea_sv_hits_by_region.csv`.""")

co(r"""V = pd.read_csv(f"{CONV}/results/variants_classified.csv", low_memory=False, usecols=K4 + ["gene", "nearest_gene"])
H = sv[sv.hit].merge(V.drop_duplicates(K4), on=K4, how="left")
H["host_gene"] = np.where(H.region.isin(["CDS", "UTR / ncRNA exon", "intron", "promoter (≤1 kb)"]), H.gene, "")
H["region"] = pd.Categorical(H.region, RORDER, ordered=True)
H = H.sort_values(["region", "min_p"])
cols = ["region", "kind", "size", "chrom", "pos", "host_gene", "nearest_gene", "MAF", "n_hit_axes", "hit_axes", "min_p"]
H[cols].to_csv(f"{CONV}/results/gea_sv_hits_by_region.csv", index=False)
print(H.groupby(["region", "kind"], observed=True).size().unstack(fill_value=0).to_string())
H[cols].assign(MAF=H.MAF.round(3), min_p=H.min_p.map(lambda v: f"{v:.1e}")).reset_index(drop=True)""")

md(r"""## 7 · Warm or cold? The direction of the GEA hits

**The sign.** The production scan (`tile_cr50`) keeps the effect size B and the z-score next to each
p-value, so the direction is the GEA's own, read from the same fit that produced the p-values — not a
new statistic. (The earlier route, a separate signed re-run of the pre-floor scan in
`genes/convergence/lfmm_signed/`, is superseded.)

**What "warm" means.** The input is the change in frequency of the ALT allele and the climate axis is
z-scored, so z > 0 means the ALT allele rose more where the axis is higher. For a deletion the ALT
allele *is* the deletion. Nine axes measure temperature: bio1, bio5, bio6, bio8, bio9, bio10, bio11,
and pc1 and pc3, which point towards higher bio1 and bio10. On all of them z > 0 = favoured **warm**.

**One sign per variant.** The temperature axes are strongly correlated, so a variant is not counted
once per axis: it takes its sign on the temperature axis where it is **most** significant
(`gea_direction_by_record.py`). *Temperature hits* are GEA hits on at least one temperature axis.

**The figures.** One per variant class — SVs first, then small indels — each with deletions and
insertions side by side. Every bar is a **count of hits**, drawn in the direction of its effect:
left (blue) = hits whose allele rises in COLDER gardens (z < 0), right (red) = hits whose allele
rises in WARMER gardens (z > 0). For a deletion the allele is the deletion itself, so a red bar
means "this deletion is favoured in the warm". The row label gives the total and the warm share;
p (binomial, against that region and type's warm share among all tested variants, which sits at
~50%) is in the table below.""")

co(r"""from scipy.stats import binomtest
DR = pd.read_csv(f"{CONV}/results/gea_direction_by_record.csv.gz")
DR["region"] = DR.tier_1kb.map(REGION)
COLD, WARM = "#2166AC", "#B2182B"          # direction is a polarity: the theme's cold/warm pair
CLSLAB = {"sv": "SV", "smallindel": "small-indel"}
rows = []
for cls in ("sv", "smallindel"):
    for k in KIND:
        d = DR[(DR.cls == cls) & (DR.kind == k)]
        for r in RORDER:
            x = d[d.region == r]; h = x[x.temp_hit]
            warm = int(h.warm.sum()); cold = len(h) - warm
            base = x.warm.mean()
            rows.append(dict(cls=cls, kind=k, region=r, hits=len(h), warm=warm, cold=cold,
                             warm_share=warm / len(h) if len(h) else np.nan, tested_warm_share=base,
                             p_vs_tested=binomtest(warm, len(h), base).pvalue if len(h) else np.nan))
W = pd.DataFrame(rows)

def direction_figure(cls):
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.0), gridspec_kw=dict(wspace=0.6))
    y = np.arange(len(RORDER))
    lim = W[W.cls == cls][["warm", "cold"]].max().max() * 1.75 + 1   # headroom for the % labels
    for ax, k, L in zip(axs, KIND, "AB"):
        w = W[(W.cls == cls) & (W.kind == k)].set_index("region").reindex(RORDER)
        ax.barh(y, w.warm, color=WARM, height=0.62, label="warmer gardens (z > 0)")
        ax.barh(y, -w.cold, color=COLD, height=0.62, label="colder gardens (z < 0)")
        ax.axvline(0, color=TH.LABEL, lw=1)
        for i, (n, ws) in enumerate(zip(w.hits, w.warm_share)):
            if n:
                ax.annotate(f"{ws:.0%} warm", xy=(lim * 0.98, i), ha="right", va="center",
                            fontsize=7.5, color=TH.MUTED)
        ax.set_yticks(y); ax.set_yticklabels([f"{r}   n={int(n)}" for r, n in zip(RORDER, w.hits)])
        ax.invert_yaxis(); ax.set_xlim(-lim, lim)
        t = np.array([t for t in ax.get_xticks() if abs(t) <= lim])
        ax.set_xticks(t); ax.set_xticklabels([f"{abs(int(v))}" for v in t])
        ax.set_xlabel(f"{CLSLAB[cls]} {k} temperature hits")
        TH.grid_only(ax, "x"); TH.panel(ax, L, f"{CLSLAB[cls]} {k}s", y=1.13)
    axs[0].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
    plt.savefig(f"{PLOTS}/temperature_hits_direction_{cls}.png", dpi=170, bbox_inches="tight"); plt.show()

th = DR[DR.temp_hit]
print(f"temperature hits: {len(th):,} ({int((th.cls == 'sv').sum())} SV, {int((th.cls == 'smallindel').sum())} small indel); "
      f"leaning warm {th.warm.mean():.3f}  (all tested: {DR.warm.mean():.3f})")
direction_figure("sv")""")

co(r"""direction_figure("smallindel")
print(W.assign(**{c: W[c].round(3) for c in ["warm_share", "tested_warm_share", "p_vs_tested"]}).to_string(index=False))""")

md(r"""**Promoter against intergenic deletions, axis by axis.** The comparison that matters is not
"do hits lean warm" — the table above shows most regions do — but whether **promoter** deletions lean
warmer than **intergenic** ones, the contrast of the burden model. Both are deletions relative to
Col-0, so anything that makes non-reference alleles as a class rise in warm gardens cancels out.
Below: the mean z of promoter deletions minus that of intergenic deletions, over **every tested**
deletion (the whole signed distribution, which is where a small shared shift would show), stratified
by size × MAF, 95% CI from resampling LD blocks (`gea_direction.py`). Above 0 = promoter deletions
lean warmer. Circles: SV deletions; squares: all deletions (small indels dominate).""")

co(r"""GDIR = pd.read_csv(f"{CONV}/results/gea_direction.csv")
TEMPAX = ["bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11", "pc1", "pc3"]
fig, ax = plt.subplots(figsize=(10.5, 3.6)); x = np.arange(len(TEMPAX))
for j, (cls, mk, lab) in enumerate([("sv", "o", "SV deletions"), ("all", "s", "all deletions")]):
    g = GDIR[(GDIR.cls == cls) & (GDIR.kind == "deletion")].set_index("axis").reindex(TEMPAX)
    xx = x + (-0.12 if j == 0 else 0.12)
    ax.errorbar(xx, g.dz_promoter_minus_intergenic, yerr=[g.dz_promoter_minus_intergenic - g.dz_lo, g.dz_hi - g.dz_promoter_minus_intergenic],
                fmt="none", ecolor=KIND["deletion"], lw=1.3)
    ax.scatter(xx, g.dz_promoter_minus_intergenic, marker=mk, s=48, color=KIND["deletion"] if j == 0 else "white",
               edgecolor=KIND["deletion"], lw=1.4, zorder=3, label=lab)
ax.axhline(0, color=TH.MUTED, lw=1)
ax.set_xticks(x); ax.set_xticklabels(TEMPAX); ax.set_ylabel("mean z, promoter − intergenic deletions")
TH.grid_only(ax, "y"); ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
plt.savefig(f"{PLOTS}/promoter_vs_intergenic_direction.png", dpi=170, bbox_inches="tight"); plt.show()
c = ["axis", "cls", "hits_promoter", "warm_share_hits_promoter", "hits_intergenic", "warm_share_hits_intergenic",
     "dz_promoter_minus_intergenic", "dz_lo", "dz_hi", "dz_p"]
print(GDIR[(GDIR.kind == "deletion") & GDIR.axis.isin(TEMPAX) & GDIR.cls.isin(["sv", "all"])][c].round(3).to_string(index=False))""")

md(r"""**Reading it.** ⚠️ *The written reading of this section has NOT been carried over from the
pre-re-run notebook, and deliberately so — it quoted counts and percentages from the record-level,
pre-call-rate scan, which no longer describe these tables. Those numbers are preserved in
`notebooks/old_lfmm/gea_hit_landscape.ipynb` if you want the comparison.*

What to read off the output above, in order:

1. the **baseline** — among all tested variants, is the warm/cold split near 50/50? (it should be;
   a departure means the sign itself is polarised and every share below inherits that);
2. the **hit** warm share against that baseline, per region and per class, with the binomial test
   each cell carries;
3. the **contrast** `promoter − intergenic` on the temperature axes, which is the quantity the
   burden model actually makes a claim about. Both arms are deletions relative to Col-0, so a
   reference-polarisation effect is shared and cancels there.

Nine temperature axes are strongly correlated, so a result on one of them is not nine chances —
read the sign across the row, not the smallest p in it.""")

md(r"""## 8 · The burden model's variant set

The burden model counts SVs as **events** (the alleles of a pangenome bubble collapsed with
`truvari collapse`), uses its own context — CDS > UTR > intron > reference TE > promoter (≤ 2 kb) >
intergenic — and keeps only **shared** variants (carried by both cold- and warm-origin founders),
called in ≥ 90% of founders with ≥ 2 carriers. Those labels exist per event in
`sv_event_catalog.csv`, and since the re-run the GEA counts SVs as events too, so the two are now
on the same unit and the burden's labels apply to the GEA's own signed z unchanged
(`event_gea/event_direction_burden_defs.py`, run on `event_gea/cr50/`). Note the call-rate bars
still differ — the burden set asks for ≥ 90%, the GEA for ≥ 50% — so the filtered arm is a strict
subset. Same contrast and CI as section 7; open = all tested SV deletion events, filled = the
burden's set.""")

co(r"""BD = pd.read_csv(f"{CONV}/event_gea/cr50/event_direction_burden_defs.csv")
fig, ax = plt.subplots(figsize=(10.5, 3.6)); x = np.arange(len(TEMPAX))
for j, (setname, lab, fill) in enumerate([("all tested", "all tested SV deletion events", False),
                                          ("shared + passes (burden filter)", "shared + well-called (burden filter)", True)]):
    g = BD[(BD.kind == "deletion") & (BD.set == setname)].set_index("axis").reindex(TEMPAX)
    xx = x + (-0.12 if j == 0 else 0.12)
    ax.errorbar(xx, g.dz, yerr=[g.dz - g.lo, g.hi - g.dz], fmt="none", ecolor=KIND["deletion"], lw=1.3)
    ax.scatter(xx, g.dz, s=48, color=KIND["deletion"] if fill else "white", edgecolor=KIND["deletion"], lw=1.4, zorder=3, label=lab)
ax.axhline(0, color=TH.MUTED, lw=1)
ax.set_xticks(x); ax.set_xticklabels(TEMPAX); ax.set_ylabel("mean z, promoter − intergenic deletions")
TH.grid_only(ax, "y"); ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
plt.savefig(f"{PLOTS}/promoter_direction_burden_filter.png", dpi=170, bbox_inches="tight"); plt.show()
print(BD.round(3).to_string(index=False))""")

md(r"""**Reading it.** ⚠️ *As in section 7, the previous written reading is not carried over — it
described the pre-re-run run. Compare `notebooks/old_lfmm/gea_hit_landscape.ipynb` for it.*

The question this section asks: does the GEA's own signed z reproduce the burden model's
warm-favoured promoter deletions **once the burden's own variant definitions are applied**? Read
the open markers (all tested SV deletion events) against the filled ones (shared + well-called)
axis by axis, and ask whether the filled series sits consistently above zero — consistency across
the nine correlated temperature axes is the signal here, not any single p. The filter removes
variants carried only by cold- or only by warm-origin founders, whose frequencies track the fate
of their lineage rather than the variant; if the contrast moves when they go, that is where it
came from. Both panels are per-variant tests, and a burden test sums over all variants, so a
shift too small to clear per-variant significance is not evidence against the burden result.""")

md(r"""## 9 · What the record-level scan would have given

*The comparison arm, kept at `results/multiaxis/lfmm_cr50` and
`genes/convergence/event_gea/cr50/`.*

Everything above counts SVs as **truvari events**, which is the production unit since the
2026-09-30 re-run and the unit the burden and purging work already used. The arm compared here
is the same LFMM (ridge K = 16, raw p, the same Δp input and climate scaling, the same MAF 0.05
and call rate 0.50) run instead on arch3 **records** — one row per assembly path, so one
insertion or deletion recurs as several "alleles" differing by a base or two.

Two things follow from splitting a variant across alleles, and both are visible below. An SV
whose alleles are each below MAF 0.05 is never tested at all, though together they clear it —
those events are **testable only after collapsing**. And an SV whose alleles do clear it is
tested several times over, so its evidence is divided among near-duplicate tests.

Collapsing is not free in the other direction either: a record-arm hit can fall below the bar
once its alleles are summed, because the summed frequency is a different (and usually less
extreme) trajectory than the single allele that happened to carry the signal. The counts for
both directions are printed below — read them together, not just the gains.""")

co(r"""EV = pd.read_csv(f"{CONV}/event_gea/cr50/event_sv_hits.csv")
EV["region"] = EV.tier.map(REGION)
EV["stratum"] = (pd.cut(EV["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                 + pd.cut(EV.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
new = ~EV.rec_tested
print(f"SV events tested by the production (event) GEA: {len(EV):,}")
print(f"  testable ONLY as events (no member record passes MAF 0.05 + call rate 0.50 alone): {int(new.sum()):,}")
print(f"  records the record arm spends per event: mean {EV.n_rec.mean():.2f}, max {int(EV.n_rec.max())}")
print(f"SV event hits: {int(EV.hit.sum())}  ({int((EV.hit & (EV.kind == 'deletion')).sum())} deletions, "
      f"{int((EV.hit & (EV.kind == 'insertion')).sum())} insertions)")
print(f"  from events only testable after collapsing: {int((EV.hit & new).sum())}  "
      f"(hit rate {100 * (EV.hit & new).sum() / new.sum():.2f}% "
      f"vs {100 * (EV.hit & ~new).sum() / (~new).sum():.2f}% among record-testable events)")
print("\ngained and lost by collapsing:")
print(f"  events the record arm also called: {int(EV.rec_hit.sum())}")
print(f"    still a hit as an event:        {int((EV.hit & EV.rec_hit).sum())}")
print(f"    record-only (lost by collapsing): {int((~EV.hit & EV.rec_hit).sum())}")
print(f"  event-only (not called by any member record): {int((EV.hit & ~EV.rec_hit).sum())}")
S_ev = region_stats(EV)
region_figure(EV, S_ev, "SV event", "sv_event_hits_by_region.png")
print(S_ev.assign(**{c: S_ev[c].round(3) for c in ["hit_share", "tested_share", "OR", "lo", "hi", "p"]}).to_string(index=False))""")

co(r"""ED = pd.read_csv(f"{CONV}/event_gea/cr50/event_direction.csv")
g = ED.set_index("axis").reindex(TEMPAX)
fig, ax = plt.subplots(figsize=(10.5, 3.3)); x = np.arange(len(TEMPAX))
ax.errorbar(x, g.dz, yerr=[g.dz - g.dz_lo, g.dz_hi - g.dz], fmt="o", ms=7, color=KIND["deletion"], lw=1.3)
ax.axhline(0, color=TH.MUTED, lw=1)
ax.set_xticks(x); ax.set_xticklabels(TEMPAX); ax.set_ylabel("mean z, promoter − intergenic\nSV deletion events")
TH.grid_only(ax, "y")
plt.savefig(f"{PLOTS}/event_promoter_direction.png", dpi=170, bbox_inches="tight"); plt.show()
print(ED[ED.temperature].round(3).to_string(index=False))""")

md(r"""**Reading it.** Counting SVs as events recovers **5,345 SV events** the record-level GEA never
tested, and they are not noise: 88 of them are hits, a hit rate (1.65%) a little above that of events
the GEA could already test (1.23%). The unit also **reshuffles which SVs are hits**: of the 344 events
that carry a record-level SV hit, only 171 are still hits when their alleles are summed, and 184 event
hits are new. Where the event hits fall looks like section 2 with **one exception: coding (CDS)
deletions are depleted among the event hits** — 11 against about 23 expected (odds ratio 0.40, 95% CI
0.19–0.67), which survives correction across the 14 region × type cells. At the record level the same
direction was there but not significant (0.68), so counting SVs as events sharpens it. It fits coding
deletions rarely being the variants that track climate, though this analysis cannot say why (they
are adjusted for size and frequency here). The promoter-vs-intergenic direction at the event level
matches the record-level one (slightly cold on most temperature axes, none significant). The unit, then, is not what separates the GEA from
the burden model's direction; the variant set is (section 8). Because the event unit changes which
individual SVs are hits, it bears on the candidate list, which other lines of work are handling — it
is kept here as a parallel check.""")

nb = new_notebook(); nb["cells"] = C
# history off: the shared ~/.ipython history.sqlite is often locked by another open kernel, which
# made a build die at startup and made IPython move the user's history file aside
ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180,
                         extra_arguments=["--HistoryManager.enabled=False"])
ep.preprocess(nb, {"metadata": {"path": NBDIR}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
