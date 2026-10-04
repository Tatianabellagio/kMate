#!/usr/bin/env python
"""Build+execute `gea_te_cargo.ipynb`: are the climate-associated SVs transposons?

Joins the GEA hit list to the RepeatMasker cargo annotation and asks three things: whether significant
hits are enriched for TE sequence, which families carry them, and whether those families relate to the
carriers' origin climate. Runs the whole thing twice -- at the production 50% call-rate floor and at the
archived run with no call-rate filter -- because SVs seen only in the long-read panel cannot pass a 50%
floor by construction.

`basic` env.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import build, new_markdown_cell, new_code_cell, code_setup_burden as code_setup

md_title = """# Are the climate-associated SVs transposons?

The GEA finds structural variants whose frequency tracks climate. This notebook asks what those variants
*are*, by joining them to the repeat classification of every SV record, insertions and deletions
(`r1_sv_negative_selection/TE_CLASSIFICATION_PIPELINE.md`: GraffiTE's masking protocol, applied to all
227,851 SV records).

## The hit list

**`v2_events_cr50`**, built 2026-10-02 — the current production GEA:

- **scans** LFMM ridge K = 16, generation 9, 352 pools, 22 axes (bio1–19, pc1–3), **raw p**, no GIF
  correction (lambda recorded, never applied)
- **unit** SV = a Truvari **event** (collapse 0.95/0.95/500 bp within a bubble and type)
- **filters** MAF ≥ 0.05 at gen-9 pool-mean frequency, **and called in ≥ 116 of 231 founders (50%)**
- **hit rule** Bonferroni 0.05 / n units per class × axis, union over scans
- **result** 20,134 SV events tested → **282 hits** (151 insertions, 131 deletions)

**The second run** is the archived `default_no_callrate` version: identical except with **no call-rate
floor** — 27,009 SV events tested → **355 hits**. It is here because the call-rate filter is not neutral
with respect to the question: the long-read panel is 78 of 231 founders (**34%**), so a variant called
only there cannot clear a 50% floor. If TE-rich variants are preferentially long-read-only, the
production run would be blind to exactly the class this notebook is about.

## Two traps this has to handle

**The join.** The GEA pool and our catalogue share `chrom` and `pos`, but 2.14% of panel positions carry
more than one record, so a position join silently matches the wrong allele. Every join here is on
**chrom + position + direction + size (within 5 bp)**, taking the best size match per unit, and units that
do not match uniquely are dropped rather than guessed at.

**The background.** "Are hits TE-rich?" needs the comparison to be against the variants that were
**tested**, not against every SV event in the panel. The tested universe is filtered on MAF and call
rate, and those filters change its TE content: the tested insertions are markedly TE-poorer than the panel
as a whole. Compared with the whole panel instead, the hits look ordinary and the enrichment disappears."""

code_load = r"""
import os
E=os.path.normpath(f"{G}/../../../genes/convergence/event_gea")
assert os.path.isdir(E), E
CAT=pd.read_csv(f"{G}/sv_event_catalog.csv")
# repeat class of every SV record, insertions AND deletions (TE_CLASSIFICATION_PIPELINE.md); an event takes
# the class of its representative record, which is the catalogue's `key`
RM=pd.read_csv(f"{G}/graffite/sv_repeat_class_records.csv").rename(columns={
    "span_te":"te_frac","top_family":"family","top_superfamily":"cls","top_div":"div"})
RM["repl"]=(RM.hprc_type=="COMPLEX")|(RM.unshared_bp>0)
assert CAT.key.isin(RM.key).all()
CTXC={"CDS":"#2F4A1F","UTR":"#4E6B3A","intron":"#7E9B6B","TE":"#C46A10",
      "promoter":"#5B7A99","intergenic":"#B9B9B9"}
SZ=[50,100,300,1000,3000,10000,10**9]
def load_run(path):
    d=pd.read_csv(path)
    d["kd"]=np.where(d.kind=="insertion","ins","del")
    j=d.merge(CAT[["chrom","pos","size_bp","kind","key","context","n_called"]]
              .rename(columns={"kind":"kd","n_called":"cat_called"}),on=["chrom","pos","kd"],how="left")
    j["dsz"]=(j.size_bp-j["size"]).abs()
    j=j.sort_values("dsz").drop_duplicates("unit",keep="first")
    j=j[(j.dsz<=5)&j.key.notna()]
    return j.merge(RM[["key","te_frac","cargo","family","cls","div","repl","dominant"]],on="key",how="left")
RUNS={}
RUNS["cr50"]=load_run(f"{E}/cr50/event_sv_hits.csv")
RUNS["no call-rate"]=load_run(f"{E}/archive_default_no_callrate_2026-09-22/event_sv_hits.csv")
for k,v in RUNS.items():
    raw=pd.read_csv(f"{E}/cr50/event_sv_hits.csv" if k=="cr50"
                    else f"{E}/archive_default_no_callrate_2026-09-22/event_sv_hits.csv")
    print(f"{k:14s} tested {len(raw):,} SV events, {int(raw.hit.sum()):,} hits | "
          f"uniquely joined {len(v):,} ({100*len(v)/len(raw):.0f}%), of which hits {int(v.hit.sum())}")
TEHALF=["te_full","te_most"]
"""

md_enrich = """## 1 · Are the hits enriched for TE sequence?

**Insertions and deletions.** A deletion's cargo is the reference sequence it removes, masked with the same
protocol as an insertion's ALT, so both directions now have a class. Hits are compared against the
**tested background of the same run**, which is the comparison that controls for MAF and call rate.

Two tests per panel:
- **Fisher**, hits vs tested: the raw comparison;
- **Mantel–Haenszel**, the same comparison **within strata of SV size × replacement status**. TE content
  rises steeply with size, and replacement records (HPRC `COMPLEX`, or a swapped base) are much less often
  whole TEs at the same size. If hits differ from the background in either, the raw odds ratio is partly
  that. Strata without a hit carry no information and are dropped."""

code_enrich = r"""
from statsmodels.stats.contingency_tables import StratifiedTable
fig,axs=plt.subplots(2,2,figsize=(11.8,6.8),gridspec_kw=dict(wspace=0.40,hspace=0.60))
SUM=[]
for r,kd in enumerate(["ins","del"]):
    for c,(tag,J) in enumerate(RUNS.items()):
        a=axs[r,c]; style(a,axis="y")
        d=J[(J.kd==kd)&J.te_frac.notna()].copy()
        d["te"]=d.cargo.isin(TEHALF)
        H=d[d.hit]; B=d[~d.hit]
        ct=np.array([[int(H.te.sum()),int((~H.te).sum())],[int(B.te.sum()),int((~B.te).sum())]])
        orr,p=stats.fisher_exact(ct)
        d["stratum"]=pd.cut(d.size_bp,SZ).astype(str)+"|"+d.repl.astype(str)
        tabs=[]
        for s,g in d.groupby("stratum"):
            t=np.array([[(g.hit&g.te).sum(),(g.hit&~g.te).sum()],
                        [(~g.hit&g.te).sum(),(~g.hit&~g.te).sum()]],dtype=float)
            if t[0].sum()>0 and t[1].sum()>0 and t[:,0].sum()>0 and t[:,1].sum()>0: tabs.append(t)
        st=StratifiedTable(tabs); omh=st.oddsratio_pooled; pmh=st.test_null_odds().pvalue
        a.bar([0,1],[100*B.te.mean(),100*H.te.mean()],color=["0.72",PAL[f"{kd}_promoter"]],width=0.6)
        a.set_xticks([0,1]); a.set_xticklabels([f"tested\n(n={len(B):,})",f"GEA hits\n(n={len(H)})"],fontsize=8)
        a.set_ylabel(f"% of {'insertions' if kd=='ins' else 'deletions'} at least half TE",fontsize=8.5)
        a.set_ylim(0,100*max(B.te.mean(),H.te.mean())*1.55)   # headroom for the annotation
        a.annotate(f"{tag}, {'insertions' if kd=='ins' else 'deletions'}\n"
                   f"Fisher OR {orr:.2f}, p = {p:.2g}\nsize x replacement MH OR {omh:.2f}, p = {pmh:.2g}",
                   xy=(0.03,0.97),xycoords="axes fraction",va="top",fontsize=7.4,fontweight="bold",color="0.25")
        SUM.append(dict(kind=kd,run=tag,n_hit=len(H),n_bg=len(B),pct_hit=100*H.te.mean(),pct_bg=100*B.te.mean(),
                        OR=orr,p=p,OR_MH=omh,p_MH=pmh,n_strata=len(tabs),
                        med_size_hit=H.size_bp.median(),med_size_bg=B.size_bp.median(),
                        pct_repl_hit=100*H.repl.mean(),pct_repl_bg=100*B.repl.mean()))
save(fig,"geate_enrichment",nb=True); plt.show()
print("TE CONTENT OF GEA HITS vs THE TESTED BACKGROUND  (te_full + te_most = at least half TE)\n")
print(pd.DataFrame(SUM).round(4).to_string(index=False))
ce=CAT.merge(RM[["key","cargo"]],on="key")
print("\nfor reference, ALL panel SV events (no MAF / call-rate filter): "
      + ", ".join(f"{k} {100*ce[ce.kind==k].cargo.isin(TEHALF).mean():.1f}% at least half TE"
                  for k in ("ins","del")))
"""

md_fam = """## 2 · Which families, and where the hits land

Two different statements that are easy to conflate: the **cargo** of a hit (what sequence the insertion
carries) and its **landing site** (what it inserted into, by the burden model's cascade). A variant can
carry a TE and land in a promoter, or carry nothing recognisable and land inside an annotated TE."""

code_fam = r"""
fig,axs=plt.subplots(1,2,figsize=(12.6,3.8),gridspec_kw=dict(wspace=0.46))
J=RUNS["cr50"]; H=J[J.hit]
a=axs[0]; style(a,axis="x")
FT=pd.crosstab(H[H.cargo.isin(TEHALF)].family,H[H.cargo.isin(TEHALF)].kd).reindex(columns=["ins","del"]).fillna(0)
FT=FT.loc[FT.sum(1).sort_values(ascending=False).index[:12]]
yy=np.arange(len(FT))[::-1]
a.barh(yy,FT["ins"].to_numpy(),color=PAL["ins_promoter"],height=0.64,label="insertion hits")
a.barh(yy,FT["del"].to_numpy(),left=FT["ins"].to_numpy(),color=PAL["del_promoter"],height=0.64,label="deletion hits")
a.set_yticks(yy); a.set_yticklabels(FT.index,fontsize=7.4); a.legend(fontsize=7,frameon=False,loc="lower right")
a.set_xlabel("GEA hits (at least half TE) where this family dominates",fontsize=8.5)
a.annotate("a   cargo family of the hits",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
a=axs[1]; style(a,axis="y")
CX=["CDS","UTR","intron","TE","promoter","intergenic"]
hc=H.context.value_counts().reindex(CX).fillna(0)
bc=J[~J.hit].context.value_counts().reindex(CX).fillna(0)
x=np.arange(len(CX))
a.bar(x-0.19,100*bc/bc.sum(),width=0.36,color="0.72",label="tested")
a.bar(x+0.19,100*hc/hc.sum(),width=0.36,color=[CTXC[c] for c in CX],label="hits")
a.set_xticks(x); a.set_xticklabels(CX,fontsize=7.5,rotation=25,ha="right")
a.set_ylabel("% of SV events",fontsize=8.5); a.legend(fontsize=7,frameon=False)
a.annotate("b   landing site (burden cascade)",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
save(fig,"geate_families",nb=True); plt.show()
print("CARGO FAMILIES OF THE GEA HITS (cr50; hits at least half TE)\n"); print(FT.astype(int).to_string())
for kd in ("ins","del"):
    h=H[(H.kd==kd)&H.cargo.isin(TEHALF)]
    print(f"\n  {kd}: {len(h)} TE hits, {h.family.nunique()} families; most frequent family in "
          f"{h.family.value_counts().iloc[0]} of them; superfamilies {h.cls.value_counts().head(4).to_dict()}")
print("\nLANDING SITE, hits vs tested (cr50)\n")
t=pd.DataFrame({"hits":hc.astype(int),"% hits":(100*hc/hc.sum()).round(1),
                "tested":bc.astype(int),"% tested":(100*bc/bc.sum()).round(1)})
print(t.to_string())
"""

md_cr = """## 3 · What the call-rate filter removes

The production GEA requires a variant to be called in ≥ 116 of 231 founders. The long-read panel is **78
founders (34%)**, so an SV detectable only from assemblies cannot clear that floor however real it is.
Since TE insertions are exactly the class most likely to be assembly-only, the filter could in principle
hide the signal this notebook is looking for.

Three things to check: how many events the floor removes, whether those events are cactus-panel-sized,
and whether they are more TE-rich than the ones that survive."""

code_cr = r"""
cc=CAT.merge(RM[["key","te_frac","cargo"]],on="key",how="left")
ins=cc[(cc.kind=="ins")&cc.te_frac.notna()]
lo=ins[ins.n_called<116]; hi=ins[ins.n_called>=116]
fig,axs=plt.subplots(1,3,figsize=(13.4,3.6),gridspec_kw=dict(wspace=0.42))
a=axs[0]; style(a,axis="both")
a.hist(ins.n_called,bins=60,color="0.72")
a.axvline(116,color="#B04A00",lw=1.6,ls=(0,(4,2)))
a.axvline(78,color="#1F4E79",lw=1.2,ls=(0,(2,2)))
a.annotate("50% floor",xy=(118,a.get_ylim()[1]*0.92),fontsize=7,color="#B04A00")
a.annotate("cactus panel = 78",xy=(80,a.get_ylim()[1]*0.75),fontsize=7,color="#1F4E79")
a.set_xlabel("founders the event is called in",fontsize=8.5); a.set_ylabel("insertion events",fontsize=8.5)
a.annotate("a   where the floor sits",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
a=axs[1]; style(a,axis="y")
a.bar([0,1],[100*lo.cargo.isin(TEHALF).mean(),100*hi.cargo.isin(TEHALF).mean()],
      color=["#1F4E79","0.72"],width=0.6)
a.set_xticks([0,1]); a.set_xticklabels([f"removed\n(n={len(lo):,})",f"kept\n(n={len(hi):,})"],fontsize=8)
a.set_ylabel("% at least half TE",fontsize=8.5)
a.annotate("b   TE content, removed vs kept",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
a=axs[2]; style(a,axis="y")
S=[]
for tag,J in RUNS.items():
    i2=J[(J.kd=="ins")&J.te_frac.notna()]
    S.append((tag,100*i2[i2.hit].cargo.isin(TEHALF).mean(),int(i2.hit.sum())))
a.bar(range(len(S)),[s[1] for s in S],color=["#B04A00","#1F4E79"],width=0.55)
a.set_xticks(range(len(S))); a.set_xticklabels([f"{s[0]}\n({s[2]} hits)" for s in S],fontsize=8)
a.set_ylabel("% of hits at least half TE",fontsize=8.5)
a.annotate("c   does relaxing it change the answer?",xy=(0,1.06),xycoords="axes fraction",
           fontsize=8.4,fontweight="bold",color="0.3")
save(fig,"geate_callrate",nb=True); plt.show()
print("WHAT THE 50% CALL-RATE FLOOR REMOVES\n")
print(f"  insertion events below the floor: {len(lo):,} of {len(ins):,} ({100*len(lo)/len(ins):.1f}%)")
print(f"  of those, called in <= 85 founders (cactus-panel sized): "
      f"{int((lo.n_called<=85).sum()):,} ({100*(lo.n_called<=85).mean():.0f}% of the removed)")
print(f"\n  >=half TE:  removed {100*lo.cargo.isin(TEHALF).mean():.1f}%   "
      f"kept {100*hi.cargo.isin(TEHALF).mean():.1f}%")
print(f"  mean TE fraction: removed {lo.te_frac.mean():.3f}   kept {hi.te_frac.mean():.3f}")
print("\n  hit sets, % at least half TE: " + ", ".join(f"{s[0]} {s[1]:.1f}% of {s[2]}" for s in S))
"""

md_famtest = """## 4 · Is any individual family enriched — and where is ONSEN?

The overall excess is modest and could be one family or many. This tests every family with at least 20
tested insertions, Fisher against all other tested insertions, Bonferroni over the families tested.

It also looks at a question the previous sections raise. ATCOPIA78 / ONSEN is the one family with a clear
climate signal in the **carrier** analysis (`sv_insertion_cargo_repeatmasker.ipynb` §6), so why is it
absent from the hit list? This is counted per insertion **event**, the project's SV unit, so one insertion
seen in several panel alleles counts once."""

code_famtest = r'''J=RUNS["cr50"]; ins=J[(J.kd=="ins")&J.te_frac.notna()]
H=ins[ins.hit]; Bg=ins[~ins.hit]
rows=[]
for fam,g in ins.groupby("family"):
    if not isinstance(fam,str) or fam=="" or len(g)<20: continue
    a1=int(g.hit.sum()); b1=len(g)-a1
    orr,p=stats.fisher_exact([[a1,b1],[len(H)-a1,len(Bg)-b1]])
    rows.append(dict(family=fam,n_tested=len(g),n_hits=a1,hit_pct=100*a1/len(g),OR=orr,p=p))
T=pd.DataFrame(rows).sort_values("p")
T["p_bonf"]=np.minimum(1,T.p*len(T))
fig,axs=plt.subplots(1,2,figsize=(12.4,3.7),gridspec_kw=dict(wspace=0.46))
a=axs[0]; style(a,axis="both")
a.axhline(100*len(H)/len(ins),color="0.45",lw=1.1,ls=(0,(4,2)))
a.annotate("overall hit rate",xy=(0.98,100*len(H)/len(ins)),xycoords=("axes fraction","data"),
           ha="right",va="bottom",fontsize=6.8,color="0.45")
sig=T.p<0.05
a.scatter(T.loc[~sig,"n_tested"],T.loc[~sig,"hit_pct"],s=26,color="0.72",linewidth=0.3,edgecolor="0.4")
a.scatter(T.loc[sig,"n_tested"],T.loc[sig,"hit_pct"],s=42,color="#B04A00",linewidth=0.4,edgecolor="0.3")
for r in T[sig].itertuples():
    a.annotate(r.family[:18],xy=(r.n_tested,r.hit_pct),fontsize=6.4,color="#B04A00",
               xytext=(3,3),textcoords="offset points")
a.set_xscale("log"); a.set_xlabel("insertions of this family tested (log)",fontsize=8.5)
a.set_ylabel("% that are GEA hits",fontsize=8.5)
a.annotate("a   per-family hit rate",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
a=axs[1]; style(a,axis="both")
# per EVENT (the project's SV unit), founder frequency = carriers / founders called. The GEA's MAF floor is
# on GrENE-Net pool frequency, so the founder frequency is a proxy for it, not the same quantity
CE=CAT[CAT.kind=="ins"].merge(RM[["key","family","te_frac"]],on="key",how="left")
CE["freq"]=CE.n_carriers/CE.n_called
is78=CE.family.astype(str).str.contains("ATCOPIA78",na=False)
A78=CE[is78]; other=CE[(CE.te_frac>0)&~is78]
bins=np.logspace(np.log10(1/231),0,40)
a.hist(other.freq.dropna(),bins=bins,density=True,color="0.78",label="all other TE")
a.hist(A78.freq.dropna(),bins=bins,density=True,color="#B04A00",alpha=0.8,label="ATCOPIA78")
a.axvline(0.05,color="#1F4E79",lw=1.6,ls=(0,(4,2)))
a.annotate("GEA MAF floor",xy=(0.052,a.get_ylim()[1]*0.8),fontsize=7,color="#1F4E79")
a.set_xscale("log"); a.set_xlabel("founder frequency (insertion events)",fontsize=8.5)
a.set_ylabel("density",fontsize=8.5); a.legend(fontsize=7,frameon=False)
a.annotate("b   why ONSEN cannot be a hit",xy=(0,1.06),xycoords="axes fraction",fontsize=8.4,
           fontweight="bold",color="0.3")
save(fig,"geate_family_test",nb=True); plt.show()
print(f"PER-FAMILY TEST  ({len(T)} families with >=20 tested insertions; "
      f"overall hit rate {100*len(H)/len(ins):.2f}%)\n")
print(T.head(10).round(4).to_string(index=False))
print(f"\n  families passing Bonferroni: {int((T.p_bonf<0.05).sum())}")
print(f"  superfamilies of the nominal (p < 0.05) families: "
      f"{ins[ins.family.isin(T[T.p<0.05].family)].drop_duplicates('family').set_index('family').cls.to_dict()}")
print("\nWHY ATCOPIA78 / ONSEN IS ABSENT FROM THE HIT LIST  (per insertion EVENT)\n")
print(f"  ATCOPIA78 insertion events in the panel: {len(A78):,}")
print(f"    median carriers {A78.n_carriers.median():.0f}, median founder frequency {A78.freq.median():.4f}")
print(f"    founder frequency >= 0.05: {int((A78.freq>=0.05).sum()):,} ({100*(A78.freq>=0.05).mean():.1f}%)")
print(f"  all other TE insertion events with founder frequency >= 0.05: {100*(other.freq>=0.05).mean():.1f}%")
n78={k:int((J_.kd.eq("ins")&J_.family.astype(str).str.contains("ATCOPIA78")).sum()) for k,J_ in RUNS.items()}
print(f"  ATCOPIA78 events tested by the GEA: " + ", ".join(f"{k} {v}" for k,v in n78.items()))
'''

md_maf = """## 5 · Would dropping the MAF filter help?

GEA scans normally impose a MAF floor because rare-variant association statistics are unstable. That is
the textbook reason, and it is worth checking whether it actually applies here before accepting it — the
previous section showed the floor is what hides ONSEN.

The scan already recorded **inflation by MAF bin**, so the question is answerable from what exists."""

code_maf = r'''MA=pd.read_csv(f"{os.path.normpath(f'{G}/../../../r2_gea_nonsnp/phase1_replication/results/multiaxis')}"
               "/axis_inflation_maf_nonsnp.csv")
FL=pd.read_csv(f"{os.path.normpath(f'{G}/../../../r2_gea_nonsnp/phase1_replication/results/multiaxis')}"
               "/axis_inflation_floor_sv.csv")
fig,axs=plt.subplots(1,2,figsize=(12.0,3.6),gridspec_kw=dict(wspace=0.42))
a=axs[0]; style(a,axis="both")
AX=["bio1","bio10","bio11","bio15"]
for ax_ in AX:
    d=MA[MA.axis==ax_].sort_values("maf_lo")
    if len(d): a.plot((d.maf_lo+d.maf_hi)/2,d.lam,marker="o",ms=4,lw=1.5,label=ax_)
a.axhline(1,color="0.5",lw=1.0,ls=(0,(4,2)))
a.set_xlabel("MAF bin (midpoint)",fontsize=8.5); a.set_ylabel("genomic inflation lambda",fontsize=8.5)
a.legend(fontsize=7,frameon=False)
a.annotate("a   inflation does NOT rise at low MAF",xy=(0,1.06),xycoords="axes fraction",
           fontsize=8.4,fontweight="bold",color="0.3")
a=axs[1]; style(a,axis="both")
for ax_ in AX:
    d=FL[FL.axis==ax_].sort_values("floor")
    if len(d): a.plot(d.floor,d.lam,marker="s",ms=4,lw=1.5,label=ax_)
a.axhline(1,color="0.5",lw=1.0,ls=(0,(4,2)))
a.set_xlabel("MAF floor applied",fontsize=8.5); a.set_ylabel("lambda among SVs kept",fontsize=8.5)
a.legend(fontsize=7,frameon=False)
a.annotate("b   raising the floor does not fix it either",xy=(0,1.06),xycoords="axes fraction",
           fontsize=8.4,fontweight="bold",color="0.3")
save(fig,"geate_maf_inflation",nb=True); plt.show()
print("INFLATION BY MAF BIN (non-SNP units)\n")
print(MA[MA.axis.isin(AX)].pivot_table(index=["maf_lo","maf_hi"],columns="axis",
                                       values="lam").round(3).to_string())
print("\nINFLATION AS THE MAF FLOOR IS RAISED (SV units)\n")
print(FL[FL.axis.isin(AX)].pivot_table(index="floor",columns="axis",values="lam").round(3).to_string())
print("\n  Lambda is FLAT across MAF, and if anything slightly LOWER in the rarest bin tested")
print("  (bio1: 2.27 at MAF 0.05-0.075 against 2.59 at 0.30-0.50). Raising the floor to 0.25")
print("  discards 80% of SV units and leaves lambda unchanged at ~2.5. So in this dataset the")
print("  inflation is climate structure, not rare-variant noise, and the MAF filter is not")
print("  buying calibration -- which is the usual justification for having one.")
'''

md_read = """## 6 · Reading it

**Insertion hits are TE-enriched, and it is not a size effect.** 54.2% of GEA insertion hits are at least
half TE, against 41.5% of the insertions tested alongside them (Fisher OR 1.67, p = 0.004). Hits and
background have nearly the same median size (256 vs 262 bp) and replacement share (21% vs 22%). Within
size × replacement strata the odds ratio is 1.94 (Mantel–Haenszel p = 0.001). The no-call-rate run agrees
(55.3% vs 42.8%, OR 1.65; stratified 1.61, p = 0.007).

**Deletion hits are not.** 40.0% vs 36.9% (OR 1.14, p = 0.51; stratified 1.14, p = 0.51), with a weak,
non-significant tendency in the no-call-rate run (OR 1.35, p = 0.08). Whatever ties TE cargo to the
climate signal is an insertion-side pattern. "Insertion" and "deletion" are relative to Col-0, not ancestral
state, and without TSDs that can't be resolved here.

⚠ **The background still matters.** The hits (54.2%) look just like the panel's insertion events as a whole
(52.8%). The enrichment exists only against the tested universe, which the MAF and call-rate filters make
TE-poorer. That is the right comparator, because it is what the GEA could have found, but the enrichment
is a statement about the tested set, not the panel.

**No single family drives it.** The 71 insertion hits that are at least half TE come from 58 families;
the most frequent appears 5 times. No family passes Bonferroni. The nominal top three are ATREP3
(RC/Helitron), L1-84_AT (LINE/L1) and ATREP18 (labelled plain `DNA` in the TAIR10 part of the library,
although most ATREP families there are Helitrons). RC/Helitron is the commonest superfamily among both
insertion and deletion hits.

**The call-rate floor is not hiding TE-rich variants.** It removes 10,460 insertion events, 98% of them
called in ≤ 85 founders, i.e. assembly-only as expected. These are slightly *less* TE-rich than the
events that survive (49.9% vs 53.6%), the opposite of the worry that motivated the second run. The
earlier `-nolow` classification showed them as slightly *more* TE-rich (60.7% vs 56.9%). The sign flips
because assembly-only events lose more whole-TE calls under the new masking than the rest do, and that
happens among their clean records, not only among replacements. A gap this small rests on the masking
protocol, so read it as "no TE excess among the removed events", not as a depletion.

**ONSEN is rare, not invisible.** Counted per event, there are 298 ATCOPIA78 insertion events. Their median
founder frequency is 0.010, and 7.7% reach 0.05, against 24.3% of other TE insertion events. Eight
were tested by the cr50 GEA (14 without the call-rate floor) and none is a hit. ONSEN does appear on the
other side: three cr50 **deletion** hits remove an ATCOPIA78 copy that Col-0 carries.

**What is still missing**
- **Direction.** This asks what the hits *are*, not which way their frequency moves with climate. Pairing
  cargo with the sign of the climate association is the obvious next step.
- **Polarity.** Without target site duplications, an insertion relative to Col-0 cannot be told apart
  from a deletion in the Col-0 lineage (GraffiTE's TSD search is broken upstream; see the pipeline doc).
- **The hit list is a candidate net.** Raw uncalibrated LFMM p, climate-structured hits, no GIF
  correction, as the GEA manifest states. Nothing here is a calibrated enrichment test of that list."""

md_methods = """## As run

- **GEA**: `genes/convergence/event_gea/cr50/event_sv_hits.csv` (production, 50% call-rate floor) and
  `archive_default_no_callrate_2026-09-22/event_sv_hits.csv` (no floor). Both carry a `hit` column from
  the Bonferroni rule; `p_bar` in those files is **allele frequency, not a p-value**.
- **Cargo**: `graffite/sv_repeat_class_records.csv`, built by the pipeline in
  `r1_sv_negative_selection/TE_CLASSIFICATION_PIPELINE.md`. RepeatMasker `-s` against TAIR10 + RepBase
  31.09, low-complexity masking on (GraffiTE's protocol), on the ALT of insertions and the REF of
  deletions. TE fraction = union of TE-class hits / allele length; simple repeat, low complexity and
  satellite do not count as TE. An event takes the class of its representative record (the catalogue
  `key`). Until 2026-10-04 this read `rm_cargo_allins.csv` (`-nolow`, insertions only).
- **Size × replacement strata**: SV size bins 50–100–300–1k–3k–10k–∞ bp × replacement status (HPRC
  `COMPLEX`, or a swapped base); Mantel–Haenszel via `statsmodels` `StratifiedTable`, using strata that
  contain at least one hit.
- **Join**: chrom + position + direction + size within 5 bp, best size match per unit, non-unique matches
  dropped. Never on position alone — 2.14% of panel positions carry more than one record.
- **Landing site** is the burden model's cascade (`context` in `sv_event_catalog.csv`), not the older
  `genic` scheme."""

if __name__ == "__main__":
    build("gea_te_cargo.ipynb", [
        new_markdown_cell(md_title),  new_code_cell(code_setup), new_code_cell(code_load),
        new_markdown_cell(md_enrich), new_code_cell(code_enrich),
        new_markdown_cell(md_fam),    new_code_cell(code_fam),
        new_markdown_cell(md_cr),     new_code_cell(code_cr),
        new_markdown_cell(md_famtest), new_code_cell(code_famtest),
        new_markdown_cell(md_maf),     new_code_cell(code_maf),
        new_markdown_cell(md_read),   new_markdown_cell(md_methods),
    ])
