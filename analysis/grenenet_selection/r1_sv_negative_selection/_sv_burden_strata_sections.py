md_freq = """## 7 · Frequency strata — is the effect carried by rare variants?

Purifying selection predicts that variants are rare *because* they are deleterious, so a coding-SV burden
split by allele frequency should be more damaging in the rare stratum. Frequency is ALT-carrier count over
the full 231-founder panel, a property of the variant.

Bins are **rare** (MAC 2–12, ≤ 5.2 % of the panel) and **common** (MAC 13+), fitted on all 231 founders.
The bottom row, in orange, is the **unstratified** coding-SV coefficient with its own 95 % CI — where the
effect falls when the two strata are held together — and the shaded band carries that interval across the
panel, so each stratum can be read against the undivided estimate rather than only against its point value.

⚠ **Private variants (MAC 1) are excluded.** A variant carried by exactly one founder can only be detected
in a founder that has its own assembly: cactus founders carry a median 52 private events against 3 for
short-read (17×, p = 2 × 10⁻³³), 98 of 151 short-read founders carry **zero** private coding events, and
assembly route alone explains 37 % of the variance in the private burden rate. A three-bin fit restricted
to the 80 cactus founders, where detection is homogeneous, is in the model registry; it is not shown here
because n = 80, the fit is singular, and the private bin still behaves like a detection measure.

⚠ The test is the **contrast**, not two separate p-values. "Rare is significant, common is not" is a
difference in significance, not a significant difference."""

code_freq = r'''
FQ=f"{PA}/coef_M2p_reml_x_obio1cobio1cbio1ccallratesv_re0obio1cgarden_refid_freqCDSrarecommon.csv"
def fload(coefp,scalep,vcovp):
    co=pd.read_csv(coefp).set_index("term"); sc=pd.read_csv(scalep).set_index("context")
    pv={f"sv_{c}":100/(sc.loc[c,"sd_rate"]*sc.loc[c,"med_called"]) for c in sc.index}
    return co,pv,sc,pd.read_csv(vcovp,index_col=0)
co,pv,sc,V=fload(FQ,f"{PA}/pooled_rate_scale_x_obio1cobio1cbio1ccallratesv_freqCDSrarecommon.csv",
                 FQ.replace("coef_","vcov_"))
def contrast(co,pv,V,a,b):
    ka,kb=pv[a],pv[b]
    d=co.loc[a,"estimate"]*ka-co.loc[b,"estimate"]*kb
    v=V.loc[a,a]*ka**2+V.loc[b,b]*kb**2-2*V.loc[a,b]*ka*kb
    z=d/np.sqrt(v); return d,z,2*stats.norm.sf(abs(z))
# the undivided coding-SV effect, with its interval, so the strata can be compared against the whole
UNS=CO.loc["sv_CDS","estimate"]*PV("sv_CDS")
UNSE=CO.loc["sv_CDS","se"]*PV("sv_CDS")
UNLO,UNHI=UNS-1.96*UNSE,UNS+1.96*UNSE
BINS=["rare","common"]; FC={"rare":"#6B5A8E","common":"#4A3B6B"}
fig,a=plt.subplots(figsize=(7.6,3.2)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
a.axvspan(UNLO,UNHI,color="#B04A00",alpha=0.10,lw=0,zorder=1)     # 95% CI of the undivided estimate
a.axvline(UNS,color="#B04A00",lw=1.3,ls=":",zorder=2)
for i,b in enumerate(BINS):
    t=f"sv_CDS_{b}"; e=co.loc[t,"estimate"]*pv[t]; se=co.loc[t,"se"]*pv[t]
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=FC[b],mfc=FC[b] if co.loc[t,"p"]<0.05 else "white",
               mew=1.2,ms=7,capsize=2.5,zorder=3)
    a.annotate(f"p = {co.loc[t,'p']:.3g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=7,color="0.35")
    if i<len(BINS)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.axhline(len(BINS)-0.5,color="0.75",lw=0.9,ls="--")               # separator: strata above, whole below
a.errorbar(UNS,len(BINS),xerr=1.96*UNSE,fmt="o",color="#B04A00",mew=1.2,ms=7,capsize=2.5,zorder=3)
a.annotate(f"p = {CO.loc['sv_CDS','p']:.3g}",xy=(0.985,len(BINS)),xycoords=("axes fraction","data"),
           ha="right",va="center",fontsize=7,color="#B04A00")
a.set_yticks(range(len(BINS)+1))
a.set_yticklabels(["CDS rare\n(MAC 2-12)","CDS common\n(MAC 13+)","all frequencies\ntogether"],fontsize=8)
a.get_yticklabels()[-1].set_color("#B04A00")
a.invert_yaxis(); a.set_ylim(len(BINS)+0.6,-0.6)
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.30*(x1-x0))
a.set_xlabel("% fitness per carried variant  (± 95% CI)",fontsize=8.5)
save(fig,"model_freq_strata",nb=True); plt.show()
print("FREQUENCY-STRATIFIED CDS BURDEN, all 231 founders\n")
print(pd.DataFrame([{"bin":b,
    "%/carried variant":f"{co.loc[f'sv_CDS_{b}','estimate']*pv[f'sv_CDS_{b}']:+.4f}",
    "95% CI":f"[{(co.loc[f'sv_CDS_{b}','estimate']-1.96*co.loc[f'sv_CDS_{b}','se'])*pv[f'sv_CDS_{b}']:+.3f}, "
             f"{(co.loc[f'sv_CDS_{b}','estimate']+1.96*co.loc[f'sv_CDS_{b}','se'])*pv[f'sv_CDS_{b}']:+.3f}]",
    "Wald p":f"{co.loc[f'sv_CDS_{b}','p']:.4g}",
    "median carried":f"{sc.loc[f'CDS_{b}','med_carried']:.0f}",
    "total load":f"{co.loc[f'sv_CDS_{b}','estimate']*pv[f'sv_CDS_{b}']*sc.loc[f'CDS_{b}','med_carried']:+.1f}%"}
    for b in BINS]).to_string(index=False))
print(f"\nall frequencies together (unstratified): {UNS:+.4f}%  [{UNLO:+.3f}, {UNHI:+.3f}]  "
      f"p = {CO.loc['sv_CDS','p']:.4g}")
d,z,p=contrast(co,pv,V,"sv_CDS_rare","sv_CDS_common")
print(f"\nCONTRAST rare vs common (per carried variant): difference {d:+.4f}   z = {z:+.2f}   p = {p:.4g}")
print("\nRare coding SVs are ~2.6x more deleterious per variant than common ones -- the predicted direction --")
print("but the contrast is NOT significant, and common variants carry ~3x more of the TOTAL load because")
print("founders carry far more of them. Both statements are true and answer different questions.")
'''

md_size = """## 8 · Size strata — is the effect carried by large SVs?

The same question along the other axis. Insertions and deletions stay **pooled**: which one a variant is
called depends on the Col-0 reference and is not a biological distinction — the same reason the
site-frequency spectra in Figure 1 are folded.

Bins are 50–99 / 100–499 / 500+ bp, the most balanced three-way split available in CDS (444 / 481 / 482
events). MAC ≥ 2 throughout, matching the model, so the singleton class cannot leak in. As in §7 the bottom
row and shaded band are the **unstratified** coefficient and its 95 % CI — where the effect falls with all
sizes held together.

⚠ Size is partly a proxy for insertion vs deletion at the *event* level — the largest CDS bin is 80 %
insertions, the smallest 38 %, because CDS insertions are ~3× larger here. At the *carried* level, which is
what the burden counts, every bin is 50–53 % insertions, so pooling is not hiding a class contrast.

Because three bins each hold a third of the variants, each is far less precise than the undivided term. The
cell therefore reports an **omnibus test** of whether the bins differ at all, and a **power calculation** —
both needed to read "only the large bin is significant" correctly."""

code_size = r'''
SZ=f"{PA}/coef_M2p_reml_x_obio1cobio1cbio1c_re0obio1cgarden_refid_freqCDSsmallmediumlarge.csv"
co,pv,sc,V=fload(SZ,f"{PA}/pooled_rate_scale_x_obio1cobio1cbio1c_freqCDSsmallmediumlarge.csv",
                 SZ.replace("coef_","vcov_"))
BINS=["small","medium","large"]; SC3={"small":"#C6BCD8","medium":"#6B5A8E","large":"#4A3B6B"}
RNG={"small":"50-99 bp","medium":"100-499 bp","large":"500+ bp"}
e=np.array([co.loc[f"sv_CDS_{b}","estimate"]*pv[f"sv_CDS_{b}"] for b in BINS])
se=np.array([co.loc[f"sv_CDS_{b}","se"]*pv[f"sv_CDS_{b}"] for b in BINS])
fig,a=plt.subplots(figsize=(7.6,3.5)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
a.axvspan(UNLO,UNHI,color="#B04A00",alpha=0.10,lw=0,zorder=1)     # 95% CI of the undivided estimate
a.axvline(UNS,color="#B04A00",lw=1.3,ls=":",zorder=2)
for i,b in enumerate(BINS):
    t=f"sv_CDS_{b}"
    a.errorbar(e[i],i,xerr=1.96*se[i],fmt="o",color=SC3[b],mfc=SC3[b] if co.loc[t,"p"]<0.05 else "white",
               mew=1.2,ms=7,capsize=2.5,zorder=3)
    a.annotate(f"p = {co.loc[t,'p']:.3g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=7,color="0.35")
    if i<len(BINS)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.axhline(2.5,color="0.75",lw=0.9,ls="--")                        # separator: strata above, whole below
a.errorbar(UNS,3,xerr=1.96*UNSE,fmt="o",color="#B04A00",mew=1.2,ms=7,capsize=2.5,zorder=3)
a.annotate(f"p = {CO.loc['sv_CDS','p']:.3g}",xy=(0.985,3),xycoords=("axes fraction","data"),
           ha="right",va="center",fontsize=7,color="#B04A00")
a.set_yticks(range(4)); a.set_yticklabels([f"CDS {RNG[b]}" for b in BINS]+["all sizes\ntogether"],fontsize=8.5)
a.get_yticklabels()[-1].set_color("#B04A00")
a.invert_yaxis(); a.set_ylim(3.6,-0.6)
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.30*(x1-x0))
a.set_xlabel("% fitness per carried variant  (± 95% CI)",fontsize=8.5)
save(fig,"model_size_strata",nb=True); plt.show()
print("SIZE-STRATIFIED CDS BURDEN (pooled ins+del, MAC>=2)\n")
print(pd.DataFrame([{"bin":b,"range":RNG[b],"%/carried variant":f"{e[i]:+.4f}",
    "95% CI":f"[{e[i]-1.96*se[i]:+.3f}, {e[i]+1.96*se[i]:+.3f}]",
    "Wald p":f"{co.loc[f'sv_CDS_{b}','p']:.4g}",
    "median carried":f"{sc.loc[f'CDS_{b}','med_carried']:.0f}"} for i,b in enumerate(BINS)]
    ).to_string(index=False))
print(f"\nall sizes together (unstratified): {UNS:+.4f}%  [{UNLO:+.3f}, {UNHI:+.3f}]  "
      f"p = {CO.loc['sv_CDS','p']:.4g}")
print("\nCONTRASTS (per carried variant):")
for a_,b_ in (("small","medium"),("medium","large"),("small","large")):
    d,z,p=contrast(co,pv,V,f"sv_CDS_{a_}",f"sv_CDS_{b_}")
    print(f"  {a_:<7} vs {b_:<7} difference {d:+.4f}   z = {z:+.2f}   p = {p:.4g}")
k=np.array([pv[f"sv_CDS_{b}"] for b in BINS]); T3=[f"sv_CDS_{b}" for b in BINS]
Vs=V.loc[T3,T3].to_numpy()*np.outer(k,k); Lc=np.array([[1,-1,0],[0,1,-1]])
dd=Lc@e; W=dd@np.linalg.solve(Lc@Vs@Lc.T,dd)
print(f"\nOMNIBUS TEST -- do the three bins differ at all?  chi2 = {W:.3f} on 2 df, p = {stats.chi2.sf(W,2):.4g}")
print("\nPOWER -- if all three bins truly had the unstratified effect, how often would each")
print("         individually reach p < 0.05?")
pw=[stats.norm.sf(1.96-abs(UNS)/s)+stats.norm.cdf(-1.96-abs(UNS)/s) for s in se]
for b,s,w in zip(BINS,se,pw): print(f"  {b:<7} SE {s:.3f}  ->  power {w:.0%}")
print(f"  expected number of bins reaching p<0.05 under that null: {sum(pw):.1f} of 3 -- and 1 was observed.")
print("\nEvery bin's CI contains the other bins' point estimates and the unstratified value, the omnibus")
print("test is far from significant, and one-of-three significant is the EXPECTED outcome even if size is")
print("irrelevant. The estimates order as theory predicts, but the study cannot resolve them.")
'''
