"""Plot cactus 82 vs xwu 231 SNP overlap (Venn-style summary)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

DATA  = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/cactus_panel_overlap/data")
PLOTS = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/plots")
PLOTS.mkdir(exist_ok=True)


def load_pos(path):
    s = set()
    with open(path) as f:
        for line in f:
            c, p = line.rstrip("\n").split("\t")
            s.add((c, p))
    return s


print("Loading positions...")
xwu231 = load_pos(DATA / "xwu231_positions.tsv")
xwu80  = load_pos(DATA / "xwu80_positions.tsv")
xwu82  = load_pos(DATA / "xwu82_positions.tsv")
cactus = load_pos(DATA / "cactus82_snp_positions.tsv")
print(f"  xwu231: {len(xwu231):,}, xwu82: {len(xwu82):,}, cactus SNPs: {len(cactus):,}")

# 1) cactus vs xwu231 Venn
fig, ax = plt.subplots(figsize=(7, 6))
v = venn2([cactus, xwu231], set_labels=("cactus 82\nSNPs", "xwu 231\nSNPs"), ax=ax)
for t in v.set_labels: t.set_fontsize(13)
for t in v.subset_labels:
    if t: t.set_fontsize(11)
shared = cactus & xwu231
ax.set_title(f"Cactus 82-founder SNPs ↔ xwu GrENE-Net 231 SNPs\n"
             f"shared: {len(shared):,}  |  cactus-only: {len(cactus - xwu231):,}  |  "
             f"xwu-only: {len(xwu231 - cactus):,}", fontsize=11)
plt.tight_layout()
plt.savefig(PLOTS / "cactus_vs_xwu231_snp_venn.png", bbox_inches="tight", dpi=200)
plt.close()
print(f"Wrote {PLOTS}/cactus_vs_xwu231_snp_venn.png")

# 2) cactus vs xwu82 (subset to overlapping accessions) Venn
fig, ax = plt.subplots(figsize=(7, 6))
v = venn2([cactus, xwu82], set_labels=("cactus 82\nSNPs", "xwu 82-subset\nSNPs"), ax=ax)
for t in v.set_labels: t.set_fontsize(13)
for t in v.subset_labels:
    if t: t.set_fontsize(11)
shared = cactus & xwu82
ax.set_title(f"Apples-to-apples: same 82 founders\n"
             f"shared: {len(shared):,}  |  cactus-only: {len(cactus - xwu82):,}  |  "
             f"xwu82-only: {len(xwu82 - cactus):,}", fontsize=11)
plt.tight_layout()
plt.savefig(PLOTS / "cactus_vs_xwu82_snp_venn.png", bbox_inches="tight", dpi=200)
plt.close()
print(f"Wrote {PLOTS}/cactus_vs_xwu82_snp_venn.png")

# 3) Three-way: cactus, xwu80 (current ref panel SNPs), xwu231 (target panel includes these)
fig, ax = plt.subplots(figsize=(8, 7))
v = venn3([cactus, xwu80, xwu231 - xwu80],
          set_labels=("cactus SNPs", "xwu in ref_80\n(current backbone)",
                      "xwu polymorphic\nonly in 151 missing"), ax=ax)
for t in v.set_labels: t.set_fontsize(11)
for t in v.subset_labels:
    if t: t.set_fontsize(10)
ax.set_title("Imputation reference-panel design\n"
             "Current ref_80 SNPs = xwu 80-subset only (yellow + green ∩ red)\n"
             "Adding cactus would add the green-only region to the backbone", fontsize=10)
plt.tight_layout()
plt.savefig(PLOTS / "imputation_panel_design_3way.png", bbox_inches="tight", dpi=200)
plt.close()
print(f"Wrote {PLOTS}/imputation_panel_design_3way.png")

print("\nSummary table:")
print(f"  Current ref_80 SNP backbone  :  {len(xwu80):,} positions")
print(f"  Adding cactus SNPs (union)   :  {len(xwu80 | cactus):,} positions  "
      f"(+{100*len(xwu80 | cactus)/len(xwu80) - 100:.1f}%)")
print(f"  Cactus-only sites added      :  {len(cactus - xwu80):,}")
print(f"  Cactus pangenome covers      :  {100*len(cactus & xwu80)/len(xwu80):.1f}% "
      f"of current ref_80 backbone")
