#!/usr/bin/env python
"""Render the CAM5 founder pangenome graph coloured by exon, via graphviz sfdp.

Maps every graph node to a Chr2 position (ref nodes by walking the reference path;
alt/variant nodes by nearest reference node-id), classifies it as exon A/B/C /
intron / UTR / flank, writes a coloured DOT, and lays it out with sfdp (compact,
Bandage-like). Established tool = graphviz (system /usr/bin/sfdp).
"""
import re, subprocess

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/cam5_replication/pang_cam5"
GFA = f"{HERE}/cam5_full.gfa"
REG0 = 11531800
# TAIR10 AT2G27030 features (representative 3-exon + UTR split)
FEAT = [("exon A (5′)", 11532004, 11532068, "#9ecae1"),   # 5'UTR
        ("exon A (5′)", 11532069, 11532144, "#1f77b4"),   # CDS
        ("exon B (mid)", 11532687, 11533060, "#2ca02c"),  # CDS
        ("exon B (mid)", 11533061, 11533276, "#98df8a"),  # 3'UTR(.1)
        ("exon C (3′)", 11534077, 11534333, "#ff7f0e")]   # 3' exon
INTRON, FLANK = "#e0e0e0", "#fafafa"
LEG = [("exon A 5′UTR", "#9ecae1"), ("exon A CDS", "#1f77b4"),
       ("exon B CDS", "#2ca02c"), ("exon B 3′UTR", "#98df8a"),
       ("exon C (3′, GEA)", "#ff7f0e"), ("intron", "#e0e0e0"), ("flank", "#fafafa")]

def colour(pos):
    if pos is None: return FLANK
    for _, s, e, c in FEAT:
        if s <= pos <= e: return c
    return INTRON if 11532004 <= pos <= 11534333 else FLANK

# parse GFA
S, L, refpath = {}, [], []
for ln in open(GFA):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "S": S[f[1]] = len(f[2])
    elif f[0] == "L": L.append((f[1], f[3]))
    elif f[0] == "P" and f[1] == "Chr2":
        refpath = [x[:-1] for x in f[2].split(",")]

# ref node -> genomic pos
pos = {}; off = 0
for n in refpath:
    pos[n] = REG0 + off; off += S[n]
# alt/other nodes -> nearest ref node by integer id (vg numbers ~positionally)
ref_ids = sorted((int(n), pos[n]) for n in refpath)
import bisect
rk = [r[0] for r in ref_ids]
for n in S:
    if n in pos: continue
    i = bisect.bisect_left(rk, int(n))
    cand = [ref_ids[j][1] for j in (i-1, i) if 0 <= j < len(ref_ids)]
    pos[n] = cand[0] if cand else None

# write coloured DOT
with open(f"{HERE}/cam5_exons.dot", "w") as d:
    d.write('digraph G {\n  graph [bgcolor=white];\n')
    d.write('  node [shape=box, style="filled", fixedsize=true, height=0.18, '
            'penwidth=0.4, label=""];\n  edge [arrowsize=0.3, color="#888888", penwidth=0.6];\n')
    for n, l in S.items():
        c = colour(pos.get(n))
        w = max(0.10, min(l, 60) / 60 * 0.9)            # width ~ node length (capped)
        d.write(f'  {n} [fillcolor="{c}", width={w:.2f}];\n')
    for a, b in L:
        d.write(f'  {a} -> {b};\n')
    # legend as a cluster
    d.write('  subgraph cluster_leg { label="exon colour"; fontsize=11; color=white;\n')
    for i, (nm, c) in enumerate(LEG):
        d.write(f'    L{i} [label="{nm}", fillcolor="{c}", shape=box, '
                f'width=1.4, height=0.25, fontsize=9];\n')
    for i in range(len(LEG)-1):
        d.write(f'    L{i} -> L{i+1} [style=invis];\n')
    d.write('  }\n}\n')

n_by = {}
for n in S: n_by[colour(pos.get(n))] = n_by.get(colour(pos.get(n)), 0) + 1
print("nodes by colour:", {dict(LEG+[("?","?")]).get(c, c): v for c, v in n_by.items()} if False else
      {nm: sum(1 for n in S if colour(pos.get(n)) == c) for nm, c in LEG})

# render with sfdp (compact, Bandage-like) — system graphviz
out = f"{HERE}/cam5_exons_graph.png"
r = subprocess.run(["/usr/bin/sfdp", "-Tpng", "-Gdpi=130", "-Goverlap=prism",
                    "-Gsplines=true", f"{HERE}/cam5_exons.dot", "-o", out],
                   capture_output=True, text=True)
print("sfdp:", "OK ->", out if r.returncode == 0 else r.stderr[:300])
