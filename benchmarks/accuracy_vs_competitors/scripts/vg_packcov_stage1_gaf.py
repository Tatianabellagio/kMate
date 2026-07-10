#!/usr/bin/env python3
"""Stage 1 of the raw-pack-coverage vg-SV estimator, built from `vg paths` output
(NOT vg deconstruct, which is too slow / treats _alt_ paths as samples).

Inputs (all in the construct graph's node space = the pack's node space):
  --alt-gaf  : `vg paths -x g -a -A`  -> per _alt_ path: name, len, NODE WALK (col6)
  --alt-fa   : `vg paths -x g -a -F`  -> per _alt_ path: name, SEQUENCE
  --chr1-gaf : `vg paths -x g -p Chr1 -A` -> the reference path's full node walk
  --panel-sv-vcf : panel SV VCF (gives svidx in SV-mask order + ALT sequences)

For each _alt_ path whose sequence matches a panel SV ALT allele:
  alt_nodes = the path's node walk
  anchors   = alt_nodes that lie on the Chr1 reference walk (shared boundary nodes)
  ref_nodes = the Chr1 walk segment spanning [min..max anchor]  (the ref allele)
Emits rec.json {svidx: {ref:[...], alt:[...]}} + node list, consumed by
vg_cov_sv.py `compute-af` (boundary-normalized AF, works for INS and DEL).
"""
import argparse, gzip, json, re, sys
from collections import defaultdict
from pathlib import Path

NODE = re.compile(r"\d+")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def build_panel_altseq(panel_sv_vcf, svlen):
    """ALT-sequence -> ordered list of svidx (SV-mask order == truth order).
    Position (from VCF) kept for greedy disambiguation of repeated sequences."""
    by_seq = defaultdict(list)
    idx = 0
    with opener(panel_sv_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            for alt in alts:
                if abs(len(alt) - len(ref)) < svlen:
                    continue
                by_seq[alt.upper()].append((idx, pos))
                idx += 1
    print(f"[stage1] panel SV alleles: {idx:,} ({len(by_seq):,} distinct ALT seqs)",
          file=sys.stderr)
    return by_seq, idx


def parse_fasta(fa):
    name, seq, out = None, [], {}
    with opener(fa) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if name is not None:
                    out[name] = "".join(seq).upper()
                name, seq = ln[1:].split()[0], []
            else:
                seq.append(ln)
    if name is not None:
        out[name] = "".join(seq).upper()
    return out


def parse_gaf_nodes(gaf):
    """name -> ordered node-id list (from GAF path col 6)."""
    out = {}
    with opener(gaf) as fh:
        for ln in fh:
            if ln.startswith("@"):
                continue
            f = ln.split("\t")
            if len(f) < 6:
                continue
            out[f[0]] = [int(x) for x in NODE.findall(f[5])]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alt-gaf", required=True)
    ap.add_argument("--alt-fa", required=True)
    ap.add_argument("--chr1-gaf", required=True)
    ap.add_argument("--panel-sv-vcf", required=True)
    ap.add_argument("--svlen", type=int, default=50)
    ap.add_argument("--out-nodelist", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    by_seq, n_sv = build_panel_altseq(a.panel_sv_vcf, a.svlen)
    fa = parse_fasta(a.alt_fa)
    print(f"[stage1] {len(fa):,} _alt_ path sequences", file=sys.stderr)
    gaf = parse_gaf_nodes(a.alt_gaf)
    print(f"[stage1] {len(gaf):,} _alt_ path node walks", file=sys.stderr)

    # Chr1 reference walk -> node -> order index
    chr1 = parse_gaf_nodes(a.chr1_gaf)
    chr1_nodes = max(chr1.values(), key=len)  # the single Chr1 path
    chr1_idx = {n: i for i, n in enumerate(chr1_nodes)}
    print(f"[stage1] Chr1 walk: {len(chr1_nodes):,} nodes", file=sys.stderr)

    # candidate _alt_ paths whose seq matches a panel SV ALT, grouped by seq
    # so we can greedily disambiguate repeated sequences by genomic order.
    cand = defaultdict(list)  # seq -> list of (anchor_order, name)
    for name, seq in fa.items():
        if len(seq) < a.svlen or seq not in by_seq:
            continue
        nodes = gaf.get(name, [])
        anchors = [chr1_idx[n] for n in nodes if n in chr1_idx]
        ord_key = min(anchors) if anchors else 1 << 60
        cand[seq].append((ord_key, name))

    rec, nodes_needed = {}, set()
    n_match = n_collide = n_noanchor = 0
    for seq, svlist in by_seq.items():
        paths = sorted(cand.get(seq, []))
        svlist_sorted = sorted(svlist, key=lambda t: t[1])  # by pos
        if len(paths) > 1 or len(svlist_sorted) > 1:
            n_collide += 1
        for (svidx, _pos), (_ord, name) in zip(svlist_sorted, paths):
            alt_nodes = gaf[name]
            anchors = [n for n in alt_nodes if n in chr1_idx]
            if anchors:
                lo = min(chr1_idx[n] for n in anchors)
                hi = max(chr1_idx[n] for n in anchors)
                ref_nodes = chr1_nodes[lo:hi + 1]
            else:
                ref_nodes = []
                n_noanchor += 1
            rec[int(svidx)] = {"ref": list(map(int, ref_nodes)),
                               "alt": list(map(int, alt_nodes))}
            nodes_needed.update(ref_nodes); nodes_needed.update(alt_nodes)
            n_match += 1

    Path(a.out_nodelist).write_text("\n".join(map(str, sorted(nodes_needed))) + "\n")
    Path(a.out_json).write_text(json.dumps(rec))
    print(f"[stage1] matched {n_match:,}/{n_sv:,} panel SVs to _alt_ paths "
          f"({100*n_match/n_sv:.1f}%); {n_collide:,} seqs w/ repeats; "
          f"{n_noanchor:,} no-anchor; {len(nodes_needed):,} nodes for pack -d",
          file=sys.stderr)


if __name__ == "__main__":
    main()
