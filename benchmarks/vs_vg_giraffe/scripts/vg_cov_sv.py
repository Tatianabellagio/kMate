#!/usr/bin/env python3
"""Option 2 vg-SV allele frequency from pack node-COVERAGE (no remap; reuses the
existing cactus/d2 pack). Two stages:

  stage1 (--emit-nodelist): parse a `vg deconstruct -a` VCF of the d2 graph (AT =
    ref+alt node traversals, in the pack's node-ID space). Subset to SVs, bridge each
    to the panel `svidx` by REF/ALT SEQUENCE (node-ID-version-independent), and write
    (a) a node-list for `vg pack -d -N`, (b) a JSON of per-svidx node sets.

  stage2 (--compute-af): read the `vg pack -d` node-coverage table + the JSON, compute
    per-SV AF and write (svidx, est) for score_sv to consume via --vg-cov.

Unified boundary-normalized AF (works for INS and DEL):
  boundary = ref_nodes ∩ alt_nodes (shared anchors ~ local depth)
  alt_spec = alt_nodes - ref_nodes ; ref_spec = ref_nodes - alt_nodes
  alt_support = mean_cov(alt_spec) if alt_spec else max(boundary_cov - mean_cov(ref_spec), 0)
  ref_support = mean_cov(ref_spec) if ref_spec else max(boundary_cov - mean_cov(alt_spec), 0)
  AF = alt_support / (alt_support + ref_support)   # INS->altcov/boundary; DEL->1-refcov/boundary
"""
import argparse, gzip, json, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_sv import build_panel_index, seq_key  # noqa: E402

NODE = re.compile(r"\d+")


def parse_at(at_field):
    """AT=trav0,trav1,... ; each trav like '>123>124<125'. Return list of node-id sets
    (ordered list preserved for boundary detection via first/last)."""
    travs = []
    for t in at_field.split(","):
        travs.append([int(x) for x in NODE.findall(t)])
    return travs


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def stage1(decon_vcf, panel_sv_vcf, svlen, out_nodelist, out_json):
    key2idx = build_panel_index(panel_sv_vcf, svlen)
    rec = {}                      # svidx -> {ref:[...], alt:[...]}
    nodes = set()
    n_at_missing = n_unmatched = 0
    with opener(decon_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            at = info.get("AT", "")
            travs = parse_at(at) if at and at != "." else []
            for ai, alt in enumerate(alts):
                if abs(len(alt) - len(ref)) < svlen:
                    continue
                sv = key2idx.get(seq_key(chrom, pos, ref, alt))
                if sv is None:
                    n_unmatched += 1
                    continue
                if len(travs) < ai + 2:        # need ref (0) + this alt (ai+1)
                    n_at_missing += 1
                    continue
                ref_nodes, alt_nodes = travs[0], travs[ai + 1]
                rec[int(sv)] = {"ref": ref_nodes, "alt": alt_nodes}
                nodes.update(ref_nodes); nodes.update(alt_nodes)
    Path(out_nodelist).write_text("\n".join(map(str, sorted(nodes))) + "\n")
    Path(out_json).write_text(json.dumps(rec))
    print(f"[stage1] matched {len(rec):,} SVs -> svidx | unmatched(seq) {n_unmatched:,} | "
          f"AT-missing {n_at_missing:,} | {len(nodes):,} nodes for pack -d -N", file=sys.stderr)


def load_node_cov(pack_table):
    """vg pack -d table -> dict node_id -> mean per-base coverage. Header names vary;
    detect node.id and coverage columns, average coverage over each node's positions."""
    s, c = {}, {}                 # sum, count per node
    with opener(pack_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        hl = [h.lower() for h in header]
        nid = next(i for i, h in enumerate(hl) if "node" in h and "id" in h)
        cov = next(i for i, h in enumerate(hl) if "cov" in h)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            n = int(f[nid]); v = float(f[cov])
            s[n] = s.get(n, 0.0) + v; c[n] = c.get(n, 0) + 1
    return {n: s[n] / c[n] for n in s}


def mean_cov(nodes, cov):
    vals = [cov[n] for n in nodes if n in cov]
    return float(np.mean(vals)) if vals else np.nan


def stage2(pack_table, in_json, out_tsv):
    cov = load_node_cov(pack_table)
    rec = json.loads(Path(in_json).read_text())
    rows = []
    for sv, d in rec.items():
        ref_n, alt_n = set(d["ref"]), set(d["alt"])
        boundary = ref_n & alt_n
        alt_spec, ref_spec = alt_n - ref_n, ref_n - alt_n
        bcov = mean_cov(boundary, cov)
        ac = mean_cov(alt_spec, cov); rc = mean_cov(ref_spec, cov)
        alt_sup = ac if not np.isnan(ac) else (max(bcov - rc, 0.0) if not np.isnan(bcov) and not np.isnan(rc) else np.nan)
        ref_sup = rc if not np.isnan(rc) else (max(bcov - ac, 0.0) if not np.isnan(bcov) and not np.isnan(ac) else np.nan)
        if np.isnan(alt_sup) or np.isnan(ref_sup) or (alt_sup + ref_sup) <= 0:
            af = np.nan
        else:
            af = min(max(alt_sup / (alt_sup + ref_sup), 0.0), 1.0)
        rows.append((int(sv), af))
    with open(out_tsv, "w") as fh:
        fh.write("svidx\test\n")
        for sv, af in sorted(rows):
            fh.write(f"{sv}\t{af}\n")
    finite = sum(1 for _, af in rows if not np.isnan(af))
    print(f"[stage2] AF for {len(rows):,} SVs ({finite:,} finite) -> {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--svlen", type=int, default=50)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s1 = sub.add_parser("emit-nodelist")
    s1.add_argument("--decon-vcf", required=True)
    s1.add_argument("--panel-sv-vcf", required=True)
    s1.add_argument("--out-nodelist", required=True)
    s1.add_argument("--out-json", required=True)
    s2 = sub.add_parser("compute-af")
    s2.add_argument("--pack-table", required=True)
    s2.add_argument("--in-json", required=True)
    s2.add_argument("--out-tsv", required=True)
    a = ap.parse_args()
    if a.cmd == "emit-nodelist":
        stage1(a.decon_vcf, a.panel_sv_vcf, a.svlen, a.out_nodelist, a.out_json)
    else:
        stage2(a.pack_table, a.in_json, a.out_tsv)
