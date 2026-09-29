#!/usr/bin/env python3
"""vg-SV allele frequency matched to the panel by SNARL TRAVERSAL node-path (the
decomposition-correct key), requires vg run on the NON-d2 graph (panel node space).

Panel SV records (panel_sv VCF, in svidx order) carry their atomic traversal in the
ID column, e.g.  Chr1-2-INS->4410>4416>4417>...  -> node-path (4410,4416,4417,...).
`vg call -a` emits per snarl an AT=trav0(ref),trav1,trav2,... and AD per traversal;
AF(trav_i) = AD_i / sum(AD). We key every vg alt-traversal by its node-path and look
up each panel SV's traversal -> vg's AF for exactly that allele. No REF/ALT sequence
matching, no bcftools norm (which mangles SV decomposition).

Usage:
  vg_nond2_sv_af.py --vg-vcf <nond2_vgcalla.vcf> --panel-sv-vcf <panel_sv.vcf.gz> \
      --out <svidx_est.tsv> [--svlen 50]
"""
import argparse, gzip, re, sys
import numpy as np, pandas as pd

NODES = re.compile(r"\d+")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def path_key(s):
    """node-path -> tuple of ints (orientation-agnostic; order preserved)."""
    return tuple(int(x) for x in NODES.findall(s))


def panel_sv_traversals(panel_sv_vcf, svlen):
    """svidx (row in SV subset) -> traversal key, from the ID column's >node>node path."""
    out, idx = {}, 0
    with opener(panel_sv_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if abs(len(f[4]) - len(f[3])) < svlen:
                continue
            # ID like Chr1-2-INS->4410>4416>...  -> take node path after the last '-'
            trav = f[2].split("-")[-1]
            out[idx] = path_key(trav)
            idx += 1
    return out


def vg_traversal_af(vg_vcf):
    """node-path key -> AF, from vg call -a (AT alt-traversals + AD)."""
    af = {}
    with opener(vg_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            at = info.get("AT", "")
            if not at:
                continue
            travs = at.split(",")                       # trav0 = ref, 1.. = alts
            fd = dict(zip(f[8].split(":"), f[9].split(":")))
            ad = fd.get("AD", "")
            try:
                adv = [int(x) for x in ad.split(",")]
            except ValueError:
                continue
            tot = sum(adv)
            if tot <= 0 or len(adv) != len(travs):
                continue
            for i in range(1, len(travs)):              # alt traversals
                af[path_key(travs[i])] = adv[i] / tot
    return af


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-vcf", required=True)
    ap.add_argument("--panel-sv-vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--svlen", type=int, default=50)
    a = ap.parse_args()

    panel = panel_sv_traversals(a.panel_sv_vcf, a.svlen)
    vgaf = vg_traversal_af(a.vg_vcf)
    rows = [(sv, vgaf[key]) for sv, key in panel.items() if key in vgaf]
    out = pd.DataFrame(rows, columns=["svidx", "est"]).sort_values("svidx")
    out.to_csv(a.out, sep="\t", index=False)
    print(f"[vg non-d2] panel SVs={len(panel):,} | vg alt-traversals={len(vgaf):,} | "
          f"matched by traversal: {len(out):,} ({len(out)/max(len(panel),1):.1%}) -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
