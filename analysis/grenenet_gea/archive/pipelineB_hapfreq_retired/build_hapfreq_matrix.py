#!/usr/bin/env python
"""STEP 2: project per-block founder h -> per-sample haplotype-frequency matrix.

For every cohort sample (window-mode run, results/grenenet_kmate_window/*_Chr{N}.h_blocks_per_chrom.npz)
and every dynld-K500 unit, compute the frequency of each haplotype cluster:
    freq_{u,c} = sum_{f in cluster c} h_{u,f}   = bincount(global_hap_id, weights=h_block)
using the STEP-1 membership (hap_membership/{chrlc}_hapmemb_K500.npz: labels U x231, hap_offset).
Clusters partition all 231 founders so per-unit freqs sum to 1 (h on simplex).

Genome-wide hap vector = Chr1..Chr5 concatenated (67,930 haplotypes). Founder order in the
membership labels == founder order in h_blocks (var_pa meta 'founders'); verified.

Chunked for array jobs: --start/--end slice the sorted sample list; writes a chunk .npy.
Run gather mode (--gather) to stack chunks into the final matrix.

  python build_hapfreq_matrix.py --start 0 --end 271 --out-dir <chunkdir>      # one array task
  python build_hapfreq_matrix.py --gather --chunk-dir <chunkdir> --out-dir <final>

Env: kmate.
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd

WIN = "results/grenenet_kmate_window"
HM = "results/grenenet_gea/blocks_mcf90/hap_membership"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
# module-level config (overridable by callers, e.g. build_hapfreq_p0_seedmix, and by main()):
#   MEMB_TAG selects the membership file suffix ({chrlc}_hapmemb_{MEMB_TAG}.npz).
#   H_SOURCE = "window" projects the per-block window h ({ch}_h_blocks, locked to that run's
#   block grid); "global" projects the genome-global founder h ({ch}_global_h, a single
#   231-vector) onto the membership grid -> works with ANY block partition (e.g. clq0.9),
#   no window-EM re-run. window≈global here (recombination rare; corr 0.9998).
MEMB_TAG = "K500"
H_SOURCE = "window"


def load_mem(registry=False):
    """Return (mem, reg_or_None, n_haps_total). mem: per-chrom labels/hap_offset/nh.
    Global hap index runs Chr1..Chr5 concatenated. registry=True also builds the df
    (one row per haplotype) — only needed at gather; skip it for fast chunk tasks."""
    mem = {}
    reg_rows = []
    gbase = 0
    for ch in CHROMS:
        M = np.load(f"{HM}/{ch.lower()}_hapmemb_{MEMB_TAG}.npz", allow_pickle=True)
        labels = M["labels"]                       # (U, 231) int16
        off = M["hap_offset"].astype(np.int64)     # (U+1,)
        nh = int(off[-1])
        mem[ch] = dict(labels=labels, off=off, nh=nh,
                       start=M["unit_start"], end=M["unit_end"],
                       founders=M["founders"].astype(str))
        if not registry:
            gbase += nh
            continue
        U = len(labels)
        for u in range(U):
            k = int(off[u + 1] - off[u])
            fr = np.bincount(labels[u], minlength=k) / labels.shape[1]
            ne = 1.0 / np.sum((np.bincount(labels[u], minlength=k) / labels.shape[1]) ** 2)
            for c in range(k):
                reg_rows.append((gbase + int(off[u]) + c, ch, u,
                                 int(M["unit_start"][u]), int(M["unit_end"][u]), c,
                                 int((labels[u] == c).sum()), round(float(fr[c]), 5),
                                 int(M["unit_nvar"][u]), int(M["unit_kmers"][u]),
                                 bool(M["unit_covered"][u]), round(float(ne), 3)))
        gbase += nh
    reg = pd.DataFrame(reg_rows, columns=[
        "hap_id", "chrom", "unit_idx", "unit_start", "unit_end", "cluster",
        "n_founders", "panel_freq", "unit_nvar", "panel_kmers", "covered", "unit_n_eff"]) if registry else None
    return mem, reg, gbase


def project_sample(sample, mem):
    """Concatenate Chr1..Chr5 haplotype freqs for one sample -> (n_haps,) float32.

    H_SOURCE="window": project the per-block window h ({ch}_h_blocks); requires the
    membership grid to match that run's block grid (U equal).
    H_SOURCE="global": broadcast the genome-global founder h ({ch}_global_h, one 231-vector)
    across the membership's U units, then sum within cluster. Partition-independent — this
    is how a clq0.9 hapfreq is obtained from the SAME global h that produced the per-variant
    AF, with no window-EM re-run."""
    parts = []
    for ch in CHROMS:
        f = f"{WIN}/{sample}_{ch}.h_blocks_per_chrom.npz"
        H = np.load(f, allow_pickle=True)
        m = mem[ch]
        U = m["labels"].shape[0]
        if H_SOURCE == "global":
            gh = np.asarray(H[f"{ch}_global_h"], dtype=np.float64)     # (231,)
            if gh.shape[0] != m["labels"].shape[1]:
                raise ValueError(f"{sample} {ch}: global_h {gh.shape} vs labels {m['labels'].shape}")
            hb = np.broadcast_to(gh, (U, gh.shape[0]))                 # (U, 231), block-invariant
        else:
            hb = H[f"{ch}_h_blocks"]                                   # (U, 231) float32
            if hb.shape[0] != U:
                raise ValueError(f"{sample} {ch}: h_blocks {hb.shape} vs labels {m['labels'].shape}")
        # global-within-chrom hap id of each (unit, founder)
        glob = (m["off"][:-1][:, None] + m["labels"]).ravel()
        hf = np.bincount(glob, weights=np.asarray(hb).ravel(), minlength=m["nh"])
        parts.append(hf.astype(np.float32))
    return np.concatenate(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=0, help="0 = to end")
    ap.add_argument("--out-dir", default="results/grenenet_gea/hapfreq")
    ap.add_argument("--chunk-dir", default=None)
    ap.add_argument("--gather", action="store_true")
    ap.add_argument("--memb-tag", default=MEMB_TAG, help="membership suffix, e.g. clq90")
    ap.add_argument("--h-source", default=H_SOURCE, choices=["window", "global"])
    a = ap.parse_args()
    globals()["MEMB_TAG"] = a.memb_tag
    globals()["H_SOURCE"] = a.h_source
    os.makedirs(a.out_dir, exist_ok=True)

    samples = sorted(os.path.basename(p).replace("_Chr1.h_blocks_per_chrom.npz", "")
                     for p in glob.glob(f"{WIN}/*_Chr1.h_blocks_per_chrom.npz"))

    if a.gather:
        cd = a.chunk_dir or a.out_dir
        chunks = sorted(glob.glob(f"{cd}/chunk_*.npy"),
                        key=lambda p: int(os.path.basename(p).split("_")[1]))
        mats = [np.load(c) for c in chunks]
        smats = []
        for c in chunks:
            with open(c.replace(".npy", ".samples")) as fh:
                smats += [l.strip() for l in fh if l.strip()]
        mat = np.vstack(mats)
        assert mat.shape[0] == len(smats) == len(samples), (mat.shape, len(smats), len(samples))
        assert smats == samples, "sample order mismatch on gather"
        mem, reg, nh = load_mem(registry=True)
        assert mat.shape[1] == nh == len(reg), (mat.shape, nh, len(reg))
        np.save(f"{a.out_dir}/hapfreq_matrix.npy", mat.astype(np.float32))
        with open(f"{a.out_dir}/hapfreq_samples.txt", "w") as fh:
            fh.write("\n".join(samples) + "\n")
        reg.to_csv(f"{a.out_dir}/hapfreq_registry.csv", index=False)
        np.save(f"{a.out_dir}/hapfreq_p0_panel.npy", reg.panel_freq.values.astype(np.float32))
        print(f"[gather] matrix {mat.shape} ({mat.nbytes/1e9:.2f} GB) -> {a.out_dir}/hapfreq_matrix.npy")
        print(f"[gather] registry {len(reg):,} haplotypes; samples {len(samples):,}")
        return

    end = a.end if a.end else len(samples)
    sub = samples[a.start:end]
    mem, _, nh = load_mem(registry=False)
    print(f"[chunk {a.start}:{end}] {len(sub)} samples, {nh:,} haplotypes", flush=True)
    out = np.zeros((len(sub), nh), dtype=np.float32)
    for i, s in enumerate(sub):
        out[i] = project_sample(s, mem)
        if (i + 1) % 50 == 0:
            print(f"  ...{i+1}/{len(sub)}", flush=True)
    cpath = f"{a.out_dir}/chunk_{a.start}_{end}.npy"
    np.save(cpath, out)
    with open(cpath.replace(".npy", ".samples"), "w") as fh:
        fh.write("\n".join(sub) + "\n")
    print(f"[chunk {a.start}:{end}] wrote {out.shape} -> {cpath}")


if __name__ == "__main__":
    main()
