"""`kmate selftest` — run the bundled tiny fixture end-to-end and check it.

A freshly-installed user runs `kmate selftest` to confirm the install works
(Python deps + jellyfish/samtools on PATH + the full count→EM→projection path)
*before* pointing kMate at their own data. It uses the small fixture shipped in
kmate/data/selftest/ (real Chr1 panel k-mers, a simulated pool with a KNOWN
founder mixture) and asserts EM recovers that mixture.

Exit code 0 = PASS, 1 = FAIL. Runs in a few seconds, no network, ~16 MB RAM.
"""
from __future__ import annotations
import argparse
import os
import sys
import tempfile

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "selftest")

# A correct install recovers the planted mixture near-perfectly; this threshold
# is deliberately loose so the check is a robust install signal, not a precision
# benchmark.
R2_PASS = 0.90


def _load_truth(founders):
    import numpy as np
    tmap = {}
    with open(os.path.join(DATA, "truth_h.tsv")) as f:
        next(f)
        for line in f:
            name, w = line.split()
            tmap[name] = float(w)
    h_true = np.array([tmap.get(str(x), 0.0) for x in founders])
    return h_true / h_true.sum(), tmap


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="kmate selftest",
        description="Run the bundled fixture end-to-end and verify the install.")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep-output", action="store_true",
                    help="keep the temporary output TSV / h npz instead of deleting")
    args = ap.parse_args(sys.argv[1:] if argv is None else argv)

    if not os.path.exists(os.path.join(DATA, "kmer_pa_Chr1.kmer_pa.npz")):
        sys.exit(f"FAIL: bundled fixture not found at {DATA} (package data missing)")

    import numpy as np
    from . import __version__, per_sample_per_chrom as driver

    print(f"kmate {__version__} selftest")
    print(f"  fixture: {DATA}")
    workdir = tempfile.mkdtemp(prefix="kmate_selftest_")
    out_tsv = os.path.join(workdir, "selftest.tsv")

    # Drive the real CLI path (same code users run), small hash for low memory.
    saved_argv = sys.argv
    sys.argv = [
        "kmate run",
        "--kmer-pa-prefix", os.path.join(DATA, "kmer_pa"),
        "--var-pa",     os.path.join(DATA, "var_Chr1.var_pa.npz"),
        "--var-called", os.path.join(DATA, "var_Chr1.var_called.npz"),
        "--var-meta",   os.path.join(DATA, "var_Chr1.meta.npz"),
        "--reads",      os.path.join(DATA, "reads.fq.gz"),
        "--sample", "selftest", "--out", out_tsv,
        "--threads", str(args.threads), "--chroms", "Chr1",
        "--block-mode", "global", "--hash-size", "100M",
    ]
    try:
        driver.main()
    except Exception as e:  # noqa: BLE001 — turn any failure into a clear FAIL
        sys.argv = saved_argv
        sys.exit(f"\nFAIL: pipeline raised {type(e).__name__}: {e}")
    sys.argv = saved_argv

    # --- checks --------------------------------------------------------------
    h_npz = out_tsv.replace(".tsv", ".h_per_chrom.npz")
    hd = np.load(h_npz, allow_pickle=True)
    founders = hd["founders"].astype(str)
    h_hat = hd["Chr1"].astype(float)
    h_hat = h_hat / h_hat.sum()
    h_true, tmap = _load_truth(founders)
    r2 = 1.0 - ((h_true - h_hat) ** 2).sum() / ((h_true - h_true.mean()) ** 2).sum()
    rmse = float(np.sqrt(np.mean((h_true - h_hat) ** 2)))

    af = np.genfromtxt(out_tsv, delimiter="\t", skip_header=1, usecols=4)
    af_finite = float(np.isfinite(af).mean())

    print(f"\n  recovered mixture (truth -> estimate):")
    for name, w in sorted(tmap.items(), key=lambda x: -x[1]):
        i = list(founders).index(name)
        print(f"    {name}: {w:.3f} -> {h_hat[i]:.3f}")
    print(f"\n  R^2={r2:.3f}  RMSE={rmse:.4f}  AF finite={af_finite:.0%}")

    ok = (r2 >= R2_PASS) and (af_finite > 0.99)
    if args.keep_output:
        print(f"  output kept in {workdir}")
    else:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)

    if ok:
        print("\nPASS — kMate is correctly installed and working.")
        sys.exit(0)
    print(f"\nFAIL — recovery below threshold (R^2>={R2_PASS}, AF finite>99%).")
    sys.exit(1)


if __name__ == "__main__":
    main()
