"""Build the bundled `kmate selftest` fixture (run once; outputs are committed).

Produces a tiny, fully deterministic end-to-end fixture under
src/kmate/data/selftest/ so a freshly-installed user can run `kmate selftest`
and confirm their install works *before* pointing kMate at their own data.

It is derived from the real Chr1 panel slice (data/test_chr1_first200) so the
k-mers are genuine Arabidopsis 31-mers and the founder IDs are real — only the
pool is simulated, with a KNOWN founder mixture (truth_h.tsv) so the self-test
can assert that EM recovers it.

Outputs (all small enough to commit):
  kmer_pa_Chr1.kmer_pa.npz / .meta.npz   founder x k-mer membership (one chrom)
  var_Chr1.var_pa.npz / .var_called.npz / .meta.npz   founder x variant, for AF
  reads.fq.gz                             simulated pool reads (known mixture)
  truth_h.tsv                             founder, weight  (the planted mixture)

The reads are emitted so that each k-mer's pool count is ~ cov * (h_true @ A),
i.e. exactly the generative model EM inverts. Reproducible: fixed RNG seed.
"""
from __future__ import annotations
import gzip
import os
import numpy as np
from scipy.sparse import load_npz, csr_matrix, save_npz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC_PANEL = os.path.join(ROOT, "data", "test_chr1_first200")
OUT = os.path.join(ROOT, "src", "kmate", "data", "selftest")

# --- knobs (kept small so the fixture commits cleanly and runs in seconds) ----
N_BUBBLES = 60          # subsample of the 200 real bubbles -> ~few k k-mers
N_RECORDS = 400         # synthetic variant records for AF projection
COVERAGE = 5.0          # simulated pool depth — low on purpose: global (chrom-
                        # wide) EM pools k-mers across the chromosome, so it
                        # recovers the mixture cleanly even at ~5x. Keeps the
                        # bundled reads.fq.gz small too.
SEED = 0
# Planted mixture: 5 of the 82 founders carry the pool (rest are 0). The five
# carriers are chosen as the most-identifiable founders in the slice (most
# private/AC=1 k-mers) so EM can cleanly separate them — otherwise a founder
# whose k-mers are all shared is unrecoverable from this tiny panel and the
# self-test would flake. (Real production panels are genome-wide and don't
# have this small-slice degeneracy.)
MIX_WEIGHTS = np.array([0.35, 0.25, 0.20, 0.12, 0.08])


def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(OUT, exist_ok=True)

    # 1) Load the real panel slice and subsample to N_BUBBLES bubbles ----------
    A_full = load_npz(SRC_PANEL + ".kmer_pa.npz").tocsr()       # F x K (int8)
    meta = np.load(SRC_PANEL + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    F = len(founders)
    bubble_id_full = np.asarray(meta["bubble_id"]).astype(np.int64)

    keep_bubbles = np.arange(N_BUBBLES)                          # first N bubbles
    kmer_mask = np.isin(bubble_id_full, keep_bubbles)
    kcols = np.where(kmer_mask)[0]
    A = A_full[:, kcols].tocsc()
    K = A.shape[1]

    kmer_index = np.asarray(meta["kmer_index"])[kcols]
    bubble_id = bubble_id_full[kcols]
    # remap bubble ids to 0..n-1 (contiguous) and slice their coords
    uniq = np.unique(bubble_id)
    remap = {b: i for i, b in enumerate(uniq)}
    bubble_id = np.array([remap[b] for b in bubble_id], dtype=np.int64)
    b_chrom = np.asarray(meta["bubble_chrom"])[uniq]
    b_start = np.asarray(meta["bubble_start"])[uniq].astype(np.int64)
    b_end = np.asarray(meta["bubble_end"])[uniq].astype(np.int64)
    # force chrom label to Chr1 (the fixture is single-chrom)
    b_chrom = np.array(["Chr1"] * len(uniq))

    np.savez(os.path.join(OUT, "kmer_pa_Chr1.meta.npz"),
             kmer_index=kmer_index, bubble_id=bubble_id,
             bubble_chrom=b_chrom, bubble_start=b_start, bubble_end=b_end,
             founders=founders)
    save_npz(os.path.join(OUT, "kmer_pa_Chr1.kmer_pa.npz"), A.tocsr())
    print(f"kmer_pa: F={F} x K={K}  ({N_BUBBLES} bubbles)")

    # 2) Plant a known mixture over 5 founders, simulate the pool reads --------
    # pick the most-identifiable founders (most private/AC=1 k-mers in the slice)
    Aint = np.asarray(A.todense()).astype(np.int32)
    ac = Aint.sum(0)
    priv_per_f = Aint[:, ac == 1].sum(1)
    mix_founders = np.argsort(-priv_per_f)[:len(MIX_WEIGHTS)]
    h_true = np.zeros(F)
    h_true[mix_founders] = MIX_WEIGHTS
    h_true /= h_true.sum()

    # expected pool count per k-mer ~ cov * (h_true @ A); Poisson draw for realism
    Adense = np.asarray(A.todense()).astype(np.float64)          # F x K
    lam = COVERAGE * (h_true @ Adense)                          # length K
    counts = rng.poisson(lam).astype(int)

    # emit each k-mer string `counts[k]` times as a 31 bp read
    n_reads = 0
    with gzip.open(os.path.join(OUT, "reads.fq.gz"), "wt") as fq:
        for k in range(K):
            seq = str(kmer_index[k])
            qual = "I" * len(seq)
            for _ in range(int(counts[k])):
                fq.write(f"@r{n_reads}\n{seq}\n+\n{qual}\n")
                n_reads += 1
    print(f"reads: {n_reads:,} reads over {(counts>0).sum()} k-mers "
          f"(planted founders {sorted(mix_founders.tolist())})")

    with open(os.path.join(OUT, "truth_h.tsv"), "w") as f:
        f.write("founder\tweight\n")
        for i in np.argsort(-h_true):
            if h_true[i] > 0:
                f.write(f"{founders[i]}\t{h_true[i]:.6f}\n")

    # 3) Synthesize a small var_pa over the same founders for AF projection ----
    # Founder ALT carriage ~ Bernoulli(0.3); a few cells ./. (uncalled).
    carrier = (rng.random((F, N_RECORDS)) < 0.30)
    called = rng.random((F, N_RECORDS)) > 0.05                  # ~5% ./.
    carrier &= called                                          # ./. can't carry
    var_pa = csr_matrix(carrier.astype(np.int8))
    var_called = csr_matrix(called.astype(np.int8))
    save_npz(os.path.join(OUT, "var_Chr1.var_pa.npz"), var_pa)
    save_npz(os.path.join(OUT, "var_Chr1.var_called.npz"), var_called)

    pos = np.sort(rng.integers(int(b_start.min()), int(b_end.max()) + 1,
                               size=N_RECORDS)).astype(np.int64)
    ref = np.array(["A"] * N_RECORDS, dtype=object)
    alt = np.array(["T"] * N_RECORDS, dtype=object)
    np.savez(os.path.join(OUT, "var_Chr1.meta.npz"),
             founders=founders, chrom=np.array(["Chr1"] * N_RECORDS),
             pos=pos, ref=ref, alt=alt,
             ref_len=np.ones(N_RECORDS, dtype=np.int64),
             alt_len=np.ones(N_RECORDS, dtype=np.int64))
    print(f"var_pa: {F} x {N_RECORDS} records  (carrier nnz={var_pa.nnz:,})")
    print(f"\nFixture written to {OUT}")


if __name__ == "__main__":
    main()
