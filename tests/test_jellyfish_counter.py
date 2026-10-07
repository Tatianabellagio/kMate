"""
build-index's genome-wide k-mer counts come from the `jellyfish` program alone
(no dna_jellyfish bindings, which bioconda ships for Python 3.9-3.12 only).

  T1 — counts match a direct count, for both orientations of each k-mer
  T2 — absent and N-containing k-mers count 0
  T3 — more queries than one pipe chunk come back complete and in order

Needs `jellyfish` (kmer-jellyfish) on PATH. Run: python tests/test_jellyfish_counter.py
"""
from __future__ import annotations
import os, random, shutil, sys, tempfile
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kmate.build_index import JellyfishCounter, canonical, revcomp

K = 11


def _setup():
    rng = random.Random(7)
    seqs = ["".join(rng.choice("ACGT") for _ in range(400)) for _ in range(3)]
    seqs.append(seqs[0][50:150])                     # repeat: these k-mers count 2
    d = tempfile.mkdtemp(prefix="jfc_")
    fa = os.path.join(d, "s.fa")
    with open(fa, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">s{i}\n{s}\n")
    truth = Counter(canonical(s[i:i + K].encode()) for s in seqs for i in range(len(s) - K + 1))
    return JellyfishCounter(fa, K, hash_size=100_000, threads=1, workdir=d), truth, d


def test_counts():
    jf, truth, d = _setup()
    kms = sorted(truth)
    try:
        assert jf.get_counts(kms) == [truth[k] for k in kms]                      # T1
        assert jf.get_counts([revcomp(k) for k in kms]) == [truth[k] for k in kms]
        print("  T1 ok")
        assert jf.get_counts([b"A" * K, b"ACGTNACGTAC"]) == [truth[b"A" * K], 0]  # T2
        print("  T2 ok")
        many = kms * (2 * JellyfishCounter._CHUNK // len(kms) + 2)               # T3
        assert len(many) > JellyfishCounter._CHUNK
        assert jf.get_counts(many) == [truth[k] for k in many]
        print("  T3 ok")
    finally:
        jf.close()
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    test_counts()
    print("PASS")
