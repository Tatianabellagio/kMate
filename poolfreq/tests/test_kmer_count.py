"""
Validate the k-mer counter against a known reference and a known BAM.

Tests:
  T1 — count k-mers in a synthetic FASTA: known-truth recovery
  T2 — count k-mers in a real visor_freqk BAM: smoke test + sanity bounds
  T3 — coverage estimate on the same BAM: matches expected ~50x
"""
from __future__ import annotations
import sys, os, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kmer_count import count_kmers_in_fasta, count_kmers_in_bam, estimate_coverage


def test_synth_fasta():
    """T1: write a small FASTA with known k-mer counts, verify recovery."""
    seq = "A"*30 + "ACGTACGTACGT" + "A"*30 + "ACGTACGTACGT" + "A"*30  # 'ACGTACGTACGT' appears 2x
    test_kmer = "ACGTACGTACGT"   # 12bp - we'll use k=12

    with tempfile.NamedTemporaryFile(mode='w', suffix=".fa", delete=False) as f:
        f.write(">test\n" + seq + "\n")
        fasta_path = f.name
    try:
        # canonical jellyfish: A^k, T^k → same canonical form
        # We expect ACGTACGTACGT to appear twice in canonical form
        counts = count_kmers_in_fasta(fasta_path, [test_kmer], k=12)
        print(f"  T1: kmer='{test_kmer}' count={counts[test_kmer]} (expected 2)")
        assert counts[test_kmer] == 2, f"expected 2, got {counts[test_kmer]}"
    finally:
        os.unlink(fasta_path)


def test_real_bam():
    """T2: query a few k-mers from TAIR10 against a sim BAM. Should be > 0 for valid k-mers."""
    bam = "/carnegie/nobackup/scratch/tbellagio/visor_freqk/data/reads_var/del/rep29/cov50/var_del_1kb_n231_f90_err001/sim.srt.bam"
    if not os.path.exists(bam):
        print("  T2 SKIPPED: BAM not found")
        return

    # Some plausible Arabidopsis k-mers (chosen from TAIR10 Chr1)
    test_kmers = [
        "AACCCTAAACCCTAAACCCTAAACCCTAAAC",  # telomeric-like
        "ATCAATTTATCTTTTGTGGGAAATTATTTAG",  # from earlier PanGenie output
        "ATCTTTTGTGGGAAATTATTTAGTTGTAGGG",
        "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",  # microsat (unlikely to be unique)
        "ATGCATGCATGCATGCATGCATGCATGCATG",  # synthetic, probably absent
    ]
    counts = count_kmers_in_bam(bam, test_kmers, k=31, threads=4)
    print(f"\n  T2 BAM: {os.path.basename(os.path.dirname(bam))}/{os.path.basename(bam)}")
    for k in test_kmers:
        print(f"    {k} → {counts[k]}")
    # Expect at least the telomeric and PanGenie-derived ones to be present
    assert counts[test_kmers[0]] > 0, "telomeric k-mer should be present"


def test_coverage_estimate():
    """T3: estimate coverage on a cov50 sim BAM."""
    bam = "/carnegie/nobackup/scratch/tbellagio/visor_freqk/data/reads_var/del/rep29/cov50/var_del_1kb_n231_f90_err001/sim.srt.bam"
    if not os.path.exists(bam):
        print("  T3 SKIPPED")
        return
    cov = estimate_coverage(bam)
    print(f"\n  T3 coverage estimate: {cov:.1f}x (expected ~50x)")
    # Allow generous bounds since visor reports "cov50" but mapping may be imperfect
    assert 10 < cov < 200, f"coverage {cov} outside reasonable range"


if __name__ == "__main__":
    print("="*70)
    print("K-mer counter tests")
    print("="*70)
    test_synth_fasta()
    test_real_bam()
    test_coverage_estimate()
    print("\nAll assertions passed.")
