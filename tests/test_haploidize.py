"""
Check `kmate.haploidize` genotype rules on hand-written VCF lines (no data needed).

Tests:
  T1 — --het missing: diploid calls collapse, hets go missing
  T2 — --het missing: haploid calls (Minigraph-Cactus on haploid assemblies) pass through
  T3 — --het split: phased diploid becomes two haplotype columns
  T4 — --het split: haploid and unphased-het input are refused, not guessed

Run: python tests/test_haploidize.py
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kmate.haploidize import haploidize_line

FIXED = "Chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT"


def _counts():
    return {"cells": 0, "het_to_missing": 0, "other_to_missing": 0}


def _gts(line):
    return line.rstrip("\n").split("\t")[9:]


def test_missing_diploid():
    """T1"""
    c = _counts()
    out = haploidize_line(FIXED + "\t0/0\t1|1\t0|1\t./.\n", c, het="missing")
    assert _gts(out) == ["0", "1", ".", "."], _gts(out)
    assert c["het_to_missing"] == 1 and c["other_to_missing"] == 0, c
    print("  T1 ok")


def test_missing_haploid():
    """T2: a haploid panel must come through unchanged, not as missing."""
    c = _counts()
    out = haploidize_line(FIXED + "\t0\t1\t.\t1\n", c, het="missing")
    assert _gts(out) == ["0", "1", ".", "1"], _gts(out)
    assert c["other_to_missing"] == 0, c
    print("  T2 ok")


def test_split_phased():
    """T3"""
    c = _counts()
    out = haploidize_line(FIXED + "\t0|1\t1|1\t.|.\n", c, het="split")
    assert _gts(out) == ["0", "1", "1", "1", ".", "."], _gts(out)
    print("  T3 ok")


def test_split_refuses():
    """T4"""
    for bad in ("1", "0/1"):
        try:
            haploidize_line(FIXED + f"\t{bad}\n", _counts(), het="split")
        except ValueError:
            continue
        raise AssertionError(f"--het split accepted {bad!r}")
    print("  T4 ok")


if __name__ == "__main__":
    test_missing_diploid()
    test_missing_haploid()
    test_split_phased()
    test_split_refuses()
    print("PASS")
