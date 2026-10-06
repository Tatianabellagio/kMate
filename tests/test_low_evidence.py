"""
A fit with no evidence is reported as undefined, not computed by dividing by zero.

  T1 — solve_em on all-zero counts returns NaN and raises no RuntimeWarning
  T2 — with evidence, solve_em still returns a proper simplex point

Run: python tests/test_low_evidence.py
"""
from __future__ import annotations
import os, sys, warnings
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kmate.em_solver import solve_em

K = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]], dtype=np.float32)   # 3 founders x 4 k-mers


def test_no_counts():
    """T1"""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        h, info = solve_em(np.zeros(4, np.float32), K, 1.0)
    assert np.isnan(h).all() and info["iterations"] == 0 and not info["converged"], (h, info)
    print("  T1 ok")


def test_with_counts():
    """T2"""
    h, _ = solve_em(np.array([3, 5, 4, 2], np.float32), K, 1.0)
    assert np.isfinite(h).all() and abs(h.sum() - 1) < 1e-6 and (h >= 0).all(), h
    print("  T2 ok")


if __name__ == "__main__":
    test_no_counts()
    test_with_counts()
    print("PASS")
