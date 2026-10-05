"""
Retired options must stay gone: `kmate run` rejects them rather than accepting
and ignoring them.

  T1 — --kmer-weight (the per-bubble 1/m_b weight, removed 2026-09-29), in full and
       abbreviated, is an unrecognized argument

Run: python tests/test_cli_flags.py
"""
from __future__ import annotations
import os, subprocess, sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
BASE = ["--kmer-pa-prefix", "x", "--var-pa-prefix", "v", "--reads", "r.fq",
        "--sample", "s", "--out", "o.tsv"]


def _run(*extra):
    env = dict(os.environ, PYTHONPATH=SRC)
    return subprocess.run([sys.executable, "-m", "kmate.cli", "run", *BASE, *extra],
                          capture_output=True, text=True, env=env)


def test_kmer_weight_rejected():
    """T1"""
    for flag in ("--kmer-weight", "--kmer-w"):
        p = _run(flag, "inv_mb")
        assert p.returncode != 0 and "unrecognized arguments" in p.stderr, (flag, p.stderr)
    print("  T1 ok")


if __name__ == "__main__":
    test_kmer_weight_rejected()
    print("PASS")
