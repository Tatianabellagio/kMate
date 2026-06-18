"""Unified `kmate` command-line interface.

Dispatches to the per-subcommand `main()`s. Each subcommand owns its own
argparse, so `kmate <cmd> --help` shows that command's full flag set
(including --block-mode / --kmer-weight on `run`, and the filter knobs).
"""
from __future__ import annotations
import sys
from . import __version__

# (subcommand -> (lazy import path, one-line help)). Lazy so `kmate --help`
# and `kmate --version` stay instant and don't import numpy/pysam.
_CMDS = {
    "run":           ("per_sample_per_chrom", "estimate founder freqs (h) + project per-record AF [the estimator]"),
    "build-kmer-pa": ("build_kmer_pa",        "build the founder x k-mer membership matrix (kmer_pa)"),
    "build-var-pa":  ("build_var_pa",         "build the founder x variant matrix (var_pa) from a panel VCF"),
    "build-kmer-db": ("build_kmer_db",        "build a Jellyfish k-mer DB from reads once (count-once path)"),
    "filter-pa":     ("filter_kmer_pa_production", "drop private (ac<=1) / invariant (ac==F) k-mer columns"),
    "selftest":      ("selftest",              "run the bundled fixture end-to-end to verify the install"),
}


def _usage() -> str:
    w = max(len(c) for c in _CMDS)
    lines = [f"  {c:<{w}}  {h}" for c, (_, h) in _CMDS.items()]
    return ("kmate — k-mer-based founder-frequency estimation for pool-seq\n\n"
            "usage: kmate <command> [options]\n\ncommands:\n" + "\n".join(lines) +
            "\n\nrun `kmate <command> --help` for a command's options.\n")


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv or argv[0] in ("-h", "--help"):
        sys.stdout.write(_usage())
        sys.exit(0 if argv else 1)
    if argv[0] in ("-V", "--version"):
        print(f"kmate {__version__}")
        sys.exit(0)
    cmd, rest = argv[0], argv[1:]
    if cmd not in _CMDS:
        sys.stderr.write(f"kmate: unknown command '{cmd}'\n\n" + _usage())
        sys.exit(2)
    module_name, _ = _CMDS[cmd]
    module = __import__(f"kmate.{module_name}", fromlist=["main"])
    # hand the remaining args to the subcommand's own argparse via sys.argv
    sys.argv = [f"kmate {cmd}"] + rest
    module.main()


if __name__ == "__main__":
    main()
