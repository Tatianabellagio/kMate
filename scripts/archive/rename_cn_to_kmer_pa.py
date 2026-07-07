#!/usr/bin/env python3
"""
One-shot project-wide rename of the legacy `cn` matrix nomenclature to the
presence/absence names (K_pa / V_pa / V_called).

  cn_full        -> kmer_pa      (founder × k-mer presence/absence; dirs, files, vars)
  cn_var_called  -> var_called   (call mask)
  cn_var         -> var_pa       (founder × variant alt-allele presence)
  cn (bare)      -> kmer_pa       (the EM-evidence matrix variable / .cn.npz suffix)
  --cn-* flags   -> --kmer-pa-* / --var-pa / --var-called / --var-meta

Renames BOTH file contents and file/dir names, and the on-disk data artifacts
under data/ and panel/arch3/ so the new code finds them.

Usage:
  python scripts/rename_cn_to_kmer_pa.py --dry-run   # report only
  python scripts/rename_cn_to_kmer_pa.py --apply     # do it

Frozen areas (never touched): .git, archive/, old_docs/, */archive/,
panel/imputation/ (deprecated), data/block_haplotype_cn, external/, papers/,
results/ plots/ logs/, and this script itself.
"""
import os
import re
import sys
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Ordered (longest / most-specific first). Each is (pattern, repl, is_regex).
RULES = [
    ("cn_var_called", "var_called", False),
    ("cn_var_meta",   "var_meta",   False),
    ("cn_var",        "var_pa",     False),
    ("cn_full",       "kmer_pa",    False),
    ("cn_production", "kmer_pa_production", False),
    ("cn_kmer",       "kmer_pa",    False),
    ("kmer_cn",       "kmer_pa",    False),
    ("cn_dense",      "kmer_pa_dense",  False),
    ("cn_em",         "kmer_pa_em",     False),
    ("cn_prefix",     "kmer_pa_prefix", False),
    ("cn_for_chrom",  "kmer_pa_for_chrom", False),
    ("cn_Chr",        "kmer_pa_Chr",  False),
    ("CN_VAR_CALLED", "VAR_CALLED",  False),
    # hyphenated CLI flags
    ("cn-var-called", "var-called",  False),
    ("cn-var-meta",   "var-meta",    False),
    ("cn-var",        "var-pa",      False),
    ("cn-kmer",       "kmer-pa",     False),
    # bare token last
    (r"\bcn\b",       "kmer_pa",     True),
]


def apply_tokens(s):
    for pat, repl, is_re in RULES:
        if is_re:
            s = re.sub(pat, repl, s)
        else:
            s = s.replace(pat, repl)
    return s


# Directory names (any path component) that freeze the whole subtree.
EXCLUDE_DIRS = {
    ".git", "archive", "old_docs", "imputation", "external", "papers",
    "results", "plots", "logs", "__pycache__", ".ipynb_checkpoints",
    "block_haplotype_cn",
    # virtualenvs / installed packages — NEVER touch
    ".venv", "venv", "env", "site-packages", "node_modules", ".tox", ".mypy_cache",
    # heavy working / output areas — deferred (rerun later if wanted)
    "scratch", "sims", "contamination_test", "preprocess_qc", "notebooks",
    "pangenie_index", "pangenie_genotyping",
}
# Specific excluded subpaths (relative to ROOT).
EXCLUDE_PATHS = {"src/archive", "scripts/archive", "panel/imputation"}

TEXT_EXT = {".py", ".sh", ".sbatch", ".md", ".txt", ".json", ".cfg",
            ".yaml", ".yml", ".tsv", ".csv"}
# Extensions we rename but never edit contents of.
BINARY_RENAME_EXT = {".npz", ".gz", ".vcf", ".bgz", ".tbi", ".csi"}

THIS_FILE = os.path.abspath(__file__)


def is_excluded(path):
    rel = os.path.relpath(path, ROOT)
    parts = rel.split(os.sep)
    if any(p in EXCLUDE_DIRS for p in parts):
        return True
    for ex in EXCLUDE_PATHS:
        if rel == ex or rel.startswith(ex + os.sep):
            return True
    return False


def git_tracked(path):
    r = subprocess.run(["git", "ls-files", "--error-unmatch", path],
                       cwd=ROOT, capture_output=True)
    return r.returncode == 0


def do_rename(src, dst, apply):
    print(f"  RENAME {os.path.relpath(src, ROOT)}  ->  {os.path.relpath(dst, ROOT)}")
    if not apply:
        return
    if git_tracked(src):
        subprocess.run(["git", "mv", src, dst], cwd=ROOT, check=True)
    else:
        os.rename(src, dst)


def main():
    apply = "--apply" in sys.argv
    if not apply and "--dry-run" not in sys.argv:
        print("specify --dry-run or --apply"); sys.exit(2)

    content_changed = 0
    renamed = 0

    # --- Pass 1: file contents ---
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames
                       if not is_excluded(os.path.join(dirpath, d))]
        for fn in filenames:
            fp = os.path.join(dirpath, fn)
            if os.path.abspath(fp) == THIS_FILE or is_excluded(fp):
                continue
            if os.path.splitext(fn)[1] not in TEXT_EXT:
                continue
            try:
                with open(fp, encoding="utf-8") as fh:
                    old = fh.read()
            except (UnicodeDecodeError, OSError):
                continue
            new = apply_tokens(old)
            if new != old:
                content_changed += 1
                print(f"  EDIT   {os.path.relpath(fp, ROOT)}")
                if apply:
                    with open(fp, "w", encoding="utf-8") as fh:
                        fh.write(new)

    # --- Pass 2: file + dir names ---
    # Collect rename ops during a PRUNED topdown walk (so we never recurse into
    # .venv/.git/scratch), then apply deepest-first so children rename before parents.
    ops = []  # (src, dst)
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames
                       if not is_excluded(os.path.join(dirpath, d))]
        for fn in filenames:
            fp = os.path.join(dirpath, fn)
            if os.path.abspath(fp) == THIS_FILE or is_excluded(fp):
                continue
            ext = os.path.splitext(fn)[1]
            if ext not in TEXT_EXT and ext not in BINARY_RENAME_EXT:
                continue
            new_fn = apply_tokens(fn)
            if new_fn != fn:
                ops.append((fp, os.path.join(dirpath, new_fn)))
        if dirpath != ROOT:
            base = os.path.basename(dirpath)
            new_base = apply_tokens(base)
            if new_base != base:
                ops.append((dirpath, os.path.join(os.path.dirname(dirpath), new_base)))

    # deepest paths first (longest string = deepest), so a dir is renamed only
    # after all its contents have been renamed.
    ops.sort(key=lambda p: p[0].count(os.sep), reverse=True)
    for src, dst in ops:
        renamed += 1
        do_rename(src, dst, apply)

    print(f"\n{'APPLIED' if apply else 'DRY-RUN'}: "
          f"{content_changed} files edited, {renamed} paths renamed.")


if __name__ == "__main__":
    main()
