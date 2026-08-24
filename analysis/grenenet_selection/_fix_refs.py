#!/usr/bin/env python
"""Rewrite path references after the GEA restructure. Disposable.

Reads the authoritative old->new mapping from git's staged renames instead of
walking the filesystem (the walk stalls on the multi-TB output trees), and
rewrites only TRACKED text files plus a short explicit list of untracked ones.

  python _fix_refs.py           # dry run
  python _fix_refs.py --apply
"""
import os, subprocess, sys

GEA = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(GEA))
APPLY = "--apply" in sys.argv
EXT = (".py", ".sh", ".sbatch", ".md", ".ipynb", ".R", ".r")

def git(*a):
    return subprocess.run(["git", "-C", REPO, *a], capture_output=True,
                          text=True, check=True).stdout

# ---- old -> new, straight from the staged renames
pairs = []
for line in git("diff", "--cached", "--name-status", "-M").splitlines():
    p = line.split("\t")
    if p[0].startswith("R") and len(p) == 3:
        pairs.append((p[1], p[2]))

# Directory moves, stated EXPLICITLY. Deriving these by stripping the common
# suffix off file renames is unsafe: `gea/build_kendall.py` ->
# `gea/r2_gea_nonsnp/build_kendall.py` reduces to the prefix rule
# `analysis/grenenet_selection/` -> `analysis/grenenet_selection/r2_gea_nonsnp/`, which
# would rewrite every reference to the GEA root repo-wide.
B = "analysis/grenenet_selection/"
DIRS = {
    "phase1_replication":  "r2_gea_nonsnp",
    "wza_investigation":   "r2_gea_nonsnp",
    "cam5_replication":    "r2_gea_nonsnp",
    "driver_passenger":    "extras",
    "rerun_kfw_hb":        "common",
    "hap_blocks":          "blocks",
    "bigld_env":           "blocks",
    "seedmix_validation":  "qc",
}
dirmap = {f"{B}{d}/": f"{B}{dest}/{d}/" for d, dest in DIRS.items()}

# sanity: a rule must never rewrite the GEA root itself
for k in dirmap:
    assert k.rstrip("/") != B.rstrip("/"), f"unsafe rule {k}"

rules = sorted(set(list(dirmap.items()) + pairs), key=lambda kv: -len(kv[0]))
print(f"{len(pairs)} file renames + {len(dirmap)} explicit directory rules "
      f"-> {len(rules)} total")

# ---- candidate files: tracked + the few untracked top-level scripts.
# archive/ and old_docs/ are the historical record: their scripts were written
# against the old layout, so rewriting them would falsify it.
SKIP = ("analysis/grenenet_selection/archive/", "old_docs/", "archive/")
files = [f for f in git("ls-files").splitlines()
         if f.endswith(EXT) and not f.startswith(SKIP)]
files += [os.path.relpath(os.path.join(GEA, f), REPO)
          for f in os.listdir(GEA)
          if f.endswith(EXT) and os.path.isfile(os.path.join(GEA, f))]

hits, changed = 0, []
for rel in sorted(set(files)):
    p = os.path.join(REPO, rel)
    if not os.path.isfile(p):
        continue
    try:
        t = open(p, encoding="utf-8").read()
    except (UnicodeDecodeError, OSError):
        continue
    o = t
    for old, new in rules:
        if old in t:
            t = t.replace(old, new)
    if t != o:
        n = sum(1 for a, b in rules if a in o)
        hits += n
        changed.append(rel)
        if APPLY:
            open(p, "w", encoding="utf-8").write(t)

print(f"{'rewrote' if APPLY else 'would rewrite'} {len(changed)} files "
      f"({hits} rule hits)")
for c in changed[:40]:
    print("   ", c)
if len(changed) > 40:
    print(f"    ... +{len(changed)-40} more")
if not APPLY:
    print("\n(dry run -- pass --apply)")
