#!/usr/bin/env python
"""Back-compat shim -> kmate.build_index. Prefer the installed CLI (`kmate build-index`).

The k-mer index builder moved into the package so that a conda/PyPI install can build a
panel without cloning the repo. Kept so existing
`python panel/pangenie_index/scripts/build_kmers_tsv.py ...` callers keep working.
"""
import os, sys
_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src")
sys.path.insert(0, os.path.abspath(_SRC))       # src/ -> 'import kmate' works
from kmate.build_index import main

if __name__ == "__main__":
    main()
