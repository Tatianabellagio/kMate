#!/usr/bin/env python
"""Back-compat shim -> kmate.filter_kmer_pa_production. Prefer the installed CLI (`kmate run`, etc.).

Kept so existing `python src/filter_kmer_pa_production.py ...` callers (grenenet/, panel/arch3/) work
without an install: adds src/ to sys.path and forwards to the package's main().
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # src/ -> 'import kmate' works
from kmate.filter_kmer_pa_production import main

if __name__ == "__main__":
    main()
