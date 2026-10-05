#!/usr/bin/env python
"""v2 step 1 -- GEA unit pool for the NEW candidate list (events + call-rate 0.50, GWAS arm excluded).

Runs `build_gea_pool.py` unchanged on `tile_cr50` + the lambda table, but writes into
results/candidate_list_v2_events_cr50/ so nothing in the old / hybrid files under results/ is touched.
env: kmate. Compute node.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, CONV)
import build_gea_pool as m
m.OUT = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
m.main()
