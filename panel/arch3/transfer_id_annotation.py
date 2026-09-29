#!/usr/bin/env python3
"""Back-compat shim -> kmate.transfer_id (also reachable as `kmate transfer-id`).

Moved into the package so `kmate decompose` can use it without a repo checkout.
Kept so existing `python ../transfer_id_annotation.py ...` callers in chr*/jobA*.sh work.
"""
import os, sys
_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
sys.path.insert(0, os.path.abspath(_SRC))
from kmate.transfer_id import main

if __name__ == "__main__":
    main()
