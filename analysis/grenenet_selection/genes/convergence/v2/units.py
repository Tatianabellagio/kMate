#!/usr/bin/env python
"""Shared unit helpers for the NEW (v2) candidate list: every evidence layer keys on the integer `unit` id.

A unit is a short-indel / MNP record or, for SVs, a Truvari EVENT (nonsnp_event_units.npz, built by
r1_sv_negative_selection/_build_event_units.py). `unit` is the row index in that file and the `unit` column of tile_cr50 / gea_pool.

  info(units)            DataFrame: unit, chrom, pos, ref_len, alt_len, dlen, n_members, n_carriers, n_called, rep_row, rep_col
  member_rows(unit)      af_store (index_nonsnp) rows of the unit's member records, in record order
  panel_cols(rows)       af_store rows -> (chrom, column index in that chrom's var_pa / meta arrays)
  rep_alleles(units)     REF / ALT strings of the representative member (meta.npz)
  carrier_matrix(units)  founders x units, True where the founder carries ANY member (union), plus the called matrix (any member)

The representative member = the one with most carriers (ties: first), as in _build_event_units.py. Never join a unit to anything on
(chrom, pos, ref_len, alt_len): 14% of non-SNP records share such a key.
env: kmate.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GEA = os.path.abspath(os.path.join(CONV, "..", ".."))
sys.path.insert(0, GEA)
import lib                                                      # noqa: E402

UNITS_NPZ = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
_cache: dict = {}


def _units():
    if "U" not in _cache:
        _cache["U"] = np.load(UNITS_NPZ, allow_pickle=True)
    return _cache["U"]


def _store_index():
    """af_store non-SNP index (chrom, pos, ref_len, alt_len per row) + per-chrom panel non-SNP column arrays."""
    if "S" not in _cache:
        idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
        ch = idx["chrom"].astype(str)
        off, cols, meta = {}, {}, {}
        o = 0
        for c in CHROMS:
            n = int((ch == c).sum())
            base = f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}"
            m = np.load(f"{base}.meta.npz", allow_pickle=True)
            prl = m["ref_len"].astype(np.int64); pal = m["alt_len"].astype(np.int64)
            panel = np.where(~((prl == 1) & (pal == 1)))[0]
            assert panel.size == n and np.array_equal(prl[panel], idx["ref_len"][o:o + n].astype(np.int64)), f"{c}: store/panel order differs"
            off[c] = o; cols[c] = panel; meta[c] = m
            o += n
        _cache["S"] = (idx, ch, off, cols, meta)
    return _cache["S"]


def panel_cols(rows):
    """af_store rows -> list of (chrom, panel column) pairs (vectorised per chrom)."""
    idx, ch, off, cols, _ = _store_index()
    rows = np.asarray(rows, np.int64)
    out_c = ch[rows]
    out_k = np.empty(len(rows), np.int64)
    for c in CHROMS:
        m = out_c == c
        if m.any():
            out_k[m] = cols[c][rows[m] - off[c]]
    return out_c, out_k


def member_rows(unit: int) -> np.ndarray:
    U = _units()
    return U["indices"][U["indptr"][unit]:U["indptr"][unit + 1]]


def info(units) -> pd.DataFrame:
    U = _units(); idx, ch, off, cols, _ = _store_index()
    units = np.asarray(units, np.int64)
    rep = U["rep"][units]
    d = pd.DataFrame({"unit": units, "chrom": U["chrom"][units].astype(str), "pos": U["pos"][units].astype(np.int64),
                      "ref_len": idx["ref_len"][rep].astype(np.int64), "alt_len": idx["alt_len"][rep].astype(np.int64),
                      "dlen": U["dlen"][units].astype(np.int64), "n_members": U["indptr"][units + 1] - U["indptr"][units],
                      "n_carriers": U["n"][units].astype(int), "n_called": U["ncal"][units].astype(int), "rep_row": rep})
    assert (d.chrom.to_numpy() == ch[rep]).all() and (d.pos.to_numpy() == idx["pos"][rep]).all(), "rep row disagrees with unit chrom/pos"
    d["rep_col"] = panel_cols(rep)[1]
    return d


def rep_alleles(units) -> pd.DataFrame:
    """REF / ALT strings of each unit's representative member record."""
    d = info(units); _, _, _, _, meta = _store_index()
    d["ref"] = ""; d["alt"] = ""
    for c in d.chrom.unique():                                   # NpzFile re-reads an array on every access: load once per chromosome
        R, A = meta[c]["ref"], meta[c]["alt"]
        m = d.chrom == c
        d.loc[m, "ref"] = [str(R[k]) for k in d.rep_col[m]]; d.loc[m, "alt"] = [str(A[k]) for k in d.rep_col[m]]
    assert (d.ref.str.len() == d.ref_len).all() and (d.alt.str.len() == d.alt_len).all(), "ref/alt strings disagree with stored lengths"
    return d


def carrier_matrix(units):
    """(founders, carriers[F x n_units] bool, called[F x n_units] bool): union over each unit's member records.

    Fast path: no full-matrix format conversion. For each chromosome the member columns of the requested units are mapped to a small
    index, the stored (row, col) entries of var_pa / var_called are filtered by that index in one vectorised pass, and member columns are
    folded into units with a sparse membership matrix.
    """
    units = list(map(int, units)); _, _, _, _, meta = _store_index()
    founders = [str(x) for x in meta["Chr1"]["founders"]]; NF = len(founders)
    car = np.zeros((NF, len(units)), bool); cal = np.zeros((NF, len(units)), bool)
    per_chrom: dict = {}
    for j, u in enumerate(units):
        c, k = panel_cols(member_rows(u))
        per_chrom.setdefault(str(c[0]), []).append((j, k))
    for c, items in per_chrom.items():
        base = f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}"
        N = len(meta[c]["pos"])
        need = np.unique(np.concatenate([k for _, k in items]))
        remap = np.full(N, -1, np.int64); remap[need] = np.arange(need.size)
        # membership: needed columns -> requested units of this chromosome
        mj, mc = [], []
        for j, k in items:
            mj += [j] * len(k); mc += list(remap[k])
        Mem = sp.csr_matrix((np.ones(len(mj), np.int32), (mc, mj)), shape=(need.size, len(units)))
        for kind, out in (("var_pa", car), ("var_called", cal)):
            z = np.load(f"{base}.{kind}.npz", allow_pickle=False)
            ip, ix = z["indptr"], z["indices"]
            hit = remap[ix] >= 0
            rows = np.repeat(np.arange(NF), np.diff(ip))[hit]
            X = sp.csr_matrix((np.ones(int(hit.sum()), np.int32), (rows, remap[ix][hit])), shape=(NF, need.size))
            out |= ((X @ Mem) > 0).toarray()
    return founders, car, cal
