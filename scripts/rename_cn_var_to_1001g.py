"""Rename cn_var_82 founder IDs from cactus 6-digit (e.g. '100042') to 1001G
IDs (e.g. '6911') using sample_rename.txt.

Cactus founders without a 1001G counterpart are dropped. Output cn_var has
a subset of cactus founders that DO map to 1001G IDs.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import load_npz, save_npz


def main():
    _root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument('--cn-var', default=str(_root / 'data/cn_var_82.cn_var.npz'))
    ap.add_argument('--cn-var-meta', default=str(_root / 'data/cn_var_82.meta.npz'))
    ap.add_argument('--rename-map', default=str(_root / 'panel/imputation/work/sample_rename.txt'))
    ap.add_argument('--out-prefix', default=str(_root / 'data/cn_var_82_renamed_to_1001g'))
    args = ap.parse_args()

    cn = load_npz(args.cn_var)        # F × N
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    cactus_ids = [str(f) for f in meta['founders']]
    print(f'cn_var_82: {cn.shape}, founders={len(cactus_ids)}')

    rename = pd.read_csv(args.rename_map, sep='\t', header=None,
                         names=['cactus', 'g1001'], dtype=str)
    rename_dict = dict(zip(rename['cactus'], rename['g1001']))
    print(f'rename map: {len(rename_dict)} cactus->1001G entries')

    # Pick which cactus founders we can rename
    keep_idx = []
    new_names = []
    for i, c in enumerate(cactus_ids):
        if c in rename_dict:
            keep_idx.append(i)
            new_names.append(rename_dict[c])
    print(f'kept {len(keep_idx)} cactus founders (dropping {len(cactus_ids) - len(keep_idx)} without 1001G mapping)')
    print(f'  dropped: {[c for c in cactus_ids if c not in rename_dict]}')

    cn_sub = cn[keep_idx, :].tocsr()
    save_npz(args.out_prefix + '.cn_var.npz', cn_sub)

    np.savez(args.out_prefix + '.meta.npz',
             founders=np.asarray(new_names, dtype=object),
             chrom=meta['chrom'],
             pos=meta['pos'],
             ref_len=meta['ref_len'],
             alt_len=meta['alt_len'])
    print(f'wrote {args.out_prefix}.cn_var.npz ({cn_sub.shape})')
    print(f'wrote {args.out_prefix}.meta.npz')


if __name__ == '__main__':
    main()
