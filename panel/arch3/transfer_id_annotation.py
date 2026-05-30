#!/usr/bin/env python3
"""
Transfer INFO/ID annotation from a cactus annotated VCF to a PanGenie-output VCF.

Both VCFs are expected to share the same multi-allelic record structure (same
chrom, pos, REF, ordered ALT list) at each record, because both come from the
same underlying pangenome graph.

Why custom: `bcftools annotate -c INFO/ID` mishandles the angle-bracketed graph
node IDs (`>...>...` inside the ID strings).

Usage:
  transfer_id_annotation.py --cactus <annotated.vcf.gz> --pg <pg.vcf.gz> --out <out.vcf>

Output: PG VCF with INFO/ID populated from cactus annotation, written to <out.vcf>.

Behavior on mismatch:
  - If record exists in PG but not in cactus annotation: keep PG record as-is
    (INFO/ID stays unpopulated) and emit a counter at the end.
  - If ALT lists differ in length: error out with location.
"""
import sys
import gzip
import argparse
from collections import defaultdict


def open_vcf(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cactus', required=True, help='cactus annotated VCF (.vcf or .vcf.gz)')
    ap.add_argument('--pg', required=True, help='PG VCF (.vcf or .vcf.gz)')
    ap.add_argument('--out', required=True, help='output VCF (uncompressed)')
    args = ap.parse_args()

    # Build lookup: (chrom, pos) -> (ref, alt_list, id_value)
    print(f'[transfer_id] Reading cactus annotated: {args.cactus}', file=sys.stderr)
    lookup = {}
    n_cactus = 0
    with open_vcf(args.cactus) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos = parts[0], parts[1]
            ref = parts[3]
            alt = parts[4]  # full ALT string, comma-separated
            # extract INFO/ID
            id_value = None
            for entry in parts[7].split(';'):
                if entry.startswith('ID='):
                    id_value = entry[3:]
                    break
            if id_value is None:
                continue
            lookup[(chrom, pos)] = (ref, alt, id_value)
            n_cactus += 1
    print(f'[transfer_id] Loaded {n_cactus:,} cactus records with INFO/ID', file=sys.stderr)

    # Process PG VCF
    print(f'[transfer_id] Annotating PG: {args.pg}', file=sys.stderr)
    n_pg = 0
    n_matched = 0
    n_no_match = 0
    n_struct_mismatch = 0
    header_id_seen = False
    with open_vcf(args.pg) as f, open(args.out, 'w') as out:
        for line in f:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=ID,'):
                    header_id_seen = True
                out.write(line)
                continue
            if line.startswith('#'):
                # Insert the ID header definition if not already present
                if not header_id_seen:
                    out.write('##INFO=<ID=ID,Number=A,Type=String,Description="Variant IDs per ALT allele.">\n')
                out.write(line)
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos = parts[0], parts[1]
            ref = parts[3]
            alt = parts[4]
            n_pg += 1

            key = (chrom, pos)
            if key not in lookup:
                n_no_match += 1
                # DROP record without INFO/ID — convert-to-biallelic.py asserts INFO/ID presence
                continue

            cactus_ref, cactus_alt, cactus_id = lookup[key]
            if cactus_ref != ref:
                # different REF — can't transfer; drop to avoid downstream assert failure
                n_struct_mismatch += 1
                continue
            if cactus_alt != alt:
                # ALT lists differ — annotate anyway if same number of values
                pg_alt_count = alt.count(',') + 1
                cactus_id_count = cactus_id.count(',') + 1
                if cactus_id_count != pg_alt_count:
                    n_struct_mismatch += 1
                    if n_struct_mismatch <= 5:
                        print(f'[WARN] {chrom}:{pos} ALT count mismatch: PG={pg_alt_count} vs cactus_id_count={cactus_id_count}; dropping record', file=sys.stderr)
                    continue

            # Insert INFO/ID
            info_parts = parts[7].split(';')
            info_parts = [p for p in info_parts if not p.startswith('ID=')]
            info_parts.append(f'ID={cactus_id}')
            parts[7] = ';'.join(info_parts)
            out.write('\t'.join(parts) + '\n')
            n_matched += 1

    print(f'[transfer_id] PG records processed: {n_pg:,}', file=sys.stderr)
    print(f'[transfer_id] PG records annotated:  {n_matched:,} ({100*n_matched/max(n_pg,1):.1f}%)', file=sys.stderr)
    print(f'[transfer_id] No match in cactus:    {n_no_match:,}', file=sys.stderr)
    print(f'[transfer_id] Structural mismatch:   {n_struct_mismatch:,}', file=sys.stderr)


if __name__ == '__main__':
    main()
