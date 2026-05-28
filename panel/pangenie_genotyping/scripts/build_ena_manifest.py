"""
Query ENA for SRA/ENA Run accessions matching each of the 151 missing
GrENE-Net ecotype IDs.

Strategy: ecotype IDs are used as `sample_alias` in 1001 Genomes Project
submissions. We query ENA's data warehouse with `tax_id=3702 AND sample_alias=<id>`
and prefer PRJNA273563 (1001 Genomes original, paired-end Illumina, ~25× cov).

Output:
  data/ena_manifest.tsv — ecotype_id, run_accession, study_accession, fastq_ftp,
                          base_count, lib_strategy, lib_source
  data/ena_missing.tsv  — ecotypes for which no PRJNA273563 paired-end run found
"""
from __future__ import annotations
import sys, time, urllib.parse
from pathlib import Path
import urllib.request


BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = ["run_accession", "sample_alias", "study_accession",
          "fastq_ftp", "fastq_md5", "base_count", "read_count",
          "library_strategy", "library_source", "library_layout",
          "instrument_platform", "instrument_model"]


def query_ena(ecotype_id: str, retries: int = 3) -> list[dict]:
    """Return list of ENA records matching this ecotype as sample_alias."""
    q = f'tax_id=3702 AND sample_alias="{ecotype_id}"'
    params = {
        "result": "read_run",
        "query": q,
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": "100",
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                txt = r.read().decode()
            break
        except Exception as e:
            time.sleep(2)
    else:
        return []
    rows = []
    lines = txt.strip().split("\n")
    if not lines or len(lines) < 2:
        return []
    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        d = dict(zip(header, parts))
        rows.append(d)
    return rows


def pick_best(rows: list[dict]) -> dict | None:
    """Pick the best run for each ecotype — prefer PRJNA273563 paired-end Illumina."""
    if not rows:
        return None
    # Filter: Illumina, genomic, paired or single (jellyfish doesn't care about
    # pairing for k-mer counting; some PRJNA273563 runs are stored as a single
    # interleaved fastq with library_layout="PAIRED" but no ";" in fastq_ftp).
    candidates = [r for r in rows
                  if r.get("instrument_platform", "").upper() == "ILLUMINA"
                  and r.get("library_strategy", "").upper() in ("WGS", "WGA", "OTHER")
                  and r.get("library_source", "").upper() == "GENOMIC"
                  and r.get("fastq_ftp", "")]               # any fastq URL
    if not candidates:
        return None
    # Prefer 1001 Genomes (PRJNA273563); otherwise the run with highest base_count
    p1001g = [r for r in candidates if r.get("study_accession") == "PRJNA273563"]
    pool = p1001g if p1001g else candidates
    pool.sort(key=lambda r: int(r.get("base_count", 0) or 0), reverse=True)
    return pool[0]


def main():
    import argparse
    DATA = Path(__file__).parent.parent / "data"
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", default=str(DATA/"missing_151_ecotypes.txt"),
                    help="newline-delimited list of ecotype IDs")
    ap.add_argument("--out-manifest", default=str(DATA/"ena_manifest.tsv"))
    ap.add_argument("--out-missing", default=str(DATA/"ena_missing.tsv"))
    args = ap.parse_args()

    eco_list = [l.strip() for l in open(args.in_path) if l.strip()]
    print(f"Querying ENA for {len(eco_list)} ecotypes from {args.in_path}...")

    manifest = []
    missing = []
    for i, eid in enumerate(eco_list, 1):
        rows = query_ena(eid)
        best = pick_best(rows)
        if best is None:
            print(f"  [{i:3d}/{len(eco_list)}] {eid} → NO suitable run", flush=True)
            missing.append({"ecotype_id": eid, "n_total_rows": len(rows)})
        else:
            print(f"  [{i:3d}/{len(eco_list)}] {eid} → {best['run_accession']} "
                  f"({best['study_accession']}, {int(best['base_count'])/1e9:.1f} Gb)",
                  flush=True)
            manifest.append({
                "ecotype_id": eid,
                "run_accession": best["run_accession"],
                "study_accession": best["study_accession"],
                "fastq_ftp": best["fastq_ftp"],
                "fastq_md5": best.get("fastq_md5", ""),
                "base_count": best.get("base_count", "0"),
                "read_count": best.get("read_count", "0"),
                "library_strategy": best.get("library_strategy", ""),
                "library_layout": best.get("library_layout", ""),
                "instrument_model": best.get("instrument_model", ""),
            })

    cols = ["ecotype_id", "run_accession", "study_accession", "fastq_ftp",
            "fastq_md5", "base_count", "read_count", "library_strategy",
            "library_layout", "instrument_model"]
    with open(args.out_manifest, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in manifest:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nWrote manifest: {args.out_manifest}  ({len(manifest)} ecotypes)")

    if missing:
        with open(args.out_missing, "w") as f:
            f.write("ecotype_id\tn_total_rows\n")
            for r in missing:
                f.write(f"{r['ecotype_id']}\t{r['n_total_rows']}\n")
        print(f"Wrote missing list: {args.out_missing}  ({len(missing)} ecotypes)")

    # Summary stats
    if manifest:
        total_bytes = sum(int(r["base_count"]) for r in manifest)
        # base_count → roughly 1/3 the bytes when gzipped (lossy estimate)
        approx_gb = total_bytes / 1e9 / 3
        cov30_gb = sum((int(r["base_count"]) / 1e9) for r in manifest) / 3
        print(f"\nTotal: {sum(int(r['base_count']) for r in manifest)/1e12:.2f} Tb base count")
        print(f"Approx download size (paired-end gzipped): {approx_gb:.0f} GB")
        print(f"Studies represented:")
        from collections import Counter
        for s, c in Counter(r["study_accession"] for r in manifest).most_common():
            print(f"  {s}: {c} ecotypes")


if __name__ == "__main__":
    main()
