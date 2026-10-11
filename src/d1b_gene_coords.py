#!/usr/bin/env python
"""Deliverable 1b: resolve GRCh37 gene-body coordinates for every Minikel target.

Source is the pinned Ensembl GRCh37 release-87 GTF, verified against Ensembl's own
published CHECKSUMS file. (The Ensembl GRCh37 REST server returned HTTP 500 during this
run, and a pinned release file is more reproducible than a live API regardless.)

GRCh37 is used deliberately: eQTLGen phase I and the 1000G phase 3 LD panel are both
GRCh37, so the census needs no liftover. The GRCh37-vs-GRCh38 reconciliation required by
the full pipeline is a separate, still-unresolved pre-registration item.

MAPPING RULE: a Minikel symbol maps only if it matches exactly one `gene` feature on a
primary chromosome (1-22, X, Y, MT) in the GTF's `gene_name` attribute. Zero matches, or
more than one, are recorded as UNMAPPED with the reason. Nothing is guessed.
"""
import collections
import csv
import gzip
import re

FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
GTF = r"D:\causal2drug\data\raw\Homo_sapiens.GRCh37.87.gtf.gz"
OUTDIR = r"D:\causal2drug\deliverables\d1b"
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

RE_NAME = re.compile(r'gene_name "([^"]+)"')
RE_ID = re.compile(r'gene_id "([^"]+)"')
RE_BIO = re.compile(r'gene_biotype "([^"]+)"')


def main():
    with open(FROZEN, encoding="utf-8") as fh:
        targets = sorted({r["gene"] for r in csv.DictReader(fh, delimiter="\t")})
    tset = set(targets)
    print("Minikel unique targets:", len(targets))

    by_symbol = collections.defaultdict(list)
    n_gene_rows = 0
    with gzip.open(GTF, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            n_gene_rows += 1
            chrom = f[0]
            if chrom not in CHROMS:
                continue
            m = RE_NAME.search(f[8])
            if not m:
                continue
            sym = m.group(1)
            if sym not in tset:
                continue
            by_symbol[sym].append(dict(
                symbol=sym, ensembl_id=RE_ID.search(f[8]).group(1), chrom=chrom,
                start=int(f[3]), end=int(f[4]),
                biotype=(RE_BIO.search(f[8]).group(1) if RE_BIO.search(f[8]) else "")))
    print("gene features in GTF:", n_gene_rows)

    resolved, unmapped = {}, {}
    for s in targets:
        hits = by_symbol.get(s, [])
        if not hits:
            unmapped[s] = "not_in_GRCh37_87_GTF"
        elif len(hits) > 1:
            unmapped[s] = "ambiguous_%d_gene_records:%s" % (
                len(hits), ";".join("%s:%s:%d-%d" % (h["ensembl_id"], h["chrom"],
                                                     h["start"], h["end"])
                                    for h in hits))
        else:
            resolved[s] = hits[0]

    with open(OUTDIR + r"\gene_coords_grch37.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "ensembl_id", "chrom", "start", "end", "biotype"])
        for s in sorted(resolved):
            v = resolved[s]
            w.writerow([v["symbol"], v["ensembl_id"], v["chrom"], v["start"],
                        v["end"], v["biotype"]])
    with open(OUTDIR + r"\gene_coords_unmapped.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "reason"])
        for s in sorted(unmapped):
            w.writerow([s, unmapped[s]])

    print("resolved (exactly one gene record):", len(resolved))
    print("UNMAPPED:", len(unmapped))
    kinds = collections.Counter(v.split(":")[0] for v in unmapped.values())
    for k, n in kinds.most_common():
        print("   %-30s %d" % (k, n))
    print("\nfirst 15 unmapped:")
    for s in sorted(unmapped)[:15]:
        print("   %-12s %s" % (s, unmapped[s][:90]))


if __name__ == "__main__":
    main()
