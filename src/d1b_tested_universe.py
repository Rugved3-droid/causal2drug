#!/usr/bin/env python
"""Deliverable 1b: extract the eQTLGen phase I TESTED-gene universe.

This is what separates the two failure modes the census must not merge:
  - "absent from eQTLGen entirely"  = symbol never appears as a tested gene
  - "tested but no instrument"      = tested, but no cis-eQTL at p < 5e-8

Requires the FULL results file, because the FDR<0.05 file only contains genes that
already have a significant eQTL and therefore cannot answer the first question.
"""
import csv
import gzip
import time

FULL = (r"D:\causal2drug\data\raw"
        r"\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
OUT = r"D:\causal2drug\deliverables\d1b\eqtlgen_tested_universe.tsv"


def main():
    t0 = time.time()
    genes = {}
    n = 0
    with gzip.open(FULL, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_gene, i_sym = hdr.index("Gene"), hdr.index("GeneSymbol")
        stop = max(i_gene, i_sym) + 1
        for line in fh:
            n += 1
            f = line.split("\t", stop)
            g = f[i_gene]
            if g not in genes:
                genes[g] = f[i_sym]
            if n % 50_000_000 == 0:
                print("  %d rows, %d genes, %.0fs" % (n, len(genes), time.time() - t0),
                      flush=True)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "GeneSymbol"])
        for g in sorted(genes):
            w.writerow([g, genes[g]])

    print("total rows in full file      :", n)
    print("distinct ENSG tested         :", len(genes))
    print("distinct GeneSymbol tested   :", len(set(genes.values())))
    print("seconds                      : %.1f" % (time.time() - t0))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
