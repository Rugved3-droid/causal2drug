#!/usr/bin/env python
"""Deliverable 1b: re-extract eQTLGen hits using the STABLE identifier join.

MAPPING METHOD (primary): Minikel symbol -> Ensembl gene ID via the pinned GRCh37
release-87 GTF -> eQTLGen `Gene` (ENSG) column. ENSG IDs are stable across Ensembl
versions; eQTLGen's `GeneSymbol` column is Ensembl v71 (~2013) nomenclature and drifts
from the newer HGNC symbols Minikel uses.

Nothing is aliased or guessed: a target that does not resolve to exactly one GRCh37 gene
record is UNMAPPED and is counted as such.
"""
import csv
import gzip
import time

SIG = (r"D:\causal2drug\data\raw"
       r"\2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
D = r"D:\causal2drug\deliverables\d1b"
P_INSTRUMENT = 5e-8


def main():
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    ensg2sym = {v["ensembl_id"]: k for k, v in coords.items()}
    print("target ENSGs (uniquely resolved in GRCh37 GTF):", len(ensg2sym))

    t0 = time.time()
    rows = []
    with gzip.open(SIG, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i["Gene"]] not in ensg2sym:
                continue
            if float(f[i["Pvalue"]]) >= P_INSTRUMENT:
                continue
            rows.append((ensg2sym[f[i["Gene"]]], f[i["Gene"]], f[i["SNP"]],
                         f[i["SNPChr"]], f[i["SNPPos"]], f[i["AssessedAllele"]],
                         f[i["OtherAllele"]], f[i["Zscore"]], f[i["Pvalue"]],
                         f[i["NrSamples"]], f[i["NrCohorts"]], f[i["GeneChr"]],
                         f[i["GenePos"]]))

    out = D + r"\eqtlgen_hits_by_ensg.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["MinikelSymbol", "Gene", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                    "OtherAllele", "Zscore", "Pvalue", "NrSamples", "NrCohorts",
                    "GeneChr", "GenePos"])
        w.writerows(rows)
    print("hit rows:", len(rows))
    print("distinct targets with >=1 hit (eQTLGen native +-1Mb window):",
          len({r[0] for r in rows}))
    print("seconds: %.1f" % (time.time() - t0))
    print("wrote", out)


if __name__ == "__main__":
    main()
