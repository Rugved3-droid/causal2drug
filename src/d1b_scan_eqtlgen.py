#!/usr/bin/env python
"""Deliverable 1b step 1-3: single streaming pass over the eQTLGen phase I
significant cis-eQTL file.

Writes:
  - eqtlgen_sig_gene_universe.tsv : every (ENSG, GeneSymbol) with >=1 FDR<0.05 cis-eQTL
  - eqtlgen_hits_targets.tsv      : every SNP-gene row with p < 5e-8 for a Minikel target
Also reports file-level facts (row count, max p, min p) used to verify that the
FDR<0.05 file is a superset of the p<5e-8 set.
"""
import collections
import csv
import gzip
import time

SIG = (r"D:\causal2drug\data\raw"
       r"\2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
OUT = r"D:\causal2drug\deliverables\d1b"

P_INSTRUMENT = 5e-8


def load_targets():
    with open(FROZEN, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return sorted({r["gene"] for r in rd})


def main():
    targets = load_targets()
    tset = set(targets)
    print("Minikel unique targets (frozen dedup set):", len(targets))

    t0 = time.time()
    n_rows = 0
    max_p = 0.0
    min_p = 1.0
    n_below = 0
    sig_genes = {}                       # ENSG -> symbol
    sym_to_ensg = collections.defaultdict(set)
    hits = []                            # rows for Minikel targets at p<5e-8

    with gzip.open(SIG, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_rows += 1
            p = float(f[i["Pvalue"]])
            if p > max_p:
                max_p = p
            if p < min_p:
                min_p = p
            ensg = f[i["Gene"]]
            sym = f[i["GeneSymbol"]]
            sig_genes[ensg] = sym
            sym_to_ensg[sym].add(ensg)
            if p < P_INSTRUMENT:
                n_below += 1
                if sym in tset:
                    hits.append((sym, ensg, f[i["SNP"]], f[i["SNPChr"]], f[i["SNPPos"]],
                                 f[i["AssessedAllele"]], f[i["OtherAllele"]],
                                 f[i["Zscore"]], f[i["Pvalue"]], f[i["NrSamples"]],
                                 f[i["NrCohorts"]], f[i["GeneChr"]], f[i["GenePos"]]))
    dt = time.time() - t0

    print("rows in significant file      :", n_rows)
    print("  min Pvalue                  : %.4g" % min_p)
    print("  max Pvalue (FDR<0.05 edge)  : %.4g" % max_p)
    print("  rows with p < 5e-8          :", n_below)
    print("  distinct ENSG in sig file   :", len(sig_genes))
    print("  distinct GeneSymbol in sig  :", len(set(sig_genes.values())))
    print("  scan seconds                : %.1f" % dt)

    with open(OUT + r"\eqtlgen_sig_gene_universe.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "GeneSymbol"])
        for e, s in sorted(sig_genes.items()):
            w.writerow([e, s])

    with open(OUT + r"\eqtlgen_hits_targets.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["GeneSymbol", "Gene", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                    "OtherAllele", "Zscore", "Pvalue", "NrSamples", "NrCohorts",
                    "GeneChr", "GenePos"])
        w.writerows(hits)
    print("hit rows for Minikel targets  :", len(hits))
    print("distinct target symbols hit   :", len({h[0] for h in hits}))

    # symbols mapping to more than one ENSG in eQTLGen -> ambiguous, counted unmapped
    amb = {s: v for s, v in sym_to_ensg.items() if len(v) > 1 and s in tset}
    print("ambiguous target symbols (>1 ENSG in eQTLGen):", len(amb))
    for s, v in sorted(amb.items()):
        print("   ", s, sorted(v))


if __name__ == "__main__":
    main()
