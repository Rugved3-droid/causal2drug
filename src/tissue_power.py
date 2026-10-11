#!/usr/bin/env python
"""Measure empirical instrument power per candidate Lane B tissue.

For each eQTL Catalogue dataset: BH-FDR eGene rate from its own permutation output, and
the count of genes whose TOP nominal cis p-value clears p < 5e-8 -- the latter being the
quantity that actually matters, because DECISION 6 freezes the instrument threshold at
5e-8 regardless of a tissue's internal FDR calibration.
"""
import csv
import gzip
import sys

TISSUES = [
    ("QTD000266", "GTEx_v10 liver", 269),
    ("QTD000136", "GTEx_v10 artery (coronary)", 279),
    ("QTD000251", "GTEx_v10 heart (atrial appendage)", 473),
    ("QTD000256", "GTEx_v10 heart (left ventricle)", 481),
    ("QTD000131", "GTEx_v10 artery (aorta)", 492),
    ("QTD000141", "GTEx_v10 artery (tibial)", 736),
]
RAW = r"D:\causal2drug\data\raw\%s.permuted.tsv.gz"
CONTROLS = {"ENSG00000113161": "HMGCR", "ENSG00000169174": "PCSK9",
            "ENSG00000160712": "IL6R"}


def analyse(qtd):
    rows = list(csv.DictReader(gzip.open(RAW % qtd, "rt", encoding="utf-8"),
                              delimiter="\t"))
    vals = sorted((float(r["p_beta"]), r["molecular_trait_id"]) for r in rows)
    m = len(vals)
    thr = 0.0
    for i, (p, _g) in enumerate(vals, 1):
        if p <= 0.05 * i / m:
            thr = p
    n_egene = sum(1 for p, _g in vals if p <= thr)
    n_5e8 = sum(1 for r in rows if float(r["pvalue"]) < 5e-8)
    byid = {r["molecular_trait_id"]: r for r in rows}
    return m, n_egene, thr, n_5e8, byid


def main():
    print("%-34s %5s %7s %8s %7s %8s %7s"
          % ("tissue", "n", "genes", "eGenes", "eGene%", "top_p<5e-8", "%"))
    tables = {}
    for qtd, label, n in TISSUES:
        m, n_eg, thr, n_5e8, byid = analyse(qtd)
        tables[qtd] = (label, n, byid)
        print("%-34s %5d %7d %8d %6.1f%% %8d %6.1f%%"
              % (label, n, m, n_eg, 100.0 * n_eg / m, n_5e8, 100.0 * n_5e8 / m))
    print()
    print("Control genes: top nominal cis p-value in each tissue")
    print("%-34s %12s %12s %12s" % ("tissue", "HMGCR", "PCSK9", "IL6R"))
    for qtd, label, n in TISSUES:
        _l, _n, byid = tables[qtd]
        cells = []
        for ensg in ("ENSG00000113161", "ENSG00000169174", "ENSG00000160712"):
            r = byid.get(ensg)
            cells.append("%.3g" % float(r["pvalue"]) if r else "absent")
        print("%-34s %12s %12s %12s" % (label, cells[0], cells[1], cells[2]))


if __name__ == "__main__":
    main()
