#!/usr/bin/env python
"""READ-ONLY. Prints the Ensembl GRCh37 release-87 GTF gene records relevant to the 12
unmapped CV targets (current symbols, HGNC previous symbols, and the competing records).

    python manuscript_prep/verify/gtf_lookup_unmapped.py
"""
import gzip
import re

GTF = r"D:\causal2drug\data\raw\Homo_sapiens.GRCh37.87.gtf.gz"
NAMES = {"ACE", "ITGB3", "TLR9", "ANGPTL8", "C19orf80", "CCN1", "CYR61", "GRK2", "ADRBK1",
         "GUCY1A1", "GUCY1A3", "GUCY1B1", "GUCY1B3", "VEGFD", "FIGF", "FIGNL2",
         "MIR29B1", "MIR92A1", "ADRA2B"}
IDS = {"ENSG00000159640", "ENSG00000264813", "ENSG00000259207", "ENSG00000259753",
       "ENSG00000239732", "ENSG00000173366"}
g = lambda k, s: re.search(k + ' "([^"]+)"', s).group(1)
print("gene_name\tgene_id\tchrom\tstart\tend\tstrand\tbiotype\tsource")
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        gn, gi = g("gene_name", f[8]), g("gene_id", f[8])
        if (gn in NAMES or gi in IDS or gn.startswith("MIR29B") or gn.startswith("MIR92A")
                or gn.startswith("FIGNL")):
            print("\t".join([gn, gi, f[0], f[3], f[4], f[6], g("gene_biotype", f[8]),
                             g("gene_source", f[8])]))
