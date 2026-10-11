#!/usr/bin/env python
"""Item 2 -- resolve the gene-record mapping of the 12 UNMAPPED CV targets.

    python manuscript_prep/verify/resolve_unmapped.py

READ-ONLY on the project. Network: HGNC REST only (responses cached under
verify/cache/hgnc/). Writes manuscript_prep/unmapped_resolution.csv.

RESOLUTION RULE (fixed before looking at any eQTLGen state):
  1. Ask HGNC for the approved record of the Minikel symbol -> its Ensembl gene ID.
  2. If that Ensembl ID is a `gene` record in the pinned GRCh37 release-87 GTF, the
     target maps to it. This settles both failure types:
       - "ambiguous_2_gene_records": two GTF records share the symbol; keep the one
         HGNC names. The other is reported, and its eQTLGen status is shown too.
       - "not_in_GRCh37_87_GTF": the symbol is newer than the GTF; the HGNC Ensembl ID
         is in the GTF under the previous symbol.
  3. If HGNC gives no Ensembl ID, or the ID is not in the GTF, the target stays
     UNRESOLVED and is listed by name.

STATE, using the project's own definition (src/partb_classA.py:50-58), unchanged:
  ABSENT        Ensembl ID not among eQTLGen's tested genes
  INSTRUMENTED  >=1 eQTLGen cis-eQTL with p < 5e-8 inside gene body +-100 kb
  TESTED_NULL   tested, but no such variant
LD clumping was NOT re-run for these genes, so n_independent is not reported.
"""
import csv
import gzip
import json
import os
import re
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache", "hgnc")
os.makedirs(CACHE, exist_ok=True)
GTF = ROOT + r"\data\raw\Homo_sapiens.GRCh37.87.gtf.gz"
SIG = (ROOT + r"\data\raw\2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-"
       r"BonferroniAdded.txt.gz")
FLANK = 100_000
P = 5e-8


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def hgnc(sym):
    f = os.path.join(CACHE, sym + ".json")
    if not os.path.exists(f):
        req = urllib.request.Request("https://rest.genenames.org/fetch/symbol/" + sym,
                                     headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            open(f, "wb").write(r.read())
    docs = json.load(open(f, encoding="utf-8"))["response"]["docs"]
    return docs[0] if docs else None


unm = {r["symbol"]: r["reason"] for r in load(ROOT + r"\deliverables\d1b\gene_coords_unmapped.tsv")}
atlas = {r["gene"]: r for r in csv.DictReader(
    open(ROOT + r"\deliverables\cv\cardiovascular_atlas.csv", encoding="utf-8"))}
targets = sorted(g for g, r in atlas.items() if r["state"] == "UNMAPPED")
print("unmapped CV targets (%d): %s" % (len(targets), ", ".join(targets)))

gtf = {}
names = {}
g = lambda k, s: re.search(k + ' "([^"]+)"', s).group(1)
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        gid = g("gene_id", f[8])
        gtf[gid] = dict(name=g("gene_name", f[8]), chrom=f[0], start=int(f[3]),
                        end=int(f[4]), biotype=g("gene_biotype", f[8]),
                        source=g("gene_source", f[8]))
        names.setdefault(gtf[gid]["name"], []).append(gid)

tested = {r["Gene"]: r["GeneSymbol"] for r in load(
    ROOT + r"\deliverables\d1b\eqtlgen_tested_universe.tsv")}
tested_syms = set(tested.values())

res = []
for s in targets:
    h = hgnc(s)
    row = dict(gene=s, max_cv_phase=atlas[s]["max_clinical_phase"],
               original_failure=unm.get(s, "?"))
    if h is None:
        row.update(hgnc_id="", hgnc_ensembl_id="", resolution="UNRESOLVED: symbol not in HGNC")
        res.append(row)
        continue
    eid = h.get("ensembl_gene_id", "")
    prev = h.get("prev_symbol", [])
    row.update(hgnc_id=h.get("hgnc_id", ""), hgnc_locus_type=h.get("locus_type", ""),
               hgnc_prev_symbols=";".join(prev), hgnc_ensembl_id=eid)
    if not eid:
        row["resolution"] = "UNRESOLVED: HGNC lists no Ensembl gene ID"
    elif eid not in gtf:
        alt = [x for p_ in prev + [s] for x in names.get(p_, [])]
        row["resolution"] = ("UNRESOLVED: HGNC Ensembl ID not a gene record in GRCh37.87"
                             + (" (GTF has %s under a matching symbol; NOT used, rule 2 "
                                "requires the HGNC ID)" % ",".join(alt) if alt else ""))
    else:
        r_ = gtf[eid]
        others = [x for x in names.get(s, []) + [y for p_ in prev for y in names.get(p_, [])]
                  if x != eid]
        row.update(resolved_ensembl_id=eid, gtf_gene_name=r_["name"], chrom=r_["chrom"],
                   start=r_["start"], end=r_["end"], biotype=r_["biotype"],
                   gtf_source=r_["source"],
                   other_gtf_records_same_symbol=";".join(sorted(set(others))),
                   resolution=("RESOLVED: HGNC Ensembl ID selects one of two GTF records"
                               if "ambiguous" in unm.get(s, "") else
                               "RESOLVED: GTF carries the gene under previous symbol %s"
                               % r_["name"]))
    res.append(row)

# eQTLGen lookup for resolved IDs (and the competing records, for transparency)
want = {}
for r in res:
    if r.get("resolved_ensembl_id"):
        want[r["resolved_ensembl_id"]] = r
        for o in filter(None, r["other_gtf_records_same_symbol"].split(";")):
            want.setdefault(o, None)
hits = {k: [0, 1.0, ""] for k in want}       # n in window, best p, best snp
any_sig = {k: 0 for k in want}
with gzip.open(SIG, "rt", encoding="utf-8") as fh:
    hd = fh.readline().rstrip("\n").split("\t")
    ix = {c: hd.index(c) for c in ("Pvalue", "SNP", "SNPChr", "SNPPos", "Gene")}
    for line in fh:
        f = line.split("\t")
        gid = f[ix["Gene"]]
        if gid not in want:
            continue
        any_sig[gid] += 1
        p = float(f[ix["Pvalue"]])
        c = gtf[gid]
        if (p < P and f[ix["SNPChr"]] == c["chrom"]
                and c["start"] - FLANK <= int(f[ix["SNPPos"]]) <= c["end"] + FLANK):
            hits[gid][0] += 1
            if p < hits[gid][1]:
                hits[gid][1], hits[gid][2] = p, f[ix["SNP"]]

for r in res:
    eid = r.get("resolved_ensembl_id")
    if not eid:
        r["state_after_resolution"] = "UNMAPPED"
        continue
    r["in_eqtlgen_tested_genes"] = eid in tested
    r["symbol_or_prev_symbol_in_eqtlgen"] = any(
        x in tested_syms for x in [r["gene"], r["gtf_gene_name"]]
        + r["hgnc_prev_symbols"].split(";"))
    r["n_fdr05_ciseqtl_rows"] = any_sig[eid]
    r["n_hits_p5e-8_within_100kb"] = hits[eid][0]
    r["best_p_in_window"] = hits[eid][1] if hits[eid][0] else ""
    r["best_snp_in_window"] = hits[eid][2]
    r["state_after_resolution"] = ("ABSENT" if eid not in tested else
                                   "INSTRUMENTED" if hits[eid][0] else "TESTED_NULL")
    oth = []
    for o in filter(None, r["other_gtf_records_same_symbol"].split(";")):
        oth.append("%s:%s,hits=%d" % (o, "tested" if o in tested else "not_tested",
                                      hits[o][0]))
    r["other_record_eqtlgen_status"] = ";".join(oth)

cols = ["gene", "max_cv_phase", "original_failure", "hgnc_id", "hgnc_locus_type",
        "hgnc_prev_symbols", "hgnc_ensembl_id", "resolution", "resolved_ensembl_id",
        "gtf_gene_name", "chrom", "start", "end", "biotype", "gtf_source",
        "other_gtf_records_same_symbol", "other_record_eqtlgen_status",
        "in_eqtlgen_tested_genes", "symbol_or_prev_symbol_in_eqtlgen",
        "n_fdr05_ciseqtl_rows", "n_hits_p5e-8_within_100kb", "best_p_in_window",
        "best_snp_in_window", "state_after_resolution"]
out = os.path.join(os.path.dirname(HERE), "unmapped_resolution.csv")
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(res)
for r in res:
    print("%-8s %-9s -> %-16s %-9s %-13s hits=%s best_p=%s | %s | other: %s"
          % (r["gene"], r["max_cv_phase"], r.get("resolved_ensembl_id", "-"),
             r.get("gtf_gene_name", "-"), r["state_after_resolution"],
             r.get("n_hits_p5e-8_within_100kb", "-"), r.get("best_p_in_window", "-"),
             r["resolution"], r.get("other_record_eqtlgen_status", "")))
import collections
print(dict(collections.Counter(r["state_after_resolution"] for r in res)))
print("wrote", out)
