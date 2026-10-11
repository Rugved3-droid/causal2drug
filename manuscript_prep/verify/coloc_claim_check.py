#!/usr/bin/env python
"""Item 6 -- what the '96.9 % of cardiovascular target loci' figure is actually made of.

    python manuscript_prep/verify/coloc_claim_check.py

READ-ONLY. Source: deliverables/d1j/r15_results.tsv (106 screened loci, one per target,
each paired with the outcome GWAS Part C selected for one of that target's indications).
"""
import collections
import csv
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


ind = {r["indication_mesh_term"]: r for r in load(ROOT + r"\data\raw\minikel_indic.tsv")}
topl = {r["indication_mesh_id"]: r["topl"]
        for r in load(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
atlas = {r["gene"] for r in csv.DictReader(
    open(ROOT + r"\deliverables\cv\cardiovascular_atlas.csv", encoding="utf-8"))}
r15 = load(ROOT + r"\deliverables\d1j\r15_results.tsv")


def s0(t):
    return topl[ind[t]["indication_mesh_id"]] == "C14" or any(k in t.lower() for k in TOK)


def tagged(t):
    return "cardiovascular" in [a.strip() for a in ind[t]["areas"].split(",")]


def s0a(t):
    return s0(t) or (topl[ind[t]["indication_mesh_id"]] == "C23" and tagged(t))


def report(label, loci):
    sig = [r for r in loci if r["gwas_genome_wide_sig"] == "True"]
    use = [r for r in sig if r["susie_PP4"]]
    print("%-58s n=%3d  outcome GWS=%d  credible-set result=%d"
          % (label, len(loci), len(sig), len(use)))
    print("      outcomes: %s" % dict(collections.Counter(r["mesh_term"] for r in loci)
                                      .most_common()))
    if sig:
        print("      GWS loci: %s" % ", ".join("%s/%s" % (r["gene"], r["mesh_term"])
                                               for r in sig))


print("R15 screened loci: %d" % len(r15))
report("all screened loci", r15)
report("AS PUBLISHED: gene is in the 630 atlas (any outcome)", [r for r in r15
                                                               if r["gene"] in atlas])
report("outcome is an original-set (S0) CV indication", [r for r in r15 if s0(r["mesh_term"])])
report("outcome is a CV indication under proposed rule A (S0A)",
       [r for r in r15 if s0a(r["mesh_term"])])
report("outcome carries Minikel's cardiovascular tag at all",
       [r for r in r15 if s0(r["mesh_term"]) or tagged(r["mesh_term"])])
