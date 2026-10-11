#!/usr/bin/env python
"""ANALYSIS 2 -- every headline CV result under BOTH set definitions.

CV-PRIMARY  : Minikel therapeutic area == cardiovascular (MeSH C14)
CV-EXTENDED : CV-PRIMARY + indications matching a fixed lipid/atherosclerosis token list

The extension is the decision that admits PCSK9's indication, so it must be shown not to
be doing the work.
"""
import collections
import csv
import gzip
import json
import math
import re

import numpy as np
from scipy import stats

RAW = r"D:\causal2drug\data\raw"
CV = r"D:\causal2drug\deliverables\cv"
D1F = r"D:\causal2drug\deliverables\d1f"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")


def tsv(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def fisher_or(a, b, c, d):
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a+.5, b+.5, c+.5, d+.5)
    lor = math.log((aa*dd)/(bb*cc)); se = math.sqrt(1/aa+1/bb+1/cc+1/dd)
    return orr, math.exp(lor-1.96*se), math.exp(lor+1.96*se), p


def main():
    areas = {r["topl"]: r["area"] for r in tsv(RAW + r"\minikel_areas.tsv")}
    topl = {r["indication_mesh_id"]: r["topl"]
            for r in tsv(RAW + r"\minikel_indic_topl_match.tsv")}
    frozen = tsv(FROZEN)
    prim, ext = set(), set()
    pp, pe = 0, 0
    for r in frozen:
        isp = areas.get(topl.get(r["indication_mesh_id"])) == "cardiovascular"
        ise = isp or any(k in r["indication_mesh_term"].lower() for k in TOK)
        if isp:
            prim.add(r["gene"]); pp += 1
        if ise:
            ext.add(r["gene"]); pe += 1

    atlas = {r["gene"]: r for r in csv.DictReader(
        open(CV + r"\cardiovascular_atlas.csv", encoding="utf-8"))}

    def block(name, genes, npairs):
        rows = [atlas[g] for g in sorted(genes) if g in atlas]
        st = collections.Counter(r["state"] for r in rows)
        n = len(rows)
        out = {"definition": name, "n_pairs": npairs, "n_targets": n}
        for k in ("ABSENT", "TESTED_NULL", "INSTRUMENTED", "UNMAPPED"):
            out[k] = st[k]
            out[k + "_pct"] = 100.0 * st[k] / n
        out["CLASS_A"] = st["ABSENT"] + st["TESTED_NULL"]
        out["CLASS_A_pct"] = 100.0 * out["CLASS_A"] / n
        for lab, key in (("blood", "instrumented_blood"),
                         ("heart_LV", "instrumented_heart_LV"),
                         ("heart_atrial", "instrumented_heart_atrial"),
                         ("artery_aorta", "instrumented_artery_aorta"),
                         ("liver", "instrumented_liver"),
                         ("cardiac_union", "instrumented_cardiac_union")):
            c = sum(1 for r in rows if r[key] == "True")
            out[lab] = c
            out[lab + "_pct"] = 100.0 * c / n
        inh = [r for r in rows if r["tissue_restricted"] in ("YES", "NO")]
        a = sum(1 for r in inh if r["state"] == "ABSENT" and r["tissue_restricted"] == "YES")
        b = sum(1 for r in inh if r["state"] == "ABSENT" and r["tissue_restricted"] == "NO")
        c2 = sum(1 for r in inh if r["state"] == "INSTRUMENTED" and r["tissue_restricted"] == "YES")
        d = sum(1 for r in inh if r["state"] == "INSTRUMENTED" and r["tissue_restricted"] == "NO")
        orr, lo, hi, p = fisher_or(a, b, c2, d)
        out.update(or_=orr, or_lo=lo, or_hi=hi, or_p=p, or_table=[[a, b], [c2, d]])
        ca = [r for r in rows if r["untestability_class"] == "CLASS_A"]
        out["classA_launched"] = sum(1 for r in ca if r["max_phase_num"] == "5")
        out["classA_phase3"] = sum(1 for r in ca if r["max_phase_num"] == "4")
        out["classB"] = sum(1 for r in rows if r["untestability_class"] == "CLASS_B")
        mod = collections.Counter(r["required_modality"].split(" (")[0] for r in ca)
        out["modality"] = {k: v for k, v in mod.most_common()}
        out["modality_pct"] = {k: 100.0 * v / len(ca) for k, v in mod.items()}
        return out

    P = block("CV-PRIMARY", prim, pp)
    E = block("CV-EXTENDED", ext, pe)

    # ---- infrastructure attrition, both definitions, correct denominator ----
    harm = set()
    for l in open(RAW + r"\harmonised_list.txt"):
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m:
            harm.add(m.group(1))
    ANC = ["African", "East Asian", "South Asian", "Hispanic", "Latin", "Native",
           "Asian unspecified", "Oceanian", "Greater Middle Eastern", "Aboriginal",
           "Sub-Saharan", "NR"]
    eur = lambda s: ("European" in s) and not any(a in s for a in ANC)
    best = {}
    with open(RAW + r"\gwas-catalog-download-studies-v1.0.3.1.txt",
              encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd)
        it, ia, isz = (h.index("MAPPED_TRAIT"), h.index("STUDY ACCESSION"),
                       h.index("INITIAL SAMPLE SIZE"))
        for f in rd:
            if len(f) <= it or f[ia] not in harm or not eur(f[isz]):
                continue
            n = max([int(x.replace(",", "")) for x in re.findall(r"[\d,]{4,}", f[isz])] or [0])
            for t in f[it].split(","):
                t = t.strip().lower()
                if t and (t not in best or n > best[t][1]):
                    best[t] = (f[ia], n)
    skipped = {r["gwas_acc"] for r in tsv(D1F + r"\partc_skipped_accessions.tsv")}

    def attrition(pairset):
        # must be restricted to CV PAIRS, not to all indications of CV genes
        sel = [r for r in pairset
               if r["indication_mesh_term"].strip().lower() in best]
        lostn = sum(1 for r in sel
                    if best[r["indication_mesh_term"].strip().lower()][0] in skipped)
        return lostn, len(sel), (100.0 * lostn / len(sel) if sel else 0)

    cvp_pairs = [r for r in frozen
                 if areas.get(topl.get(r["indication_mesh_id"])) == "cardiovascular"]
    cve_pairs = [r for r in frozen
                 if areas.get(topl.get(r["indication_mesh_id"])) == "cardiovascular"
                 or any(k in r["indication_mesh_term"].lower() for k in TOK)]
    P["attr"] = attrition(cvp_pairs)
    E["attr"] = attrition(cve_pairs)

    hdr = "%-44s %-22s %-22s"
    print("=" * 92)
    print("ANALYSIS 2 -- every headline result under BOTH CV set definitions")
    print("=" * 92)
    print(hdr % ("", "CV-PRIMARY (C14 only)", "CV-EXTENDED (atlas)"))
    print("-" * 92)
    print(hdr % ("target-indication pairs", P["n_pairs"], E["n_pairs"]))
    print(hdr % ("distinct targets", P["n_targets"], E["n_targets"]))
    print("-" * 92)
    for k, lab in (("ABSENT", "ABSENT from eQTLGen"),
                   ("TESTED_NULL", "TESTED-BUT-NULL"),
                   ("CLASS_A", "CLASS A total"),
                   ("INSTRUMENTED", "instrumented in blood")):
        print(hdr % (lab, "%d (%.1f%%)" % (P[k], P[k+"_pct"]),
                     "%d (%.1f%%)" % (E[k], E[k+"_pct"])))
    print("-" * 92)
    for k, lab in (("blood", "yield: whole blood"), ("heart_LV", "yield: heart LV"),
                   ("heart_atrial", "yield: heart atrial"),
                   ("artery_aorta", "yield: artery aorta"), ("liver", "yield: liver"),
                   ("cardiac_union", "yield: cardiac union")):
        print(hdr % (lab, "%d (%.1f%%)" % (P[k], P[k+"_pct"]),
                     "%d (%.1f%%)" % (E[k], E[k+"_pct"])))
    print("-" * 92)
    print(hdr % ("Class A enrichment OR",
                 "%.2f [%.2f, %.2f]" % (P["or_"], P["or_lo"], P["or_hi"]),
                 "%.2f [%.2f, %.2f]" % (E["or_"], E["or_lo"], E["or_hi"])))
    print(hdr % ("  Fisher p", "%.3g" % P["or_p"], "%.3g" % E["or_p"]))
    print(hdr % ("  2x2 table", str(P["or_table"]), str(E["or_table"])))
    print("-" * 92)
    print(hdr % ("CLASS A with LAUNCHED drug", P["classA_launched"], E["classA_launched"]))
    print(hdr % ("CLASS A at Phase III", P["classA_phase3"], E["classA_phase3"]))
    print(hdr % ("CLASS B", P["classB"], E["classB"]))
    print("-" * 92)
    keys = sorted(set(P["modality"]) | set(E["modality"]),
                  key=lambda k: -E["modality"].get(k, 0))
    for k in keys:
        print(hdr % ("modality: " + k,
                     "%d (%.1f%%)" % (P["modality"].get(k, 0), P["modality_pct"].get(k, 0)),
                     "%d (%.1f%%)" % (E["modality"].get(k, 0), E["modality_pct"].get(k, 0))))
    print("-" * 92)
    print(hdr % ("infrastructure attrition (unindexed)",
                 "%d/%d = %.1f%%" % P["attr"], "%d/%d = %.1f%%" % E["attr"]))
    json.dump({"CV_PRIMARY": P, "CV_EXTENDED": E},
              open(CV + r"\cv_robustness.json", "w"), indent=2, default=str)
    print("\nwrote cv_robustness.json")


if __name__ == "__main__":
    main()
