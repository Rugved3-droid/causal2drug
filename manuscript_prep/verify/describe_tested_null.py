#!/usr/bin/env python
"""Item 3 -- what TESTED_NULL contains. Descriptive only; reclassifies nothing.

    python manuscript_prep/verify/describe_tested_null.py

READ-ONLY. Uses files already produced by Deliverable 1b:
  d1b/eqtlgen_hits_by_ensg.tsv        every eQTLGen row with p < 5e-8 for a Minikel target,
                                      anywhere in eQTLGen's own cis window
  d1b/eqtlgen_sig_gene_universe.tsv   every gene with >=1 FDR < 0.05 cis-eQTL
  d1b/eqtlgen_tested_universe.tsv     every gene eQTLGen tested
Writes manuscript_prep/tested_null_detail.csv (the 55 original-set CV targets).
"""
import collections
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
D1B = ROOT + r"\deliverables\d1b"
HERE = os.path.dirname(os.path.abspath(__file__))
FLANK = 100_000


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


coords = {r["symbol"]: r for r in load(D1B + r"\gene_coords_grch37.tsv")}
tested = {r["Gene"] for r in load(D1B + r"\eqtlgen_tested_universe.tsv")}
sig = load(D1B + r"\eqtlgen_sig_gene_universe.tsv")
print("sig_gene_universe columns:", list(sig[0].keys()))
sig_ensg = {r[list(r.keys())[0]] for r in sig}
ann = {r["symbol"]: r for r in load(ROOT + r"\deliverables\d1e\classA_annotation.tsv")}
atlas = {r["gene"]: r for r in csv.DictReader(
    open(ROOT + r"\deliverables\cv\cardiovascular_atlas.csv", encoding="utf-8"))}
clump = {r["symbol"]: r for r in load(D1B + r"\clumped_instruments_final.tsv")}

inside = collections.Counter()
outside = collections.defaultdict(list)       # (p, distance to gene body, snp)
with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    print("hits_by_ensg columns:", rd.fieldnames)
    for r in rd:
        s = r["MinikelSymbol"]
        c = coords.get(s)
        if not c:
            continue
        if r["SNPChr"] != c["chrom"]:
            outside[s].append((float(r["Pvalue"]), None, r["SNP"]))
            continue
        p = int(r["SNPPos"])
        a, b = int(c["start"]), int(c["end"])
        if a - FLANK <= p <= b + FLANK:
            inside[s] += 1
        else:
            outside[s].append((float(r["Pvalue"]), (a - p) if p < a else (p - b), r["SNP"]))

# consistency of the state rule with the shipped annotation
bad = 0
for s, r in ann.items():
    if s not in coords:
        exp = "UNMAPPED"
    elif coords[s]["ensembl_id"] not in tested:
        exp = "ABSENT"
    elif inside[s]:
        exp = "INSTRUMENTED"
    else:
        exp = "TESTED_NULL"
    bad += exp != r["state"]
print("state rule re-applied to all %d targets; disagreements with classA_annotation.tsv: %d"
      % (len(ann), bad))
print("INSTRUMENTED targets: %d ; in clump file: %d ; with n_independent == 0: %d"
      % (sum(1 for r in ann.values() if r["state"] == "INSTRUMENTED"), len(clump),
         sum(1 for r in clump.values() if r["n_independent"] == "0")))


def sub(s):
    e = coords[s]["ensembl_id"]
    if outside[s]:
        return "T1 genome-wide-significant cis-eQTL exists, but only outside gene body +-100 kb"
    if e in sig_ensg:
        return "T2 cis-eQTL at FDR<0.05 exists, none reaches p<5e-8"
    return "T3 tested, no cis-eQTL at FDR<0.05"


for label, pool in (("ALL 2,517 targets", [s for s, r in ann.items()
                                          if r["state"] == "TESTED_NULL"]),
                    ("original 630-target CV set", [g for g, r in atlas.items()
                                                    if r["state"] == "TESTED_NULL"])):
    c = collections.Counter(sub(s) for s in pool)
    print("\nTESTED_NULL in %s: %d" % (label, len(pool)))
    for k in sorted(c):
        print("   %3d  (%.1f%%)  %s" % (c[k], 100 * c[k] / len(pool), k))

rows = []
for g, r in sorted(atlas.items()):
    if r["state"] != "TESTED_NULL":
        continue
    o = sorted(outside[g])
    rows.append(dict(
        gene=g, ensembl_id=r["ensembl_id"], max_cv_phase=r["max_clinical_phase"],
        subtype=sub(g), gtex_blood_tpm=r["gtex_blood_tpm"],
        tissue_restricted=r["tissue_restricted"],
        n_hits_outside_100kb=len(o),
        best_p_outside=o[0][0] if o else "",
        best_snp_outside=o[0][2] if o else "",
        min_distance_bp_from_gene_body=min((x[1] for x in o if x[1] is not None),
                                           default=""),
        gene_length_bp=int(coords[g]["end"]) - int(coords[g]["start"])))
out = os.path.join(os.path.dirname(HERE), "tested_null_detail.csv")
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("\nwrote", out)
tp = [float(r["gtex_blood_tpm"]) for r in rows if r["gtex_blood_tpm"]]
print("CV tested-null with GTEx blood TPM: %d ; >=1 TPM: %d ; >=10 TPM: %d"
      % (len(tp), sum(1 for x in tp if x >= 1), sum(1 for x in tp if x >= 10)))
print("launched / Phase III tested-null:")
for r in rows:
    if r["max_cv_phase"] in ("Launched", "Phase III"):
        print("   %-7s %-9s blood TPM %-8s %s | outside hits %s best p %s nearest %s bp"
              % (r["gene"], r["max_cv_phase"], r["gtex_blood_tpm"], r["subtype"][:2],
                 r["n_hits_outside_100kb"], r["best_p_outside"],
                 r["min_distance_bp_from_gene_body"]))
