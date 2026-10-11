#!/usr/bin/env python
"""Counts behind two statements in the manuscript, read from the classification and the
recovery map (nothing is recomputed or reclassified).

    python manuscript_prep/manuscript/recovery_scripts/lead_findings.py \
        > manuscript_prep/manuscript/recovery_out/lead_findings_output.txt

1. The lipid targets whose human genetics is best established, and their status in eQTLGen.
2. Locked heart failure / cardiomyopathy targets, separating genes that a plasma protein panel
   could never measure (microRNA genes, mitochondrial genes) from nuclear protein-coding genes.
"""
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
P = {r["gene"]: r for r in csv.DictReader(
    open(os.path.join(MP, "analysis", "out", "per_target_primary.csv"), encoding="utf-8"))}
M = {r["gene"]: r for r in csv.DictReader(open(os.path.join(MS, "recovery_map.csv"),
                                               encoding="utf-8"))}
LIPID = ["PCSK9", "APOC3", "ANGPTL3", "APOB", "NPC1L1", "LPA"]
n_out = 0
for g in LIPID:
    p = P[g]
    n_out += p["status"] != "INSTRUMENTED"
    print("[L01] %-8s eQTLGen: %-25s phase: %-9s whole-blood TPM: %-8s recovery tier: %s"
          % (g, p["status"], p["phase_PRIMARY"], p["gtex_blood_tpm"][:6],
             M[g]["reporting_tier"] if g in M else "-"))
print("[L02] of these %d lipid targets, without an eQTLGen whole-blood instrument: %d; with no "
      "eQTL data: %d; tested without a qualifying instrument: %d"
      % (len(LIPID), n_out, sum(P[g]["status"] == "NO_EQTL_DATA" for g in LIPID),
         sum(P[g]["status"] == "NO_QUALIFYING_INSTRUMENT" for g in LIPID)))
for g in ("HMGCR", "CETP", "LDLR", "LPL"):
    print("[L03] for contrast, %-6s eQTLGen: %s" % (g, P[g]["status"]))


def kind(x):
    return ("microRNA gene" if x["gene"].startswith("MIR") else
            "mitochondrial gene" if x["chromosome"] == "MT" else "nuclear protein-coding gene")


for label, rs in (("heart failure / cardiomyopathy",
                   [x for x in M.values()
                    if P[x["gene"]]["area_phase: Heart failure / cardiomyopathy"]]),
                  ("all areas", list(M.values()))):
    lk = [x for x in rs if x["reporting_tier"].startswith("locked")]
    pc = [x for x in lk if kind(x) == "nuclear protein-coding gene"]
    na = [x for x in pc if x["reporting_tier"] == "locked: never assayed"]
    print("[L04] %s: uninstrumentable %d; locked %d = %d microRNA genes + %d mitochondrial genes "
          "+ %d nuclear protein-coding genes; of the %d protein-coding, never assayed by any "
          "proteomic study: %d (%s); assayed without a candidate: %d (%s)"
          % (label, len(rs), len(lk), sum(kind(x) == "microRNA gene" for x in lk),
             sum(kind(x) == "mitochondrial gene" for x in lk), len(pc), len(pc), len(na),
             ", ".join(sorted(x["gene"] for x in na)) if len(na) < 15 else "list in Table S15",
             len(pc) - len(na),
             ", ".join(sorted(x["gene"] for x in pc if x not in na)) if len(pc) - len(na) < 15
             else "list in Table S15"))
for g in ("DMD", "PLN", "HCN4"):
    x = M[g]
    print("[L05] %-5s phase %s; tier %s; GTEx tissues tested %s; proteomic studies that assayed "
          "it: %s" % (g, x["phase_primary"], x["reporting_tier"], x["gtex_n_tissues_tested"],
                      x["pqtl_studies_assayed"] or "none"))
