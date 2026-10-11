#!/usr/bin/env python
"""Item 3 -- one definition of 'launched', applied to the two ORIGINAL sets, and the
21/52 vs 17/44 vs 20 reconciliation. The extended set is NOT computed here.

    python manuscript_prep/verify/launched_reconciliation.py

READ-ONLY.

DEFINITION.  phase(target, S) = max Minikel `ccatnum` over the target's pairs whose
indication is in set S.  launched(target, S)  <=>  phase(target, S) == 5.
The set that defines membership is the same set that defines phase.
"""
import collections
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
HERE = os.path.dirname(os.path.abspath(__file__))
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


topl = {r["indication_mesh_id"]: r["topl"]
        for r in load(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
frozen = load(ROOT + r"\data\derived\pp_dedup_frozen.tsv")
ann = {r["symbol"]: r["state"] for r in load(ROOT + r"\deliverables\d1e\classA_annotation.tsv")}
fix = {r["gene"]: r["state_after_resolution"] for r in csv.DictReader(
    open(os.path.join(os.path.dirname(HERE), "unmapped_resolution.csv"), encoding="utf-8"))}
SETS = {"CV-PRIMARY (C14 only)": lambda r: topl[r["indication_mesh_id"]] == "C14",
        "ORIGINAL 630 (C14 + lipid tokens)": lambda r: (
            topl[r["indication_mesh_id"]] == "C14"
            or any(k in r["indication_mesh_term"].lower() for k in TOK))}

for name, fn in SETS.items():
    pairs = [r for r in frozen if fn(r)]
    ph = collections.defaultdict(int)
    hist, act = collections.defaultdict(int), collections.defaultdict(int)
    for r in pairs:
        ph[r["gene"]] = max(ph[r["gene"]], int(r["ccatnum"]))
        # NB: `acatnum` is on a different scale (Launched = 6), so the text labels are used
        if r["hcat"] == "Launched":
            hist[r["gene"]] = 5
        if r["acat"] == "Launched":
            act[r["gene"]] = 5
    L = sorted(g for g in ph if ph[g] == 5)
    print("\n=== %s: %d targets, %d launched ===" % (name, len(ph), len(L)))
    for label, st in (("as shipped (12 unmapped)", lambda g: ann[g]),
                      ("after unmapped resolution (PROVISIONAL)",
                       lambda g: fix.get(g, ann[g]) if ann[g] == "UNMAPPED" else ann[g])):
        c = collections.Counter(st(g) for g in L)
        a = c["ABSENT"] + c["TESTED_NULL"]
        print("  %-42s %s" % (label, dict(c)))
        print("  %-42s Class A (absent+tested-null) %d / %d = %.1f%% ; absent only %d = %.1f%%"
              % ("", a, len(L), 100 * a / len(L), c["ABSENT"], 100 * c["ABSENT"] / len(L)))
    print("  launched by ACTIVE label (acat) only: %d ; by HISTORICAL (hcat) only: %d ; both: %d"
          % (sum(1 for g in L if act[g] == 5 and hist[g] != 5),
             sum(1 for g in L if hist[g] == 5 and act[g] != 5),
             sum(1 for g in L if hist[g] == 5 and act[g] == 5)))

# the '20' in the existing documents
atlas = {r["gene"]: r for r in csv.DictReader(
    open(ROOT + r"\deliverables\cv\cardiovascular_atlas.csv", encoding="utf-8"))}
mixed = sorted(g for g, r in atlas.items() if r["cv_primary"] == "True"
               and r["untestability_class"] == "CLASS_A" and r["max_phase_num"] == "5")
prim = [r for r in frozen if SETS["CV-PRIMARY (C14 only)"](r)]
php = collections.defaultdict(int)
for r in prim:
    php[r["gene"]] = max(php[r["gene"]], int(r["ccatnum"]))
print("\n'20' as computed by src/a2_robustness.py:84 (CV-PRIMARY membership, phase taken over "
      "the EXTENDED set): %d" % len(mixed))
print("  of these, NOT launched in any CV-PRIMARY indication: %s"
      % [(g, php[g]) for g in mixed if php[g] != 5])
