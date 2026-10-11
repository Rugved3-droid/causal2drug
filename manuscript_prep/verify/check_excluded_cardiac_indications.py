#!/usr/bin/env python
"""READ-ONLY. Which cardiac-sounding Minikel indications does the CV-EXTENDED rule exclude,
and what would they touch?  Quantifies DISCREPANCIES D16; changes nothing.

    python manuscript_prep/verify/check_excluded_cardiac_indications.py
"""
import collections
import csv

ROOT = r"D:\causal2drug"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
# keyword screen used ONLY to surface candidates for the user to look at
KEYS = ("arrhyth", "fibrillation", "tachycard", "long qt", "bradycard", "myocardial",
        "coronary", "angina", "heart", "cardi", "aortic", "ventric")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
areas = {r["topl"]: r["area"] for r in rd(ROOT + r"\data\raw\minikel_areas.tsv")}
topl = {r["indication_mesh_id"]: r["topl"]
        for r in rd(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
frozen = rd(ROOT + r"\data\derived\pp_dedup_frozen.tsv")
area = lambda r: areas.get(topl.get(r["indication_mesh_id"]))
is_cv = lambda r: area(r) == "cardiovascular" or any(
    k in r["indication_mesh_term"].lower() for k in TOK)
atlas = {r["gene"]: r for r in csv.DictReader(
    open(ROOT + r"\deliverables\cv\cardiovascular_atlas.csv", encoding="utf-8"))}

ex = [r for r in frozen if not is_cv(r)
      and any(k in r["indication_mesh_term"].lower() for k in KEYS)]
c = collections.defaultdict(lambda: [0, 0])
for r in ex:
    k = (r["indication_mesh_term"], area(r))
    c[k][0] += 1
    c[k][1] += r["ccatnum"] == "5"
print("%-34s %-16s %6s %9s" % ("indication", "Minikel area", "pairs", "launched"))
for k, v in sorted(c.items(), key=lambda x: -x[1][0]):
    print("%-34s %-16s %6d %9d" % (k[0], k[1], v[0], v[1]))

# restrict to the clearly cardiac ones (drop the virus, the surgical procedure)
CARD = {"Myocardial Infarction", "Angina Pectoris", "Atrial Fibrillation",
        "Arrhythmias, Cardiac", "Angina, Unstable", "Cardiomyopathy, Dilated",
        "Ventricular Fibrillation", "Edema, Cardiac", "Tachycardia, Ventricular",
        "Tachycardia, Supraventricular", "Diabetic Cardiomyopathies", "Long QT Syndrome",
        "Tachycardia", "Heart Defects, Congenital", "Bradycardia"}
ex2 = [r for r in ex if r["indication_mesh_term"] in CARD]
g = {r["gene"] for r in ex2}
gl = {r["gene"] for r in ex2 if r["ccatnum"] == "5"}
print("\n15 clearly cardiac excluded indications: pairs %d, targets %d" % (len(ex2), len(g)))
print("  targets not in the 630-row atlas at all: %d" % len(g - set(atlas)))
print("  targets with a LAUNCHED label in one of them: %d" % len(gl))
print("    not in the atlas at all: %s" % sorted(gl - set(atlas)))
up = sorted(x for x in gl if x in atlas and atlas[x]["max_phase_num"] != "5")
print("    in the atlas but NOT counted as launched there: %d -> %s" % (len(up), up))
for x in up:
    print("       %-8s atlas phase=%-11s state=%-12s class=%s"
          % (x, atlas[x]["max_clinical_phase"], atlas[x]["state"],
             atlas[x]["untestability_class"] or "-"))
