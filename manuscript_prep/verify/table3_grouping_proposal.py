#!/usr/bin/env python
"""Item 7 -- PROPOSED indication-area grouping for Table 3. Pair counts only.

    python manuscript_prep/verify/table3_grouping_proposal.py

READ-ONLY. No instrument state or class is read. Writes
manuscript_prep/table3_grouping_proposal.csv. The mapping below is a PROPOSAL awaiting
approval; the script fails if any term in the candidate sets is left unmapped.

Two columns are proposed:
  area_9   the seven areas in the brief plus two that the data seem to need
           (ischaemic heart disease; cerebrovascular)
  area_7   the same terms forced into exactly the seven areas in the brief
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
LIP, IHD, HF, ARR, HTN, PH, THR, CBV, OTH = (
    "Lipids / atherosclerosis", "Ischaemic heart disease", "Heart failure / cardiomyopathy",
    "Arrhythmia", "Hypertension", "Pulmonary hypertension", "Thrombosis",
    "Cerebrovascular", "Other")

# term -> (area_9, area_7, note)
M = {
    # ---- lipids / atherosclerosis
    "Atherosclerosis": (LIP, LIP, ""),
    "Hyperlipidemias": (LIP, LIP, ""),
    "Hypercholesterolemia": (LIP, LIP, ""),
    "Hypertriglyceridemia": (LIP, LIP, ""),
    "Hyperlipoproteinemia Type II": (LIP, LIP, ""),
    "Dyslipidemias": (LIP, LIP, ""),
    "Hyperlipoproteinemia Type I": (LIP, LIP, ""),
    "Hypoalphalipoproteinemias": (LIP, LIP, ""),
    # ---- ischaemic heart disease (no home in the 7-area list)
    "Myocardial Infarction": (IHD, LIP, "7-area fallback treats coronary disease as atherosclerotic"),
    "Acute Coronary Syndrome": (IHD, LIP, "same"),
    "Angina Pectoris": (IHD, LIP, "same"),
    "Angina, Unstable": (IHD, LIP, "same"),
    "Coronary Restenosis": (IHD, LIP, "post-intervention; could be Other"),
    "Coronary Thrombosis": (IHD, THR, "could equally be Thrombosis in both schemes"),
    "Reperfusion Injury": (OTH, OTH, "not organ-specific in MeSH; includes non-cardiac reperfusion"),
    # ---- heart failure / cardiomyopathy
    "Heart Failure": (HF, HF, ""),
    "Cardiomyopathies": (HF, HF, ""),
    "Cardiomyopathy, Hypertrophic": (HF, HF, ""),
    "Endomyocardial Fibrosis": (HF, HF, ""),
    "Edema, Cardiac": (HF, HF, ""),
    "Myocarditis": (HF, HF, "could be Other"),
    # ---- arrhythmia
    "Atrial Fibrillation": (ARR, ARR, ""),
    "Arrhythmias, Cardiac": (ARR, ARR, ""),
    "Ventricular Fibrillation": (ARR, ARR, ""),
    "Tachycardia, Ventricular": (ARR, ARR, ""),
    "Tachycardia, Supraventricular": (ARR, ARR, ""),
    "Long QT Syndrome": (ARR, ARR, ""),
    "Tachycardia": (ARR, ARR, ""),
    "Bradycardia": (ARR, ARR, ""),
    # ---- hypertension
    "Hypertension": (HTN, HTN, ""),
    # ---- pulmonary hypertension
    "Hypertension, Pulmonary": (PH, PH, ""),
    # ---- thrombosis
    "Thrombosis": (THR, THR, ""),
    "Venous Thrombosis": (THR, THR, ""),
    "Pulmonary Embolism": (THR, THR, ""),
    "Thrombophlebitis": (THR, THR, ""),
    # ---- cerebrovascular (no home in the 7-area list)
    "Stroke": (CBV, OTH, "identical target list to Brain Ischemia in Minikel"),
    "Brain Ischemia": (CBV, OTH, "identical target list to Stroke in Minikel"),
    "Cerebral Hemorrhage": (CBV, OTH, ""),
    "Subarachnoid Hemorrhage": (CBV, OTH, ""),
    "Cerebral Infarction": (CBV, OTH, ""),
    "Intracranial Thrombosis": (CBV, THR, ""),
    "Vasospasm, Intracranial": (CBV, OTH, ""),
    "CADASIL": (CBV, OTH, "monogenic small-vessel disease"),
    # ---- other
    "Cardiovascular Diseases": (OTH, OTH, "generic parent term; 170 pairs"),
    "Peripheral Vascular Diseases": (OTH, OTH, ""),
    "Intermittent Claudication": (OTH, OTH, "peripheral arterial"),
    "Vasculitis": (OTH, OTH, ""),
    "Hypotension": (OTH, OTH, "not Hypertension; opposite pharmacology"),
    "Retinal Vein Occlusion": (OTH, OTH, "ophthalmic vascular"),
    "Retinal Artery Occlusion": (OTH, OTH, "ophthalmic vascular"),
    "Optic Neuropathy, Ischemic": (OTH, OTH, "ophthalmic vascular"),
    "Retinal Telangiectasis": (OTH, OTH, "ophthalmic vascular"),
    "Raynaud Disease": (OTH, OTH, ""),
    "Vascular Diseases": (OTH, OTH, "generic parent term"),
    "Venous Insufficiency": (OTH, OTH, ""),
    "Aneurysm": (OTH, OTH, ""),
    "Aortic Valve Stenosis": (OTH, OTH, "valvular"),
    "Capillary Leak Syndrome": (OTH, OTH, ""),
    "Thromboangiitis Obliterans": (OTH, OTH, ""),
    "Erythromelalgia": (OTH, OTH, ""),
    "Endocarditis": (OTH, OTH, ""),
    "Pericarditis": (OTH, OTH, ""),
    # ---- group B (only used if the CV set is widened beyond rule A)
    "Cardiomyopathy, Dilated": (HF, HF, "GROUP B"),
    "Diabetic Cardiomyopathies": (HF, HF, "GROUP B"),
    "Dementia, Vascular": (CBV, OTH, "GROUP B"),
    "Cerebral Amyloid Angiopathy": (CBV, OTH, "GROUP B"),
}


def load(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


ind = {r["indication_mesh_id"]: r for r in load(ROOT + r"\data\raw\minikel_indic.tsv")}
topl = {r["indication_mesh_id"]: r["topl"]
        for r in load(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
frozen = load(ROOT + r"\data\derived\pp_dedup_frozen.tsv")
n = collections.Counter(r["indication_mesh_id"] for r in frozen)
n5 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "5")
n4 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "4")
rows = []
for i in n:
    t = ind[i]["indication_mesh_term"]
    tag = "cardiovascular" in [a.strip() for a in ind[i]["areas"].split(",")]
    s0 = topl[i] == "C14" or any(k in t.lower() for k in TOK)
    if s0:
        src = "original set: C14 top-level" if topl[i] == "C14" else "original set: lipid token"
    elif topl[i] == "C23" and tag:
        src = "RULE A addition: filed C23, tagged cardiovascular"
    elif tag:
        src = "GROUP B (not proposed): filed %s, tagged cardiovascular" % topl[i]
    else:
        continue
    a9, a7, note = M.get(t, (OTH, OTH, "GROUP B") if src.startswith("GROUP B") else (None,) * 3)
    if a9 is None:
        sys.exit("UNMAPPED TERM: " + t)
    rows.append(dict(mesh_id=i, mesh_term=t, minikel_areas=ind[i]["areas"],
                     minikel_top_level=topl[i], source=src, pairs=n[i],
                     phase3_pairs=n4[i], launched_pairs=n5[i], proposed_area_9=a9,
                     area_if_only_7=a7, note=note))
rows.sort(key=lambda r: (r["source"].startswith("GROUP B"), r["proposed_area_9"], -r["pairs"]))
out = os.path.join(os.path.dirname(HERE), "table3_grouping_proposal.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("wrote %s (%d terms)" % (out, len(rows)))
prop = [r for r in rows if not r["source"].startswith("GROUP B")]
print("\nProposed extended set (original + rule A): %d terms, %d pairs"
      % (len(prop), sum(r["pairs"] for r in prop)))
for col in ("proposed_area_9", "area_if_only_7"):
    print("\n%-34s %6s %7s %10s" % (col, "terms", "pairs", "launched"))
    c = collections.defaultdict(lambda: [0, 0, 0])
    for r in prop:
        c[r[col]][0] += 1
        c[r[col]][1] += r["pairs"]
        c[r[col]][2] += r["launched_pairs"]
    for k, v in sorted(c.items(), key=lambda x: -x[1][1]):
        print("%-34s %6d %7d %10d" % (k, v[0], v[1], v[2]))
print("\n| area | MeSH term | pairs | Phase III | launched | source | note |")
print("|---|---|---|---|---|---|---|")
for r in prop:
    print("| %s | %s | %d | %d | %d | %s | %s |"
          % (r["proposed_area_9"], r["mesh_term"], r["pairs"], r["phase3_pairs"],
             r["launched_pairs"], "new" if r["source"].startswith("RULE A") else "orig",
             r["note"] + (" [7-area: %s]" % r["area_if_only_7"]
                          if r["area_if_only_7"] != r["proposed_area_9"] else "")))
