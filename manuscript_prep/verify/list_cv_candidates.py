#!/usr/bin/env python
"""READ-ONLY. Lists the Minikel indications relevant to extending the CV set.
Pair counts and phase labels only -- no instrument state, no class.

    python manuscript_prep/verify/list_cv_candidates.py
"""
import collections
import csv
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
KEYS = ("cardi", "heart", "vascul", "arter", "aort", "venous", "vein", "thromb", "embol",
        "isch", "infarct", "hemorrhag", "aneurysm", "hypertens", "hypotens", "shock",
        "angio", "edema", "lipid", "cholesterol", "valve", "syncope", "atri", "ventric",
        "claudication", "phleb", "varic", "stroke")


def load(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


ind = {r["indication_mesh_id"]: r for r in load(ROOT + r"\data\raw\minikel_indic.tsv")}
topl = {r["indication_mesh_id"]: r["topl"]
        for r in load(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
frozen = load(ROOT + r"\data\derived\pp_dedup_frozen.tsv")
n = collections.Counter(r["indication_mesh_id"] for r in frozen)
n3 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "4")
n5 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "5")
tg = collections.defaultdict(set)
for r in frozen:
    tg[r["indication_mesh_id"]].add(r["gene"])


def term(i):
    return ind[i]["indication_mesh_term"]


def areas(i):
    return [a.strip() for a in ind[i]["areas"].split(",")]


def current(i):
    return topl[i] == "C14" or any(k in term(i).lower() for k in TOK)


def show(title, ids):
    print("\n== %s (%d terms, %d pairs) ==" % (title, len(ids), sum(n[i] for i in ids)))
    print("mesh_id|term|minikel_areas|topl|pairs|phaseIII_pairs|launched_pairs")
    for i in sorted(ids, key=lambda i: -n[i]):
        print("%s|%s|%s|%s|%d|%d|%d" % (i, term(i), ind[i]["areas"], topl[i], n[i], n3[i],
                                         n5[i]))


cur = [i for i in n if current(i)]
A = [i for i in n if not current(i) and topl[i] == "C23" and "cardiovascular" in areas(i)]
B = [i for i in n if not current(i) and topl[i] != "C23" and "cardiovascular" in areas(i)]
C = [i for i in n if not current(i) and topl[i] == "C23" and "cardiovascular" not in areas(i)]
D = [i for i in n if not current(i) and topl[i] != "C23" and "cardiovascular" not in areas(i)
     and any(k in term(i).lower() for k in KEYS)]
show("CURRENT atlas set (C14 top-level, or lipid token)", cur)
print("current-set terms whose Minikel areas do NOT contain 'cardiovascular': %s"
      % [term(i) for i in cur if "cardiovascular" not in areas(i)])
show("A: filed under C23 signs/symptoms AND Minikel areas contains cardiovascular", A)
show("B: filed elsewhere (not C14, not C23) AND Minikel areas contains cardiovascular", B)
show("C: every other C23 signs/symptoms term (areas lacks cardiovascular) - for eyeballing", C)
show("D: filed elsewhere, areas lacks cardiovascular, term matches a vascular keyword", D)

for name, ids in (("current", cur), ("current + A", cur + A), ("current + A + B", cur + A + B)):
    g = set().union(*[tg[i] for i in ids])
    print("\n%s: terms %d, pairs %d, distinct targets %d"
          % (name, len(ids), sum(n[i] for i in ids), len(g)))
