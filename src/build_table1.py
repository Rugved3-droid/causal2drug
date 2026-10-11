#!/usr/bin/env python
"""Table 1 -- the 21 launched and 18 Phase III CLASS A cardiovascular targets.

Example agents come from Open Targets 26.06 drug_molecule + drug_mechanism_of_action:
for each target, the approved (or highest-stage) drug acting on it, preferring an approved
small molecule / biologic with a recognisable name. Indication is the Minikel MeSH term of
the CV pair that reached the highest phase.
"""
import collections
import csv
import glob
import json

import pyarrow.parquet as pq

OT = r"D:\causal2drug\data\raw\ot2606"
CV = r"D:\causal2drug\deliverables\cv"
D1B = r"D:\causal2drug\deliverables\d1b"
RAW = r"D:\causal2drug\data\raw"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
STAGE = {"APPROVAL": 5, "PHASE_3": 4, "PHASE_2_3": 4, "PHASE_2": 3, "PHASE_1_2": 2,
         "PHASE_1": 2, "PREAPPROVAL": 5, "IND": 1, "PRECLINICAL": 0,
         "EARLY_PHASE_1": 1, "UNKNOWN": 0}


def tsv(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def main():
    mol = pq.read_table(glob.glob(OT + r"\drug_molecule\*.parquet"),
                        columns=["id", "name", "maximumClinicalStage",
                                 "drugType"]).to_pylist()
    minfo = {m["id"]: m for m in mol}
    mech = pq.read_table(glob.glob(OT + r"\drug_mechanism_of_action\*.parquet")).to_pylist()
    by_target = collections.defaultdict(list)
    for r in mech:
        for t in (r["targets"] or []):
            for cid in (r["chemblIds"] or []):
                m = minfo.get(cid)
                if not m:
                    continue
                by_target[t].append((STAGE.get(m["maximumClinicalStage"], 0),
                                     m["name"], r["actionType"], m["drugType"]))

    atlas = list(csv.DictReader(open(CV + r"\cardiovascular_atlas.csv", encoding="utf-8")))
    coords = {r["symbol"]: r["ensembl_id"] for r in tsv(D1B + r"\gene_coords_grch37.tsv")}

    # highest-phase CV indication per target
    areas = {r["topl"]: r["area"] for r in tsv(RAW + r"\minikel_areas.tsv")}
    topl = {r["indication_mesh_id"]: r["topl"]
            for r in tsv(RAW + r"\minikel_indic_topl_match.tsv")}
    best_ind = {}
    for r in tsv(FROZEN):
        cv = (areas.get(topl.get(r["indication_mesh_id"])) == "cardiovascular"
              or any(k in r["indication_mesh_term"].lower() for k in TOK))
        if not cv or not r["ccatnum"].isdigit():
            continue
        g, n = r["gene"], int(r["ccatnum"])
        if g not in best_ind or n > best_ind[g][0]:
            best_ind[g] = (n, r["indication_mesh_term"])

    rows = []
    for r in atlas:
        if r["untestability_class"] != "CLASS_A":
            continue
        if r["max_phase_num"] not in ("4", "5"):
            continue
        ensg = coords.get(r["gene"], "")
        SALTS = (" Hydrochloride", " Sodium", " Fumarate", " Malate", " Hydrobromide",
                 " Sulfate", " Mesylate", " Tartrate", " Citrate", " Maleate",
                 " Besylate", " Acetate", " Potassium", " Calcium", " Dihydrate")
        cand = by_target.get(ensg, [])
        # keep only names that survive ASCII round-trip (a few OT names carry glyphs
        # that render as replacement characters in a journal table)
        cand = [c for c in cand
                if c[1] and c[1].encode("ascii", "ignore").decode() == c[1]]
        best_stage = max((c[0] for c in cand), default=None)
        cand = [c for c in cand if c[0] == best_stage]
        # deterministic and prefers the parent over a salt form
        def clean(n):
            for sfx in SALTS:
                if n.endswith(sfx):
                    return n[: -len(sfx)]
            return n
        cand.sort(key=lambda c: (len(clean(c[1])), clean(c[1]).lower()))
        agent = clean(cand[0][1]) if cand else ""
        action = cand[0][2] if cand else ""
        dtype = cand[0][3] if cand else ""
        n_agents = len(by_target.get(ensg, []))
        prim_t, prim_v = "", None
        for lab, key in (("heart LV", "gtex_heart_LV_tpm"), ("liver", "gtex_liver_tpm"),
                         ("aorta", "gtex_artery_aorta_tpm")):
            v = r[key]
            if v not in ("", None):
                v = float(v)
                if prim_v is None or v > prim_v:
                    prim_v, prim_t = v, lab
        rows.append(dict(
            gene=r["gene"],
            example_agent=(agent.title() if agent else "not in Open Targets"),
            n_agents_for_target=n_agents,
            drug_class=("%s %s" % (dtype, action.lower()) if action else "—"),
            indication=best_ind.get(r["gene"], (0, "—"))[1],
            highest_phase=r["max_clinical_phase"],
            blood_tpm=(round(float(r["gtex_blood_tpm"]), 2)
                       if r["gtex_blood_tpm"] not in ("", None) else None),
            primary_tissue=prim_t,
            primary_tissue_tpm=(round(prim_v, 1) if prim_v is not None else None),
            failure_class=("absent from blood eQTL resource" if r["state"] == "ABSENT"
                           else "tested, no usable signal"),
            required_modality=r["required_modality"],
            _sort=(0 if r["max_phase_num"] == "5" else 1, r["gene"])))
    rows.sort(key=lambda x: x["_sort"])
    for x in rows:
        del x["_sort"]
    cols = list(rows[0].keys())
    with open(CV + r"\table1_classA_launched_phase3.csv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    n5 = sum(1 for x in rows if x["highest_phase"] == "Launched")
    print("Table 1: %d rows (%d launched, %d Phase III)" % (len(rows), n5, len(rows) - n5))
    print("%-9s %-22s %-30s %-11s %-8s %-16s %s"
          % ("gene", "example agent", "indication", "phase", "bloodTPM", "primary tissue",
             "required modality"))
    for x in rows:
        print("%-9s %-22s %-30s %-11s %-8s %-16s %s"
              % (x["gene"], x["example_agent"][:22], x["indication"][:30],
                 x["highest_phase"],
                 ("%.2f" % x["blood_tpm"]) if x["blood_tpm"] is not None else "-",
                 "%s %s" % (x["primary_tissue"], x["primary_tissue_tpm"]),
                 x["required_modality"][:34]))


if __name__ == "__main__":
    main()
