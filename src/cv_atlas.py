#!/usr/bin/env python
"""CARDIOVASCULAR TARGET UNTESTABILITY ATLAS.

CV SET RULE -- declared before applying (see CARDIOVASCULAR_ATLAS.md):

  classification source : Minikel's own `areas.tsv` x `indic_topl_match.tsv`
                          (MeSH top-level tree code -> therapeutic area; all 1,014
                          indications map, zero misses)
  CV-PRIMARY   : therapeutic area == "cardiovascular" (MeSH C14)
  CV-EXTENDED  : CV-PRIMARY  UNION  indications whose MeSH term matches any of
                 {hyperlipid, hypercholesterol, dyslipid, hyperlipoprotein,
                  hypertriglycerid, lipoprotein, atheroscler, arterioscler}

  The extension is mechanical (a fixed token list applied to the MeSH term) and is
  declared because Minikel classifies the atherogenic lipid disorders as *metabolic*
  (C18), not cardiovascular. Hypercholesterolemia and Hyperlipoproteinemia Type II are
  the PCSK9 and statin indications; a cardiovascular atlas that excluded them would
  exclude its own anchor example. CV-EXTENDED is the atlas definition; CV-PRIMARY is
  reported alongside so the effect of the choice is visible.
"""
import collections
import csv
import gzip
import json
import math
import os

import numpy as np
from scipy import stats

RAW = r"D:\causal2drug\data\raw"
D1B = r"D:\causal2drug\deliverables\d1b"
D1E = r"D:\causal2drug\deliverables\d1e"
D1F = r"D:\causal2drug\deliverables\d1f"
OUT = r"D:\causal2drug\deliverables\cv"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"

TOKENS = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
          "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
CARDIAC = [("QTD000256", "heart_LV", "primary"), ("QTD000251", "heart_atrial", "secondary"),
           ("QTD000131", "artery_aorta", "secondary"), ("QTD000266", "liver", "secondary")]
PHASE = {"1": "Preclinical", "2": "Phase I", "3": "Phase II", "4": "Phase III",
         "5": "Launched"}


def tsv(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def main():
    os.makedirs(OUT, exist_ok=True)
    areas = {r["topl"]: r["area"] for r in tsv(RAW + r"\minikel_areas.tsv")}
    topl = {r["indication_mesh_id"]: r["topl"]
            for r in tsv(RAW + r"\minikel_indic_topl_match.tsv")}
    frozen = tsv(FROZEN)
    for r in frozen:
        r["area"] = areas.get(topl.get(r["indication_mesh_id"]), "UNMAPPED")
        t = r["indication_mesh_term"].lower()
        r["cv_primary"] = r["area"] == "cardiovascular"
        r["cv_extended"] = r["cv_primary"] or any(k in t for k in TOKENS)

    cvp = [r for r in frozen if r["cv_primary"]]
    cve = [r for r in frozen if r["cv_extended"]]
    tp = {r["gene"] for r in cvp}
    te = {r["gene"] for r in cve}
    print("CV-PRIMARY  pairs %5d   targets %4d" % (len(cvp), len(tp)))
    print("CV-EXTENDED pairs %5d   targets %4d   (+%d targets from lipid indications)"
          % (len(cve), len(te), len(te - tp)))

    # ---- per-target state + annotation ----
    ann = {r["symbol"]: r for r in tsv(D1E + r"\classA_annotation.tsv")}
    clump = {r["symbol"]: r for r in tsv(D1B + r"\clumped_instruments_final.tsv")}
    coords = {r["symbol"]: r for r in tsv(D1B + r"\gene_coords_grch37.tsv")}
    direction = json.load(open(D1F + r"\target_direction.json"))

    # ---- cardiac tissue instrumentation (top nominal cis p < 5e-8 per dataset) ----
    card = {}
    for qtd, label, _role in CARDIAC:
        hit = set()
        with gzip.open(RAW + ("\\%s.permuted.tsv.gz" % qtd), "rt", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    if float(r["pvalue"]) < 5e-8:
                        hit.add(r["molecular_trait_id"])
                except ValueError:
                    pass
        card[label] = hit
        print("  %-14s genes with top cis p<5e-8: %d" % (label, len(hit)))

    # ---- GTEx median TPM per tissue, for the modality rule ----
    tpm = {}
    with gzip.open(RAW + r"\gtex_v10_median_tpm.gct.gz", "rt", encoding="utf-8") as fh:
        fh.readline(); fh.readline()
        hdr = fh.readline().rstrip("\n").split("\t")
        want = {t: hdr.index(t) for t in
                ("Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Liver",
                 "Artery_Aorta", "Artery_Coronary", "Whole_Blood") if t in hdr}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = f[0].split(".")[0]
            try:
                tpm[g] = {k: float(f[i]) for k, i in want.items()}
                tpm[g]["_max"] = max(float(f[i]) for i in range(2, len(f)))
                tpm[g]["_argmax"] = hdr[2 + int(np.argmax([float(f[i])
                                                           for i in range(2, len(f))]))]
            except (ValueError, IndexError):
                pass

    # ---- max clinical phase per target within CV-EXTENDED ----
    maxphase = collections.defaultdict(int)
    for r in cve:
        if r["ccatnum"].isdigit():
            maxphase[r["gene"]] = max(maxphase[r["gene"]], int(r["ccatnum"]))

    def modality(sym, state, ensg):
        """Declared rule; UNKNOWN where it cannot be determined."""
        if sym == "IL6R":
            return "sQTL or pQTL (shedding variant rs2228145; total mRNA cannot represent it)"
        if state == "INSTRUMENTED":
            return "n/a - instrumented in blood"
        t = tpm.get(ensg)
        a = ann.get(sym, {})
        if t:
            if t.get("Heart_Left_Ventricle", 0) >= 1 or t.get("Heart_Atrial_Appendage", 0) >= 1:
                return "cardiac eQTL"
            if t.get("Liver", 0) >= 1:
                return "liver eQTL"
            if t.get("Artery_Aorta", 0) >= 1 or t.get("Artery_Coronary", 0) >= 1:
                return "vascular eQTL"
        if a.get("secreted") == "YES":
            return "pQTL"
        if t and t.get("_max", 0) >= 1:
            return "tissue eQTL (%s)" % t["_argmax"]
        return "UNKNOWN"

    rows = []
    for sym in sorted(te):
        a = ann.get(sym, {})
        state = a.get("state", "UNMAPPED")
        ensg = coords.get(sym, {}).get("ensembl_id", "")
        cl = clump.get(sym, {})
        t = tpm.get(ensg, {})
        card_hit = {lab: (ensg in card[lab]) for lab in card}
        rows.append(dict(
            gene=sym, ensembl_id=ensg,
            cv_primary=sym in tp, cv_extended=True,
            n_cv_pairs=sum(1 for r in cve if r["gene"] == sym),
            max_clinical_phase=PHASE.get(str(maxphase.get(sym, 0)), ""),
            max_phase_num=maxphase.get(sym, 0),
            approved_drug_direction=direction.get(ensg, ""),
            state=state,
            untestability_class=("CLASS_A" if state in ("ABSENT", "TESTED_NULL")
                                 else ("CLASS_B" if sym == "IL6R" else "")),
            n_independent_instruments=cl.get("n_independent", ""),
            best_eqtl_p=cl.get("best_p", ""),
            instrumented_blood=(state == "INSTRUMENTED"),
            instrumented_heart_LV=card_hit.get("heart_LV", False),
            instrumented_heart_atrial=card_hit.get("heart_atrial", False),
            instrumented_artery_aorta=card_hit.get("artery_aorta", False),
            instrumented_liver=card_hit.get("liver", False),
            instrumented_cardiac_union=any(card_hit.values()),
            tissue_restricted=a.get("tissue_restricted", "NA"),
            secreted=a.get("secreted", "NA"),
            gtex_blood_tpm=a.get("blood_tpm", ""),
            gtex_heart_LV_tpm=t.get("Heart_Left_Ventricle", ""),
            gtex_liver_tpm=t.get("Liver", ""),
            gtex_artery_aorta_tpm=t.get("Artery_Aorta", ""),
            required_modality=modality(sym, state, ensg)))

    cols = list(rows[0].keys())
    with open(OUT + r"\cardiovascular_atlas.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    try:
        import pyarrow as pa, pyarrow.parquet as pq
        pq.write_table(pa.Table.from_pylist(rows), OUT + r"\cardiovascular_atlas.parquet")
        pqok = True
    except Exception as e:
        pqok = False
        print("parquet failed: %r" % e)
    print("\nwrote atlas: %d CV targets (csv%s)" % (len(rows), " + parquet" if pqok else ""))

    # ---------------- summary statistics ----------------
    st = collections.Counter(r["state"] for r in rows)
    print("\n=== 1. CLASS BREAKDOWN (CV-EXTENDED, N=%d targets) ===" % len(rows))
    for k in ("ABSENT", "TESTED_NULL", "INSTRUMENTED", "UNMAPPED"):
        print("  %-14s %4d  (%.1f%%)" % (k, st[k], 100 * st[k] / len(rows)))
    nA = st["ABSENT"] + st["TESTED_NULL"]
    print("  CLASS_A total  %4d  (%.1f%%)" % (nA, 100 * nA / len(rows)))

    print("\n=== 2. TISSUE YIELD, CV targets vs genome-wide ===")
    nb = sum(1 for r in rows if r["instrumented_blood"])
    print("  whole blood (eQTLGen)     %4d / %d = %.1f%%   [genome-wide 62.6%%]"
          % (nb, len(rows), 100 * nb / len(rows)))
    for lab, key in (("heart LV (primary)", "instrumented_heart_LV"),
                     ("heart atrial appendage", "instrumented_heart_atrial"),
                     ("artery aorta", "instrumented_artery_aorta"),
                     ("liver", "instrumented_liver"),
                     ("cardiac union (4 tissues)", "instrumented_cardiac_union")):
        n = sum(1 for r in rows if r[key])
        print("  %-25s %4d / %d = %.1f%%" % (lab, n, len(rows), 100 * n / len(rows)))

    print("\n=== 3. CLASS A ENRICHMENT WITHIN CV ===")
    inh = [r for r in rows if r["tissue_restricted"] in ("YES", "NO")]
    a = sum(1 for r in inh if r["state"] == "ABSENT" and r["tissue_restricted"] == "YES")
    b = sum(1 for r in inh if r["state"] == "ABSENT" and r["tissue_restricted"] == "NO")
    c = sum(1 for r in inh if r["state"] == "INSTRUMENTED" and r["tissue_restricted"] == "YES")
    d = sum(1 for r in inh if r["state"] == "INSTRUMENTED" and r["tissue_restricted"] == "NO")
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a + .5, b + .5, c + .5, d + .5)
    lor = math.log((aa * dd) / (bb * cc)); se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    print("  tissue-restricted  ABSENT %d/%d   INSTRUMENTED %d/%d" % (a, a + b, c, c + d))
    print("  OR = %.3f  95%% CI [%.3f, %.3f]  Fisher p = %.4g"
          % (orr, math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se), p))
    print("  [genome-wide primary: OR 6.90, 95%% CI 5.11-9.31, p = 2.9e-51]")

    json.dump(dict(n_cv_primary_pairs=len(cvp), n_cv_primary_targets=len(tp),
                   n_cv_extended_pairs=len(cve), n_cv_extended_targets=len(te),
                   state_counts=dict(st),
                   blood=nb,
                   cardiac_union=sum(1 for r in rows if r["instrumented_cardiac_union"]),
                   heart_lv=sum(1 for r in rows if r["instrumented_heart_LV"]),
                   cv_classA_or=orr, cv_classA_p=p,
                   cv_classA_table=[[a, b], [c, d]]),
              open(OUT + r"\cv_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
