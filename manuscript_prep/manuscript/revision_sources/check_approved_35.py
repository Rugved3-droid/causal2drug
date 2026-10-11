#!/usr/bin/env python
"""Apply the Methods definition of "approved" to all strict-verified launched targets.

    python manuscript_prep/manuscript/revision_sources/check_approved_35.py \
        > manuscript_prep/manuscript/revision_sources/approved_35_output.txt

Definition (Methods): ever approved by the FDA or the EMA for a cardiovascular indication.
Operational rule, the same lookup as Table 2:
  * approval record  = at least one product in Drugs@FDA with a marketing status other than
                       "None (Tentative Approval)" (so discontinued products count), OR an
                       entry in the EMA medicines table whose status shows it was authorised
                       (Authorised, Withdrawn, Suspended, Lapsed, Expired, Revoked). EMA
                       "Refused", "Application withdrawn" and "Opinion" are not approvals.
  * CV indication    = the drug's approved, strictly cardiovascular indication in ChEMBL
                       (the strict rule that produced the 35).
  * A target PASSES if at least one of its strict-rule drugs has an approval record.
  * Secondary check, reported and not used for the pass: whether the agency's own indication
    text (FDA label "indications and usage"; EMA therapeutic area / indication fields)
    contains a cardiovascular term.
Writes manuscript/APPROVED_35_CHECK.csv (one row per target x drug x agency) and
manuscript/revision_sources/approved_35_by_target.csv.
"""
import collections
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
C = os.path.join(HERE, "drug_status")
ACCESSED = "2026-10-11"
EMA_APPROVED = ("authorised", "withdrawn", "suspended", "lapsed", "expired", "revoked")
SALT = [" PHOSPHORIC ACID", " MEDOXOMIL", " CILEXETIL", " ETEXILATE", " ANHYDROUS", " SODIUM",
        " HYDROCHLORIDE", " MESYLATE", " POTASSIUM", " CALCIUM"]
CVWORDS = ["hypertens", "heart failure", "cardiac", "myocardial", "coronary", "angina",
           "thromb", "embol", "stroke", "atrial fibrillation", "arrhythm", "tachycardia",
           "cholesterol", "hyperlipid", "lipid", "triglycerid", "atheroscl", "vascular",
           "ischemi", "ischaemi", "infarction", "pericarditis", "hyponatremia", "edema",
           "oedema", "anticoagul", "claudication", "ductus", "vasculitis", "polyangiitis",
           "shock", "blood pressure", "chylomicron"]


def get(url, name):
    f = os.path.join(C, name)
    if not os.path.exists(f):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                open(f, "wb").write(r.read())
        except urllib.error.HTTPError as e:
            open(f, "w").write(json.dumps({"error": {"code": e.code}}))
        time.sleep(0.35)
    try:
        return json.load(open(f, encoding="utf-8"))
    except Exception:                                            # noqa: BLE001
        return {"error": {"code": "NOT_JSON"}}


def key(s):
    return re.sub(r"[^A-Z0-9]+", "_", s.upper()).strip("_")


P = [x for x in csv.DictReader(open(os.path.join(MP, "analysis", "out", "per_target_primary.csv"),
                                    encoding="utf-8"))
     if x["phase_PRIMARY"] in ("Launched", "Phase III")
     and x["chembl_verified_strict"] == "True"]
PHASE = {x["gene"]: x["phase_PRIMARY"] for x in P}
P.sort(key=lambda x: (x["phase_PRIMARY"] != "Launched", x["gene"]))
L = list(csv.DictReader(open(os.path.join(MP, "drug_links", "chembl_target_drug_indication.csv"),
                             encoding="utf-8-sig")))
cvind = collections.defaultdict(set)
for r in L:
    if r["row_strict"] == "True":
        cvind[(r["gene"], r["drug_name"])].add(r["indication_mesh_heading"])
ema = get("https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json",
          "ema_medicines.json")
ema_rows = ema.get("data", []) if isinstance(ema, dict) else ema
print("[Q01] strict-verified launched targets: %d; target-drug pairs: %d; EMA table rows: %d"
      % (sum(1 for x in P if x["phase_PRIMARY"] == "Launched"),
         sum(len(x["chembl_strict_drugs"].split("; ")) for x in P
             if x["phase_PRIMARY"] == "Launched"), len(ema_rows)))
print("[Q01] strict-verified Phase III targets also checked: %s"
      % ", ".join(x["gene"] for x in P if x["phase_PRIMARY"] == "Phase III"))

rows, by_t = [], collections.OrderedDict()
for x in P:
    g = x["gene"]
    T = by_t.setdefault(g, dict(target=g, eqtlgen_status=x["status"], n_drugs=0,
                                drugs_with_record=[], drugs_text_confirmed=[],
                                drugs_without_record=[]))
    for drug in x["chembl_strict_drugs"].split("; "):
        T["n_drugs"] += 1
        chembl_cv = "; ".join(sorted(cvind.get((g, drug), []))) or "NOT FOUND in drug_links"
        # ---- FDA
        tried, st, apps, stem = [], {}, [], drug
        cands = [drug] + [drug.replace(s_, "") for s_ in SALT if drug.endswith(s_)]
        for q in cands:
            tried.append(q)
            u = "https://api.fda.gov/drug/drugsfda.json?" + urllib.parse.urlencode(
                {"search": 'products.active_ingredients.name:"%s"' % q, "limit": 99})
            j = get(u, "drugsfda_%s.json" % key(q))
            st, apps = {}, []
            for res in j.get("results", []):
                for p in res.get("products", []):
                    if any(q in a["name"].upper() for a in p.get("active_ingredients", [])):
                        ms = p.get("marketing_status", "?")
                        st[ms] = st.get(ms, 0) + 1
                        apps.append(res.get("application_number"))
            if st:
                stem = q
                break
        real = {k: v for k, v in st.items() if "Tentative" not in k}
        fda_rec = bool(real)
        fda_status = ("marketed" if any(k in real for k in ("Prescription", "Over-the-counter"))
                      else "discontinued" if real else "no approval record")
        lj = get("https://api.fda.gov/drug/label.json?" + urllib.parse.urlencode(
            {"search": 'openfda.generic_name:"%s"' % stem, "limit": 1}), "label_%s.json" % key(stem))
        ind_txt, setid = "", ""
        if lj.get("results"):
            ind_txt = " ".join(lj["results"][0].get("indications_and_usage", [""]))
            setid = lj["results"][0].get("set_id", "")
        fda_cv = [w for w in CVWORDS if w in ind_txt.lower()]
        rows.append(dict(
            target=g, drug=drug.title(), agency="FDA",
            approval_record=("yes: " + "; ".join("%s %d" % kv for kv in sorted(st.items()))
                             + " (" + ", ".join(sorted(set(apps))[:3]) + ")") if st else "no record",
            cv_indication_chembl=chembl_cv,
            cv_indication_agency_text=(("yes: " + ", ".join(fda_cv[:4])) if fda_cv else
                                       "label text not retrieved" if not ind_txt else
                                       "NO cardiovascular term found in the label indication"),
            agency_indication_excerpt=re.sub(r"\s+", " ", ind_txt)[:260],
            status=fda_status,
            source_link=("https://api.fda.gov/drug/drugsfda.json?search=products."
                         "active_ingredients.name:%22" + urllib.parse.quote(stem) + "%22"),
            label_link=("https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=" + setid)
            if setid else "", query_used=stem, accessed=ACCESSED))
        # ---- EMA
        em = [e for e in ema_rows
              if stem.lower() in (str(e.get("international_non_proprietary_name_common_name", ""))
                                  + " | " + str(e.get("active_substance", ""))).lower()
              and str(e.get("category", "Human")).lower().startswith("human")]
        sts = sorted({str(e.get("medicine_status", "?")) for e in em})
        ema_rec = any(s_.lower() in EMA_APPROVED for s_ in sts)
        etxt = " ".join(str(e.get(k, "")) for e in em
                        for k in ("therapeutic_area_mesh", "therapeutic_indication",
                                  "pharmacotherapeutic_group_human"))
        ema_cv = [w for w in CVWORDS if w in etxt.lower()]
        rows.append(dict(
            target=g, drug=drug.title(), agency="EMA",
            approval_record=("yes: " if ema_rec else "no approval: " if em else "no record")
            + ("; ".join("%s (%s)" % (e.get("name_of_medicine"), e.get("medicine_status"))
                         for e in em[:4])),
            cv_indication_chembl=chembl_cv,
            cv_indication_agency_text=(("yes: " + ", ".join(ema_cv[:4])) if ema_cv else
                                       "no EMA entry" if not em else
                                       "NO cardiovascular term found in the EMA fields"),
            agency_indication_excerpt=re.sub(r"\s+", " ", etxt)[:260],
            status=("authorised" if "Authorised" in sts else
                    "previously authorised (%s)" % ", ".join(s_.lower() for s_ in sts)
                    if ema_rec else "no approval record"),
            source_link="https://www.ema.europa.eu/en/medicines/download-medicine-data",
            label_link="; ".join(str(e.get("medicine_url", "")) for e in em[:2]),
            query_used=stem, accessed=ACCESSED))
        if fda_rec or ema_rec:
            T["drugs_with_record"].append(drug.title())
            if (fda_rec and fda_cv) or (ema_rec and ema_cv):
                T["drugs_text_confirmed"].append(drug.title())
        else:
            T["drugs_without_record"].append(drug.title())

# Evidence retrieved directly (2026-10-11), added to the rows it concerns. Nothing here changes
# a pass or a fail: APOB already passed on its Drugs@FDA record, and PLG still has none.
KYN = ("https://www.accessdata.fda.gov/drugsatfda_docs/appletter/2013/203568Orig1s000ltr.pdf",
       "https://www.accessdata.fda.gov/drugsatfda_docs/label/2013/203568s000lbl.pdf")
PB = open(os.path.join(HERE, "purple_book", "purplebook_latest_url.txt")).read().strip()
for r in rows:
    if (r["target"], r["drug"], r["agency"]) == ("APOB", "Mipomersen", "FDA"):
        r["approval_record"] += ("; original approval 2013-01-29, NDA 203568, new molecular "
                                 "entity (Drugs@FDA submission ORIG-1; approval letter read)")
        r["cv_indication_agency_text"] = "yes: cholesterol, lipid (approved label, 2013)"
        r["agency_indication_excerpt"] = (
            "KYNAMRO is an oligonucleotide inhibitor of apolipoprotein B-100 synthesis indicated "
            "as an adjunct to lipid-lowering medications and diet to reduce LDL-C, apo B, TC and "
            "non HDL-C in patients with homozygous familial hypercholesterolemia (HoFH)")
        r["label_link"] = "; ".join(KYN)
    if r["target"] == "PLG" and r["agency"] == "FDA":
        r["approval_record"] = (
            "no record: absent from Drugs@FDA; absent from the FDA Purple Book (September 2026 "
            "and the December 2020, 2021 and 2023 releases; searched for the proper and brand "
            "names); no DailyMed or openFDA label; no FDA approval or revocation notice in the "
            "Federal Register from 1994 onward")
        r["source_link"] += "; " + PB
    extra = {
        "Abciximab": ("yes: coronary, cardiac, angina (FDA label, ReoPro, 2019)",
                      "Abciximab is indicated as an adjunct to percutaneous coronary "
                      "intervention for the prevention of cardiac ischemic complications",
                      "https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/103575s5318lbl.pdf"),
        "Urokinase": ("yes: embol (FDA label, Kinlytic, 2008)",
                      "Kinlytic is indicated in adults for the lysis of acute massive pulmonary "
                      "emboli and of pulmonary emboli accompanied by unstable hemodynamics",
                      "https://www.accessdata.fda.gov/drugsatfda_docs/label/2008/021846s101lbl.pdf"),
        "Sildenafil": ("yes: hypertens (FDA label, Revatio)",
                       "REVATIO is indicated for the treatment of pulmonary arterial "
                       "hypertension (WHO Group I) in adults",
                       "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid="
                       "3bb9363e-b28d-4019-8aae-539233dca214"),
        "Tadalafil": ("yes: hypertens (FDA label, Adcirca)",
                      "ADCIRCA is a phosphodiesterase 5 inhibitor indicated for the treatment "
                      "of pulmonary arterial hypertension (WHO Group 1)",
                      "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid="
                      "ff61b237-be8e-461b-8114-78c52a8ad0ae"),
    }
    if r["agency"] == "FDA" and r["drug"] in extra and r["approval_record"].startswith("yes"):
        r["cv_indication_agency_text"], r["agency_indication_excerpt"], r["label_link"] = \
            extra[r["drug"]]
    if r["drug"] == "Abciximab" and r["agency"] == "FDA":
        r["status"] = "discontinued (FDA Purple Book: licence voluntarily revoked)"
    if r["target"] == "NR3C1":
        r["cv_indication_chembl"] += " [the cardiovascular indication relied on]"
    if r["target"] == "MS4A1":
        r["cv_indication_chembl"] += " [the cardiovascular indication relied on]"
for name, ph in ((os.path.join(MS, "APPROVED_35_CHECK.csv"), "Launched"),
                 (os.path.join(HERE, "approved_phase3_check.csv"), "Phase III")):
    with open(name, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows([r for r in rows if PHASE[r["target"]] == ph])
TXT = collections.defaultdict(set)      # after the direct evidence above has been added
for r in rows:
    if (r["approval_record"].startswith("yes")
            and r["cv_indication_agency_text"].startswith("yes")):
        TXT[r["target"]].add(r["drug"])
out = []
for g, T in by_t.items():
    out.append(dict(target=g, phase_primary=PHASE[g], eqtlgen_status=T["eqtlgen_status"],
                    strict_drugs=T["n_drugs"],
                    drugs_with_fda_or_ema_approval_record=len(T["drugs_with_record"]),
                    drugs_with_cv_term_in_agency_text=len(TXT[g]),
                    approved_fda_or_ema=bool(T["drugs_with_record"]),
                    approved_and_cv_text_confirmed=bool(TXT[g]),
                    drugs_with_record="; ".join(T["drugs_with_record"]),
                    drugs_without_record="; ".join(T["drugs_without_record"])))
with open(os.path.join(HERE, "approved_35_by_target.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
p3 = [o for o in out if o["phase_primary"] == "Phase III"]
print("[Q07] Phase III strict-verified targets: %s"
      % "; ".join("%s: %s" % (o["target"], "approval record" if o["approved_fda_or_ema"]
                              else "NO approval record (%s)" % o["drugs_without_record"])
                  for o in p3))
out_all = out
out = [o for o in out if o["phase_primary"] == "Launched"]
ok = [o for o in out if o["approved_fda_or_ema"]]
drop = [o for o in out if not o["approved_fda_or_ema"]]
print("[Q02] targets with an FDA or EMA approval record for >=1 strict-rule drug: %d of %d"
      % (len(ok), len(out)))
print("[Q03] targets that DROP OUT (no approval record for any strict-rule drug): %s"
      % ("; ".join("%s (%s; drugs: %s)" % (o["target"], o["eqtlgen_status"],
                                           o["drugs_without_record"]) for o in drop) or "none"))
tx = [o for o in ok if not o["approved_and_cv_text_confirmed"]]
print("[Q04] secondary check: passing targets where no agency indication text with a "
      "cardiovascular term was retrieved for any drug: %s"
      % ("; ".join("%s (%s)" % (o["target"], o["drugs_with_record"]) for o in tx) or "none"))
for lab, S in (("all 35", out), ("approved (FDA or EMA record)", ok)):
    c = collections.Counter(o["eqtlgen_status"] for o in S)
    print("[Q05] %-30s n=%d instrumented %d, no eQTL data %d, no qualifying instrument %d"
          % (lab, len(S), c["INSTRUMENTED"], c["NO_EQTL_DATA"], c["NO_QUALIFYING_INSTRUMENT"]))
print("[Q06] drugs without a record in either source, by target:")
for o in out:
    if o["drugs_without_record"]:
        print("        %-9s %d of %d drugs: %s" % (o["target"],
              len(o["drugs_without_record"].split("; ")), o["strict_drugs"],
              o["drugs_without_record"]))
