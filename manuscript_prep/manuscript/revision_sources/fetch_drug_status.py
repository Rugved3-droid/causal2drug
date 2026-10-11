#!/usr/bin/env python
"""Regulatory status and labelled mechanism for every drug named in Table 2.

    python manuscript_prep/manuscript/revision_sources/fetch_drug_status.py

Sources (cached under revision_sources/drug_status/):
  FDA  openFDA Drugs@FDA endpoint  https://api.fda.gov/drug/drugsfda.json   (marketing status)
  FDA  openFDA drug label endpoint https://api.fda.gov/drug/label.json      (mechanism of action)
  EMA  medicines table (JSON report from ema.europa.eu), if retrievable
Nothing is inferred: a drug with no record in a source is written "no record".
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, "drug_status")
os.makedirs(C, exist_ok=True)
MS = os.path.dirname(HERE)
ACCESSED = "2026-10-11"


def get(url, name):
    f = os.path.join(C, name)
    if not os.path.exists(f):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                open(f, "wb").write(r.read())
        except urllib.error.HTTPError as e:
            open(f, "w").write(json.dumps({"error": {"code": e.code}}))
        except Exception as e:                                   # noqa: BLE001
            return {"error": {"code": "FETCH_FAILED", "detail": repr(e)[:200]}}
        time.sleep(0.4)
    try:
        return json.load(open(f, encoding="utf-8"))
    except Exception:                                            # noqa: BLE001
        return {"error": {"code": "NOT_JSON"}}


rows = list(csv.DictReader(open(os.path.join(MS, "display", "table2_named_targets.csv"),
                                encoding="utf-8-sig")))
pairs = []
for r in rows:
    tgt = r["Target"].replace(" (X)", "")
    for d in r["Approved cardiovascular drug(s) in ChEMBL"].split("; "):
        pairs.append((tgt, d.strip()))

ema = get("https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json",
          "ema_medicines.json")
ema_rows = ema.get("data", ema if isinstance(ema, list) else []) if isinstance(ema, (dict, list)) \
    else []
print("EMA medicines table rows:", len(ema_rows), "| keys:",
      list(ema_rows[0].keys())[:12] if ema_rows else ema.get("error") if isinstance(ema, dict)
      else "")

out = []
for tgt, drug in pairs:
    q = drug.upper()
    stem = q.replace(" MEDOXOMIL", "").replace(" CILEXETIL", "").replace(" ETEXILATE", "") \
            .replace(" SODIUM", "")
    key = stem.replace(" ", "_")
    u = "https://api.fda.gov/drug/drugsfda.json?" + urllib.parse.urlencode(
        {"search": 'products.active_ingredients.name:"%s"' % stem, "limit": 99})
    j = get(u, "drugsfda_%s.json" % key)
    st, apps = {}, []
    for res in j.get("results", []):
        for p in res.get("products", []):
            names = [a["name"].upper() for a in p.get("active_ingredients", [])]
            if any(stem in n for n in names):
                st[p.get("marketing_status", "?")] = st.get(p.get("marketing_status", "?"), 0) + 1
                apps.append("%s (%s)" % (res.get("application_number"), p.get("brand_name")))
    fda = ("no record in Drugs@FDA (openFDA)" if not st else
           "; ".join("%s: %d product(s)" % kv for kv in sorted(st.items())))
    fda_class = ("no FDA record" if not st else
                 "marketed" if any(k in st for k in ("Prescription", "Over-the-counter"))
                 else "discontinued" if "Discontinued" in st else "other")
    u2 = "https://api.fda.gov/drug/label.json?" + urllib.parse.urlencode(
        {"search": 'openfda.generic_name:"%s"' % stem, "limit": 1})
    lj = get(u2, "label_%s.json" % key)
    moa, setid, brand = "", "", ""
    if lj.get("results"):
        L = lj["results"][0]
        moa = " ".join(L.get("mechanism_of_action", L.get("clinical_pharmacology", [""])))[:1500]
        setid = L.get("set_id", "")
        brand = "; ".join(L.get("openfda", {}).get("brand_name", [])[:2])
    # matched on the EMA table's own substance fields only, not on free text
    em = [e for e in ema_rows
          if stem.lower() in (str(e.get("international_non_proprietary_name_common_name", ""))
                              + " | " + str(e.get("active_substance", ""))).lower()
          and str(e.get("category", "Human")).lower().startswith("human")]
    ema_txt = "; ".join(sorted({"%s: %s" % (e.get("name_of_medicine") or e.get("medicine_name")
                                           or "?", e.get("medicine_status")
                                           or e.get("authorisation_status") or "?")
                                for e in em})) or ("EMA table not retrieved" if not ema_rows
                                                   else "no record in the EMA medicines table")
    out.append(dict(target=tgt, drug=drug.title(), fda_marketing_status=fda,
                    fda_class=fda_class, fda_applications="; ".join(sorted(set(apps))[:4]),
                    fda_source=u, ema_status=ema_txt, label_brand=brand,
                    label_set_id=setid,
                    label_source=("https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid="
                                  + setid) if setid else "",
                    label_mechanism_text=moa, accessed=ACCESSED))
    print("%-8s %-22s FDA: %-48s | EMA: %s" % (tgt, drug[:22], fda[:48], ema_txt[:90]))
with open(os.path.join(HERE, "drug_status.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
