#!/usr/bin/env python
"""Fetch the current HGNC record, and the Ensembl GRCh38 record it points to, for every
in-scope target that is in no GTEx v8 egenes file under its GRCh37 Ensembl ID.
Responses are cached in recovery_out/id_change/ and read by build_recovery_map.py
(RECOVERY_RULE.md, amendment of 2026-10-11 no. 8).

    python -I manuscript_prep/manuscript/recovery_scripts/fetch_id_change_records.py SYMBOL ...
"""
import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "recovery_out", "id_change")
os.makedirs(OUT, exist_ok=True)


def get(url, name):
    f = os.path.join(OUT, name)
    if not os.path.exists(f):
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            open(f, "wb").write(r.read())
        time.sleep(0.4)
    return json.load(open(f, encoding="utf-8"))


for sym in sys.argv[1:]:
    docs = get("https://rest.genenames.org/fetch/symbol/" + sym,
               "hgnc_%s.json" % sym)["response"]["docs"]
    e = docs[0].get("ensembl_gene_id") if docs else None
    print(sym, docs[0].get("hgnc_id") if docs else None, e)
    if e:
        x = get("https://rest.ensembl.org/lookup/id/%s?content-type=application/json" % e,
                "ensembl38_%s.json" % e)
        print("   ", x.get("display_name"), x.get("seq_region_name"), x.get("start"),
              x.get("end"), x.get("strand"), x.get("assembly_name"))
