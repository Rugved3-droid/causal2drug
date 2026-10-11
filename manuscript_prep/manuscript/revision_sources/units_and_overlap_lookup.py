#!/usr/bin/env python
"""Keyword-in-context extracts from the open full text of the source papers, to establish the
scale of each exposure and outcome and which cohorts each study used. Read-only.

    python manuscript_prep/manuscript/revision_sources/units_and_overlap_lookup.py \
        > manuscript_prep/manuscript/revision_sources/units_and_overlap_extracts.txt
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "prior_art_search"))
import papers                                                    # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
Q = [
 ("Graham 2021 (lipids; GCST90239658, GCST90239664)", "34887591[pmid]",
  [r"inverse.normal", r"natural log", r"UK Biobank", r"Fenland", r"mg/dL|mmol", r"residuals"]),
 ("Keaton 2024 (blood pressure; GCST90310294)", "38689001[pmid]",
  [r"mm ?Hg", r"UK Biobank", r"antihypertensive|medication", r"inverse.normal", r"ICBP|Million Veteran"]),
 ("Koskeridis 2022 (CRP; GCST90179146)", "36376304[pmid]",
  [r"MTAG", r"natural log|log.transformed|ln\(", r"UK Biobank", r"CHARGE"]),
 ("Sun 2023 (UKB-PPP)", "37794186[pmid]",
  [r"inverse.rank|inverse.normal|rank.based", r"34,557|discovery cohort", r"NPX"]),
 ("Pietzner 2021 (Fenland)", "34648354[pmid]",
  [r"inverse.normal|rank.based|rank.inverse", r"residual"]),
 ("GTEx 2020", "32913098[pmid]", [r"inverse normal|normalized expression|TMM"]),
 ("Vosa 2021 (eQTLGen)", "34475573[pmid]", [r"UK Biobank", r"Z.score|z-score"]),
]
for name, term, pats in Q:
    r = papers.lookup(term, 1)
    if not r:
        print("### %s: not found in PubMed" % name)
        continue
    r = r[0]
    print("\n### %s | PMID %s | %s | DOI %s" % (name, r["pmid"], r["pmc"] or "no PMC id", r["doi"]))
    try:
        txt = papers.fulltext(r["pmc"]) if r["pmc"] else ""
    except Exception as e:                                       # noqa: BLE001
        txt = ""
        print("    full text not retrievable: %r" % (e,))
    print("    body characters: %d; abstract: %s" % (len(txt), r["abstract"][:400]))
    for pat in pats:
        hits = list(re.finditer(pat, txt, flags=re.I))
        print("  -- /%s/ : %d hit(s)" % (pat, len(hits)))
        for m in hits[:3]:
            a, b = max(0, m.start() - 230), min(len(txt), m.end() + 230)
            print("       ..." + txt[a:b].replace("\n", " ") + "...")
