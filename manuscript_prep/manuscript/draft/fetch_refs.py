#!/usr/bin/env python
"""Resolve every cited DOI at Crossref and write the reference metadata.
    python manuscript_prep/manuscript/draft/fetch_refs.py
Writes draft/refs_crossref.json. A DOI that does not resolve is reported and not cited.
"""
import json, os, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
DOIS = """10.1038/s41586-024-07316-0 10.1038/ng.3314 10.1371/journal.pgen.1008489
10.1038/d41573-022-00120-3 10.1038/s41467-020-16969-0 10.12688/wellcomeopenres.16544.2
10.1038/s41588-020-0682-6 10.1038/s41588-021-00913-z 10.1038/s41467-021-26280-1
10.1016/j.xgen.2023.100341 10.1038/s41591-021-01310-z 10.1186/s12885-024-12449-6
10.1016/j.ebiom.2023.104859 10.1038/s41598-025-94761-0 10.1186/s12931-024-02848-5
10.3389/fcimb.2024.1386506 10.1126/science.aaz1776 10.1038/s41586-023-06592-6
10.1038/s41588-021-00978-w 10.1126/science.abj1541 10.1126/scitranslmed.aag1166
10.1038/s41588-023-01529-1 10.1016/j.tig.2020.08.009 10.1016/j.ajhg.2013.03.017
10.1186/s13073-026-01689-9 10.1016/j.xhgg.2025.100556 10.1001/jama.2021.18236
10.1093/nar/gkad1004 10.1126/science.1260419 10.1038/s43856-024-00625-5
10.1186/s12916-024-03700-9 10.1371/journal.pgen.1004383 10.1038/d41573-023-00158-x
10.1093/nar/gkac888 10.1101/2025.08.05.668745 10.12688/wellcomeopenres.15555.3 10.1002/gepi.21998 10.1371/journal.pgen.1003444
10.1056/NEJMoa054013 10.1056/NEJMoa1405386 10.1056/NEJMoa1307095 10.1056/NEJMoa1308027 10.1056/NEJMoa1612790
10.1056/NEJMoa1604304 10.1016/S0140-6736(12)60110-X 10.1056/NEJMoa0902604 10.1126/science.1081578
10.1161/01.CIR.87.6.1854""".split()
f = os.path.join(HERE, "refs_crossref.json")
have = json.load(open(f, encoding="utf-8")) if os.path.exists(f) else {}
for d in DOIS:
    if d in have:
        continue
    try:
        req = urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(d),
                                     headers={"User-Agent": "causal2drug-manuscript (mailto:none@example.org)"})
        m = json.load(urllib.request.urlopen(req, timeout=60))["message"]
        import re as _re
        _t = _re.sub(r"\s+", " ", _re.sub(r"<[^>]+>", " ", (m.get("title") or [""])[0])).strip()
        _t = _t.replace(" ,", ",").replace(" :", ":")
        have[d] = dict(title=_t,
                       journal=(m.get("short-container-title") or m.get("container-title") or [""])[0],
                       year=(m.get("published-print") or m.get("published") or m.get("issued"))["date-parts"][0][0],
                       volume=m.get("volume", ""), issue=m.get("issue", ""), page=m.get("page", "") or m.get("article-number", ""),
                       authors=[("%s %s" % (a.get("family", a.get("name", "")), "".join(w[0] for w in a.get("given", "").replace("-", " ").split()))).strip()
                                for a in m.get("author", [])])
    except Exception as e:
        print("NOT RESOLVED", d, e)
    time.sleep(0.3)
json.dump(have, open(f, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
for d in DOIS:
    r = have.get(d)
    if r:
        print("%s | %s | %s %s;%s:%s | %d authors: %s" % (d, r["title"][:90], r["journal"], r["year"], r["volume"], r["page"], len(r["authors"]), ", ".join(r["authors"][:3])))
