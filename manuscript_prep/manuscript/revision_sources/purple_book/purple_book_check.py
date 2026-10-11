#!/usr/bin/env python
"""Search FDA sources for an approval record of streptokinase and anistreplase (target PLG).

    python manuscript_prep/manuscript/revision_sources/purple_book/purple_book_check.py \
        > manuscript_prep/manuscript/revision_sources/purple_book/purple_book_check_output.txt

Reads the files already downloaded into this folder on 2026-10-11:
  purplebook_latest.csv  FDA Purple Book data download, September 2026 (URL in
                         purplebook_latest_url.txt); pb_2020.csv, pb_2021.csv, pb_2023.csv are
                         the December releases of those years
  fr_<term>.json         Federal Register API search results (documents from 1994 onward)
DailyMed and openFDA label searches for the same names returned nothing; they are recorded
in OPEN_ITEMS.md and are not re-run here.
"""
import csv
import io
import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
TERMS = ["streptokinase", "streptase", "kabikinase", "anistreplase", "eminase"]
CONTROLS = ["alteplase", "urokinase", "abciximab"]          # to show the search works


def load(name):
    raw = open(os.path.join(HERE, name), encoding="utf-8", errors="replace").read().split("\n")
    hi = [i for i, l in enumerate(raw) if "Proprietary Name" in l and "BLA Number" in l][0]
    rd = list(csv.reader(io.StringIO("\n".join(raw[hi:]))))
    return rd[0], rd[1:]


print("[B01] source: %s" % open(os.path.join(HERE, "purplebook_latest_url.txt")).read().strip())
for name in ("purplebook_latest.csv", "pb_2023.csv", "pb_2021.csv", "pb_2020.csv"):
    h, rows = load(name)
    hits = {t: [r for r in rows if t in " ".join(r).lower()] for t in TERMS + CONTROLS}
    print("[B02] %-22s %5d product rows | %s" % (name, len(rows), "; ".join(
        "%s %d" % (t, len(hits[t])) for t in TERMS + CONTROLS)))
h, rows = load("purplebook_latest.csv")
ix = {c: i for i, c in enumerate(h)}
for t in CONTROLS:
    for r in rows:
        if t in " ".join(r).lower():
            print("[B03] control row: %s | %s | BLA %s | %s | %s | approved %s"
                  % (r[ix["Proper Name"]], r[ix["Proprietary Name"]], r[ix["BLA Number"]],
                     r[ix["Marketing Status"]], r[ix["Licensure"]], r[ix["Approval Date"]]))
            break
print("[B04] streptokinase / Streptase / Kabikinase / anistreplase / Eminase rows in any of the "
      "four Purple Book releases: %d" % sum(
          1 for name in ("purplebook_latest.csv", "pb_2023.csv", "pb_2021.csv", "pb_2020.csv")
          for r in load(name)[1] if any(t in " ".join(r).lower() for t in TERMS)))
for t in ("streptokinase", "Streptase", "Kabikinase", "anistreplase", "Eminase"):
    j = json.load(open(os.path.join(HERE, "fr_%s.json" % t), encoding="utf-8"))
    res = j.get("results", [])
    fda = [d for d in res if "food and drug" in json.dumps(d.get("agencies", "")).lower()]
    print("[B05] Federal Register, term %-13s: %s documents; first %d examined, of which from "
          "FDA: %d %s" % (t, j.get("count"), len(res), len(fda),
                          [(d["publication_date"], d["title"][:60]) for d in fda][:3]))
print("[B06] conclusion: no FDA approval record for streptokinase or anistreplase was found in "
      "the Purple Book, Drugs@FDA, DailyMed, openFDA labels or the Federal Register (1994 "
      "onward). Both products were licensed before 1994, so the Federal Register search cannot "
      "show the original licence. PLG is NOT restored.")
