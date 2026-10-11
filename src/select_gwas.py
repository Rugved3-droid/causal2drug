#!/usr/bin/env python
"""Apply pre-registration rule R2: pick the outcome GWAS for a trait EFO.

Selection is mechanical -- largest parsed sample size among GWAS Catalog studies that
(a) map to the EFO, (b) have harmonised summary statistics, (c) have a sibling .tbi.
Prints the full ranked shortlist so the choice is auditable, not just the winner.
"""
import concurrent.futures as cf
import csv
import re
import sys

import requests

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
STUDIES = r"D:\causal2drug\data\raw\gwas-catalog-download-studies-v1.0.3.1.txt"
HARM = r"D:\causal2drug\data\raw\harmonised_list.txt"


def has_tbi(rel):
    url = BASE + rel.lstrip("./")
    try:
        r = requests.head(url + ".tbi", timeout=60, allow_redirects=True)
        return rel, r.status_code == 200, int(r.headers.get("Content-Length", -1))
    except Exception:
        return rel, False, -1


ANCESTRIES = ["African", "East Asian", "South Asian", "Hispanic", "Latin", "Native",
              "Asian unspecified", "Oceanian", "Greater Middle Eastern", "Aboriginal",
              "Sub-Saharan", "NR"]


def european_only(sample_str):
    """Pre-registration R2a: sample description mentions European ancestry and no other
    ancestry group, so the outcome GWAS is ancestry-matched to the EUR LD panel."""
    s = sample_str
    if "European" not in s:
        return False
    return not any(a in s for a in ANCESTRIES)


def main(efo, out_csv):
    # (a) studies mapped to the EFO
    meta = {}
    with open(STUDIES, encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        i_uri, i_acc = head.index("MAPPED_TRAIT_URI"), head.index("STUDY ACCESSION")
        i_auth, i_date = head.index("FIRST AUTHOR"), head.index("DATE")
        i_n, i_tr = head.index("INITIAL SAMPLE SIZE"), head.index("MAPPED_TRAIT")
        for f in rd:
            if len(f) <= i_uri:
                continue
            if efo in f[i_uri]:
                meta[f[i_acc]] = (f[i_auth], f[i_date], f[i_n], f[i_tr])
    print("(a) studies mapped to %s: %d" % (efo, len(meta)))

    # (b) with harmonised sumstats
    rels = {}
    for l in open(HARM):
        l = l.strip()
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m and m.group(1) in meta:
            rels[m.group(1)] = l
    print("(b) of those, with harmonised sumstats: %d" % len(rels))

    # (c) with a .tbi
    res = {}
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        for rel, ok, nb in ex.map(has_tbi, rels.values()):
            res[rel] = (ok, nb)
    rows = []
    for acc, rel in rels.items():
        ok, nb = res[rel]
        auth, date, nstr, trait = meta[acc]
        n = max([int(x.replace(",", "")) for x in re.findall(r"[\d,]{4,}", nstr)] or [0])
        rows.append({"accession": acc, "max_parsed_N": n, "tbi_present": ok,
                     "european_only": european_only(nstr),
                     "tbi_bytes": nb, "first_author": auth, "date": date,
                     "mapped_trait": trait, "initial_sample_size": nstr, "rel_path": rel})
    n_tbi = sum(1 for r in rows if r["tbi_present"])
    print("(c) of those, with .tbi: %d" % n_tbi)

    rows.sort(key=lambda r: (-r["max_parsed_N"], r["accession"]))
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_eur = sum(1 for r in rows if r["european_only"])
    print("(d) of those, European-ancestry-only (R2a): %d" % n_eur)
    elig = [r for r in rows if r["tbi_present"] and r["european_only"]]
    print("eligible after all of (a)-(d): %d" % len(elig))
    print("\nTop 5 eligible (tbi_present=True), ranked by max_parsed_N:")
    for r in elig[:5]:
        print("  %-14s N=%-9d %-14s %-11s %s" %
              (r["accession"], r["max_parsed_N"], r["first_author"], r["date"],
               r["mapped_trait"][:40]))
    print("\nSELECTED: %s  (N=%d)  %s" %
          (elig[0]["accession"], elig[0]["max_parsed_N"], elig[0]["rel_path"]))
    print("wrote", out_csv)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
