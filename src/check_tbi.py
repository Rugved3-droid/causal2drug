#!/usr/bin/env python
"""
HEAD-probe GWAS Catalog harmonised files for a sibling .tbi index.
Records HTTP status + Content-Length for both the data file and the index.
Deterministic sampling (fixed seed) so the sample is reproducible.
"""
import argparse
import concurrent.futures as cf
import csv
import hashlib
import random
import sys

import requests

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"


def probe(rel):
    url = BASE + rel.lstrip("./")
    out = {"rel": rel, "data_url": url}
    s = requests.Session()
    try:
        r = s.head(url, timeout=60, allow_redirects=True)
        out["data_status"] = r.status_code
        out["data_bytes"] = int(r.headers.get("Content-Length", -1))
    except Exception as e:
        out["data_status"] = "ERR:" + type(e).__name__
        out["data_bytes"] = -1
    try:
        r = s.head(url + ".tbi", timeout=60, allow_redirects=True)
        out["tbi_status"] = r.status_code
        out["tbi_bytes"] = int(r.headers.get("Content-Length", -1))
    except Exception as e:
        out["tbi_status"] = "ERR:" + type(e).__name__
        out["tbi_bytes"] = -1
    out["tbi_present"] = (out["tbi_status"] == 200)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=0, help="0 = all")
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()

    rels = [l.strip() for l in open(a.list) if l.strip()]
    print("population_size=%d" % len(rels))
    if a.n and a.n < len(rels):
        random.seed(a.seed)
        rels = random.sample(rels, a.n)
        print("sampled=%d seed=%d" % (len(rels), a.seed))

    rows = []
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(probe, rels), 1):
            rows.append(r)
            if i % 25 == 0:
                sys.stderr.write("  %d/%d\n" % (i, len(rels)))

    cols = ["rel", "data_url", "data_status", "data_bytes",
            "tbi_status", "tbi_bytes", "tbi_present"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    k = sum(1 for r in rows if r["tbi_present"])
    p = k / n if n else 0.0
    se = (p * (1 - p) / n) ** 0.5 if n else 0.0
    print("n_probed=%d  n_with_tbi=%d  proportion=%.4f" % (n, k, p))
    print("wald_95CI=[%.4f, %.4f]" % (max(0, p - 1.96 * se), min(1, p + 1.96 * se)))
    ok = [r["data_bytes"] for r in rows if r["data_bytes"] > 0]
    if ok:
        ok.sort()
        print("data_file_bytes: min=%d median=%d max=%d sum=%d"
              % (ok[0], ok[len(ok) // 2], ok[-1], sum(ok)))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
