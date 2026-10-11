#!/usr/bin/env python
"""Measure the total byte size of an Apache-indexed FTP-over-HTTPS directory
by summing HEAD Content-Length over every file (recursing one level into
subdirectories). Prints per-directory totals -- no size is estimated."""
import re
import sys
import concurrent.futures as cf
import requests

S = requests.Session()
LINK = re.compile(r'href="([^"?][^"]*)"')


def listing(url):
    r = S.get(url, timeout=120)
    r.raise_for_status()
    out = []
    for h in LINK.findall(r.text):
        if h.startswith("/") or h.startswith("http"):
            continue
        out.append(h)
    return out


def size_of(url):
    r = S.head(url, timeout=120, allow_redirects=True)
    if r.status_code != 200:
        return -1
    return int(r.headers.get("Content-Length", -1))


def walk(url, depth=1):
    files, dirs = [], []
    for h in listing(url):
        (dirs if h.endswith("/") else files).append(url + h)
    total, nfiles = 0, 0
    with cf.ThreadPoolExecutor(max_workers=10) as ex:
        for s in ex.map(size_of, files):
            if s > 0:
                total += s
                nfiles += 1
    if depth > 0:
        for d in dirs:
            t, n = walk(d, depth - 1)
            total += t
            nfiles += n
    return total, nfiles


if __name__ == "__main__":
    for u in sys.argv[1:]:
        t, n = walk(u if u.endswith("/") else u + "/")
        print("%-14s %15d bytes  %8.2f MiB  %6.2f GiB  n_files=%d\n    %s"
              % (u.rstrip("/").split("/")[-1], t, t / 2**20, t / 2**30, n, u))
