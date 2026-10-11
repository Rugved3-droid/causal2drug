#!/usr/bin/env python
"""Build the clean public copy of the project that is pushed to GitHub.

    python manuscript_prep/make_public_copy.py D:\\causal2drug_publish

Copies code, rule documents and derived outputs. Leaves out third-party raw data, caches,
full-text copies of papers, the unpublished manuscript text and internal notes, and any file
over 20 MB. The working project in D:\\causal2drug is not modified.
"""
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = sys.argv[1]
MAXB = 20 * 1024 * 1024
MIT = """MIT License

Copyright (c) 2026 Rugved Parmar and Md Fahim

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
TOP = ["PROTOCOL.md", "DECISIONS.md", "INSTRUMENT_AVAILABILITY.md", "src", "deliverables",
       "audit", "manuscript_prep"]
SKIP_DIRS = {"__pycache__", ".git"}
SKIP_PATHS = [                                   # relative to ROOT, forward slashes
    "data", "logs", "manuscript_prep.zip",
    "manuscript_prep/verify/cache",
    "manuscript_prep/manuscript/recovery_data",
    "manuscript_prep/manuscript/prior_art_search",
    "manuscript_prep/manuscript/PRIOR_ART.md",
    "manuscript_prep/PRIOR_ART_CHECK.md", "manuscript_prep/JOURNAL_FIT.md",
    "manuscript_prep/OUTLINE.md",
    "manuscript_prep/manuscript/manuscript.md", "manuscript_prep/manuscript/manuscript.docx",
    "manuscript_prep/manuscript/manuscript_FINAL.docx",
    "manuscript_prep/manuscript/supplement.md", "manuscript_prep/manuscript/supplement.docx",
    "manuscript_prep/manuscript/supplement_FINAL.docx",
    "manuscript_prep/manuscript/cover_letter.md", "manuscript_prep/manuscript/cover_letter.docx",
    "manuscript_prep/manuscript/OPEN_ITEMS.md", "manuscript_prep/manuscript/WORD_COUNT.md",
    "manuscript_prep/manuscript/draft/manuscript_src.md",
    "manuscript_prep/manuscript/draft/supplement_src.md",
    "manuscript_prep/originals_before_final",
    "manuscript_prep/manuscript/revision_sources/zenodo_21444599.json",
    "manuscript_prep/manuscript/revision_sources/zenodo_21444599_preprint.pdf",
    "manuscript_prep/manuscript/revision_sources/zenodo_21444599_preprint.txt",
    "manuscript_prep/manuscript/revision_sources/drug_status/ema_medicines.json",
    "manuscript_prep/manuscript/revision_sources/site_of_action",
    "manuscript_prep/manuscript/il6r_audit/eqtlgen_IL6R_rows.tsv",
    "manuscript_prep/manuscript/recovery_out/gtex_target_signif_pairs.tsv",
    "manuscript_prep/manuscript/recovery_out/gtex_target_egenes.tsv",
    "manuscript_prep/manuscript/recovery_out/gencode_v26_gene_bodies_from_gtex_egenes.tsv",
    "audit/out/eqtlgen_target_rows.tsv", "audit/out/gtex_egenes_all.tsv",
    "audit/out/gtex_signif_targets.tsv", "audit/out/hgnc_complete_set.txt", "audit/out/pqtl",
    "audit/out/eqtlgen_tested_universe_audit.tsv",
]
SKIP_SUFFIX = ("~", ".pyc")
# Superseded working material goes under archive/ so that the top level shows the final
# analysis. Longest prefix wins. Paths containing "BUGGY" are dropped: they are truncated
# checkpoints of a clumping run that was later redone.
ARCHIVE = ["deliverables", "INSTRUMENT_AVAILABILITY.md",
           "manuscript_prep/HEADLINE.md", "manuscript_prep/NUMBERS.md",
           "manuscript_prep/INVENTORY.md", "manuscript_prep/DISCREPANCIES.md",
           "manuscript_prep/LIMITATIONS.md", "manuscript_prep/TABLE3_GROUPING_PROPOSAL.md",
           "manuscript_prep/table3_grouping_proposal.csv",
           "manuscript_prep/COLOC_CLAIM_CORRECTION.md",
           "manuscript_prep/manuscript/recovery_out/recovery_map_before_amendment8.csv",
           "manuscript_prep/manuscript/recovery_out/recovery_map_summary_before_amendment8.txt",
           "manuscript_prep/manuscript/recovery_out/recovery_map_before_audit_fix.csv",
           "manuscript_prep/manuscript/recovery_out/recovery_map_summary_before_audit_fix.txt",
           "manuscript_prep/manuscript/recovery_out/build_recovery_map_before_audit_fix.py"]
DROP_IF_CONTAINS = ("BUGGY", "originals_before_", "PUBLIC_README.md")


def dest_rel(r):
    for a in sorted(ARCHIVE, key=len, reverse=True):
        if r == a or r.startswith(a + "/"):
            return "archive/" + r
    return r
skip = {p.lower() for p in SKIP_PATHS}
n, nb, big, left = 0, 0, [], []
if os.path.exists(DEST):
    for name in os.listdir(DEST):
        if name == ".git":
            continue
        q = os.path.join(DEST, name)
        shutil.rmtree(q) if os.path.isdir(q) else os.remove(q)
os.makedirs(DEST, exist_ok=True)


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def want(p):
    r = rel(p).lower()
    return not any(r == s or r.startswith(s + "/") for s in skip)


for top in TOP:
    src = os.path.join(ROOT, top)
    if os.path.isfile(src):
        o_ = os.path.join(DEST, dest_rel(top))
        os.makedirs(os.path.dirname(o_), exist_ok=True)
        shutil.copy2(src, o_)
        n += 1
        continue
    for d, ds, fs in os.walk(src):
        ds[:] = [x for x in ds if x not in SKIP_DIRS and want(os.path.join(d, x))]
        for f in fs:
            p = os.path.join(d, f)
            if (not want(p) or f.endswith(SKIP_SUFFIX) or f.startswith("~$")
                    or any(k in rel(p) for k in DROP_IF_CONTAINS)):
                left.append(rel(p))
                continue
            sz = os.path.getsize(p)
            if sz > MAXB:
                big.append((rel(p), sz))
                continue
            out = os.path.join(DEST, dest_rel(rel(p)))
            os.makedirs(os.path.dirname(out), exist_ok=True)
            shutil.copy2(p, out)
            n += 1
            nb += sz
shutil.copy2(os.path.join(ROOT, "manuscript_prep", "PUBLIC_README.md"),
             os.path.join(DEST, "README.md"))
open(os.path.join(DEST, "LICENSE"), "w", encoding="utf-8", newline="\n").write(MIT)
print("copied %d files, %.1f MB" % (n, nb / 1e6))
print("left out because over 20 MB: %s" % [(a, round(b / 1e6)) for a, b in big])
print("left out by rule: %d files" % len(left))
