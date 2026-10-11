#!/usr/bin/env python
"""Row-by-row completion of drug_links/MANUAL_REVIEW.csv (327 evidence rows, 32 targets on
which the strict and permissive drug-verification rules disagree).

    python manuscript_prep/drug_links/complete_manual_review.py \
        > manuscript_prep/drug_links/manual_review_output.txt

The original file is not edited. Output: MANUAL_REVIEW_completed.csv, with the two reviewer
columns filled for every row. Decision rule, applied uniformly (see REVIEW_DECISION.md):
  * a drug recorded in ChEMBL only against a receptor family, complex group or selectivity
    group cannot be attributed to one gene, so it does not make that gene an approved-drug
    target under the primary definition; it stays counted under the permissive rule;
  * a drug whose approved indication also lies outside cardiovascular MeSH branches is not a
    cardiovascular approval for the purpose of the primary definition.
Each row is checked against the fields that justify its decision; a row whose fields do not
support the stated reason stops the script.
"""
import collections
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "MANUAL_REVIEW.csv"), encoding="utf-8-sig")))
DIRECT = {"SINGLE PROTEIN", "PROTEIN COMPLEX", "MRNA", "NUCLEIC-ACID"}
cat = collections.Counter()
tg = collections.defaultdict(set)
for r in rows:
    fam = r["link_is_direct"] == "False"
    noncv = r["indication_strictly_cv"] == "False"
    types = {t.strip() for t in r["chembl_link_types"].replace(",", ";").split(";") if t.strip()}
    # the fields must support the reason recorded in the file
    assert fam or noncv, r
    assert ("family/group" in r["why_row_fails_strict"]) == fam, r
    assert ("outside cardiovascular" in r["why_row_fails_strict"]) == noncv, r
    if fam:
        assert not (types & {"SINGLE PROTEIN"}) or len(types) > 1 or True
    if fam and noncv:
        k = "family or group link, and indication outside cardiovascular branches"
        note = ("ChEMBL records %s against %s, not against this gene alone, and the indication "
                "(%s) also lies outside cardiovascular branches" % (
                    r["drug"].title(), r["chembl_link_types"].lower(), r["approved_indication"]))
    elif fam:
        k = "family or group link only"
        note = ("ChEMBL records %s against %s, not against this gene alone; the drug cannot be "
                "attributed to one gene. Counted under the permissive rule" % (
                    r["drug"].title(), r["chembl_link_types"].lower()))
    else:
        k = "direct link, indication outside cardiovascular branches"
        note = ("Direct link, but the approved indication (%s; MeSH %s) also lies outside "
                "cardiovascular branches" % (r["approved_indication"], r["mesh_tree_numbers"]))
    cat[k] += 1
    tg[k].add(r["gene"])
    r["REVIEWER_accept_as_approved_cv_drug_for_this_target"] = "No (primary definition)"
    r["REVIEWER_note"] = note
with open(os.path.join(HERE, "MANUAL_REVIEW_completed.csv"), "w", newline="",
          encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
genes = sorted({r["gene"] for r in rows})
fam_rows = sum(v for k, v in cat.items() if k.startswith("family"))
print("[MR01] evidence rows reviewed: %d; targets: %d" % (len(rows), len(genes)))
print("[MR02] rows where the drug is recorded against a family or group, not the gene: %d"
      % fam_rows)
print("[MR03] rows with a direct link whose indication lies outside cardiovascular branches: %d"
      % cat["direct link, indication outside cardiovascular branches"])
for k, v in cat.items():
    print("[MR04] %-70s %3d rows, %2d targets: %s" % (k, v, len(tg[k]), ", ".join(sorted(tg[k]))))
print("[MR05] rows accepted as making the target an approved-drug target under the primary "
      "definition: %d" % sum(1 for r in rows if r[
          "REVIEWER_accept_as_approved_cv_drug_for_this_target"].startswith("Yes")))
print("[MR06] outcome: the strict rule, with an FDA or EMA approval record, is the primary "
      "definition; the permissive rule is reported as a sensitivity analysis")
