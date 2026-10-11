#!/usr/bin/env python
"""Compare manuscript/recovery_map.csv, after the fixes of 2026-10-11, with the independent
audit's own output (audit/out/audit_recovery_map.csv), target by target.

    python manuscript_prep/manuscript/recovery_scripts/compare_with_audit.py \
        > manuscript_prep/manuscript/recovery_out/compare_with_audit_output.txt

Also compares with the map as it stood before the fixes
(recovery_out/recovery_map_before_audit_fix.csv), to list exactly what the fixes changed.
"""
import csv
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(MS))


def load(p, key):
    return {r[key]: r for r in csv.DictReader(open(p, encoding="utf-8"))}


new = load(os.path.join(MS, "recovery_map.csv"), "gene")
old = load(os.path.join(MS, "recovery_out", "recovery_map_before_audit_fix.csv"), "gene")
aud = load(os.path.join(ROOT, "audit", "out", "audit_recovery_map.csv"), "symbol")
print("[X01] targets: map %d, audit %d, same set: %s" % (len(new), len(aud), set(new) == set(aud)))


AUD = {"plausible " + "rescue": "candidate: CV-relevant tissue or plasma protein",
       "off-target-" + "tissue rescue": "candidate: other tissue only",
       "fixable" + "-tissue": "candidate-tissue", "fixable" + "-protein": "candidate-protein",
       "fixable" + "-both": "candidate-both"}
for _r in aud.values():           # the audit used the first tier and category names
    _r["tier"] = AUD.get(_r["tier"], _r["tier"])
    _r["category"] = AUD.get(_r["category"], _r["category"])
if "reporting_tier" in next(iter(old.values())):
    for _r in old.values():
        _r["reporting_tier"] = AUD.get(_r["reporting_tier"], _r["reporting_tier"])
        _r["category"] = AUD.get(_r["category"], _r["category"])


def rs(v):
    return v.split(" (")[0].strip()


def pos(v):
    """position from the map's 'rsID (chrN:POS)' label"""
    m = re.search(r"\(chr[^:]+:(\d+)\)", v)
    return m.group(1) if m else ""


def b(v):
    return str(v).strip().lower() == "true"


def diff(label, f):
    d = sorted(g for g in new if not f(new[g], aud[g]))
    print("[X02] %-58s differences: %d %s" % (label, len(d), ", ".join(d[:20])))
    return d


diff("category", lambda n, a: n["category"] == a["category"])
diff("reporting tier", lambda n, a: n["reporting_tier"] == a["tier"])
diff("GTEx ID used", lambda n, a: (n["gtex_ensembl_id"] or "") == (a["gtex_id_used"]
                                   if int(a["gtex_n_tissues_tested"] or 0) else ""))
diff("GTEx tissues tested (n)", lambda n, a: n["gtex_n_tissues_tested"] == a["gtex_n_tissues_tested"])
diff("GTEx non-blood tissues with a candidate",
     lambda n, a: set(filter(None, n["gtex_tissues_with_instrument_excl_blood"].split("; ")))
     == set(filter(None, a["gtex_nonblood_tissues"].split(";"))))
diff("GTEx best CV tissue and variant",
     lambda n, a: (n["gtex_best_cv_tissue"], n["gtex_best_cv_variant_b38"])
     == (a["gtex_best_cv_tissue"], a["gtex_best_cv_variant"]))
diff("UKB-PPP assayed", lambda n, a: b(n["ukbppp_assayed"]) == b(a["ukb_assayed"]))
diff("UKB-PPP candidate yes/no", lambda n, a: b(n["ukbppp_instrument"]) == b(a["ukb_instrument"]))
# the audit writes the UKB-PPP lead as the table's variant ID plus a GRCh38 position, and the
# Fenland lead as chr:pos_alleles; the map writes rsID (chr:pos). Compared on position.
diff("UKB-PPP lead variant (GRCh38 position)", lambda n, a: not b(n["ukbppp_instrument"])
     or pos(n["ukbppp_lead_variant"]) == a["ukb_pos38"])
diff("UKB-PPP sheet of the lead", lambda n, a: not b(n["ukbppp_instrument"])
     or n["ukbppp_note"] == a["ukb_sheet"])
diff("UKB-PPP -log10 p of the lead", lambda n, a: not b(n["ukbppp_instrument"])
     or abs(float(n["ukbppp_neglog10_p"]) - float(a["ukb_log10p"])) < 0.02)
diff("deCODE candidate yes/no", lambda n, a: b(n["decode_instrument"]) == b(a["decode_instrument"]))
diff("deCODE lead variant", lambda n, a: not b(n["decode_instrument"])
     or rs(n["decode_lead_variant"]) == rs(a["decode_lead"]))
diff("deCODE assayed (map) vs included (audit)",
     lambda n, a: b(n["decode_assayed"]) == b(a["decode_included"]))
diff("Fenland assayed", lambda n, a: b(n["fenland_assayed"]) == b(a["fenland_assayed"]))
diff("Fenland candidate yes/no", lambda n, a: b(n["fenland_instrument"]) == b(a["fenland_instrument"]))
diff("Fenland lead variant (GRCh37 position)", lambda n, a: not b(n["fenland_instrument"])
     or pos(n["fenland_lead_variant"]) == a["fenland_lead"].split(":")[1].split("_")[0])

print()
for col, lab in (("category", "category"), ("reporting_tier", "tier")):
    d = sorted(g for g in new if new[g][col] != old[g][col])
    print("[X03] %s changed by the fixes: %d %s" % (lab, len(d), d))
for st in ("ukbppp", "decode", "fenland"):
    d = sorted(g for g in new if rs(new[g][st + "_lead_variant"]) != rs(old[g][st + "_lead_variant"]))
    print("[X04] %s lead variant changed by the fixes: %d" % (st, len(d)))
    for g in d:
        print("        %-9s %s (p %s) -> %s (p %s)" % (g, rs(old[g][st + "_lead_variant"]),
              old[g][st + "_p"], rs(new[g][st + "_lead_variant"]), new[g][st + "_p"]))
d = sorted(g for g in new if new[g]["ukbppp_note"] != old[g]["ukbppp_note"])
print("[X05] UKB-PPP sheet of the lead changed by the fixes: %d %s" % (len(d), ", ".join(d)))
d = sorted(g for g in new if new[g]["fenland_assayed"] != old[g]["fenland_assayed"])
print("[X06] Fenland assayed changed by the fixes: %d %s" % (len(d), ", ".join(d)))
for st in ("ukbppp", "decode", "fenland"):
    d = sorted(g for g in new if new[g][st + "_instrument"] != old[g][st + "_instrument"])
    print("[X07] %s candidate yes/no changed by the fixes: %d %s" % (st, len(d), d))
