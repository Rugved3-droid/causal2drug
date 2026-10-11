#!/usr/bin/env python
"""Sensitivity result for the PCSK9 -> LDL cholesterol candidate-instrument control: the same
variant (rs11591147), with the exposure taken from the European-ancestry UKB-PPP discovery
cohort (sheet ST9, n = 34,557) instead of the combined cohort (ST10, n = 52,363) that the
protocol selects. Not a replacement for the protocol estimate.

    python manuscript_prep/manuscript/rescue_controls/sensitivity_pcsk9_discovery.py \
        > manuscript_prep/manuscript/rescue_controls/sensitivity_pcsk9_discovery_output.txt

Outcome row: the one already retrieved and harmonized by run_rescue_controls.py
(rescue_controls_results.json), unchanged.
"""
import json
import math
import os
import sys

import openpyxl

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
res = {r["control"]: r for r in json.load(open(os.path.join(HERE, "rescue_controls_results.json")))}
main = res["PCSK9_LDL"]
wb = openpyxl.load_workbook(os.path.join(MS, "recovery_data", "41586_2023_6592_MOESM3_ESM.xlsx"),
                            read_only=True, data_only=True)
it = wb["ST9"].iter_rows(values_only=True)
for row in it:
    h = [("" if c is None else str(c).strip()) for c in row]
    if "UKBPPP ProteinID" in h:
        break
hit = []
for row in it:
    r = dict(zip(h, row))
    if (str(r["UKBPPP ProteinID"]).split(":")[0] == "PCSK9" and str(r["cis/trans"]).lower() == "cis"
            and int(float(r["GENPOS (hg38)"])) == main["pos38"]):
        hit.append(r)
assert len(hit) == 1, len(hit)
r = hit[0]
vid = str(r[[k for k in r if k.startswith("Variant ID")][0]]).split(":")
assert (vid[2], vid[3]) == (main["other_allele"], main["effect_allele"])
bx, sx = float(r["BETA (discovery, wrt. A1)"]), float(r["SE (discovery)"])
by, sy = main["beta_outcome_harmonized"], main["se_outcome"]
b, se = by / bx, abs(sy / bx)
z = b / se
print("[D01] protocol estimate (combined cohort, ST10): exposure beta %+.5f (SE %.5f), F %.0f; "
      "MR beta %+.4f (95%% CI %+.4f to %+.4f)" % (main["beta_exposure"], main["se_exposure"],
                                                   main["F"], main["mr_beta"], main["mr_ci_low"],
                                                   main["mr_ci_high"]))
print("[D02] sensitivity (European discovery cohort, ST9, n = 34,557): %s allele %s, exposure "
      "beta %+.5f (SE %.5f), log10(p) %s, F %.0f" % (r["rsID"], vid[3], bx, sx,
                                                     r["log10(p) (discovery)"], (bx / sx) ** 2))
print("[D03] sensitivity MR beta %+.4f (95%% CI %+.4f to %+.4f), SE %.4f, z %.1f; outcome %s, "
      "beta %+.5f (SE %.5f) for the same allele; sign %s" % (
          b, b - 1.96 * se, b + 1.96 * se, se, z, main["outcome_study"], by, sy,
          "POSITIVE" if b > 0 else "NEGATIVE"))
print("[D04] difference, sensitivity minus protocol: %+.4f" % (b - main["mr_beta"]))
