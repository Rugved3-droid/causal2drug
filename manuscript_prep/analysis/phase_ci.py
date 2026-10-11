#!/usr/bin/env python
"""95% confidence intervals for the share of uninstrumentable targets by clinical phase.

    python manuscript_prep/analysis/phase_ci.py > manuscript_prep/analysis/phase_ci_output.txt

Proportions: Wilson score interval. Difference between two independent proportions
(launched minus preclinical): Newcombe hybrid score interval (method 10, Newcombe 1998),
built from the two Wilson intervals. Counts are read from analysis/out/key_numbers.json
(primary set), and for the recovery map from manuscript/recovery_map.csv.
Descriptive; no hypothesis was pre-specified.
"""
import csv
import json
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
K = json.load(open(os.path.join(HERE, "out", "key_numbers.json")))
PH = ["Preclinical", "Phase I", "Phase II", "Phase III", "Launched"]
Z = 1.959963984540054


def wilson(x, n):
    p = x / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return p, c - h, c + h


def newcombe(x1, n1, x2, n2):
    p1, l1, u1 = wilson(x1, n1)
    p2, l2, u2 = wilson(x2, n2)
    d = p1 - p2
    return d, d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), \
        d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def block(tag, title, counts):
    print("[%s] %s" % (tag, title))
    for ph in PH:
        x, n = counts[ph]
        p, lo, hi = wilson(x, n)
        print("  %-12s %3d / %3d = %5.1f%% (95%% CI %4.1f%% to %4.1f%%)"
              % (ph, x, n, 100 * p, 100 * lo, 100 * hi))
    d, lo, hi = newcombe(*counts["Launched"], *counts["Preclinical"])
    print("  Launched minus Preclinical: %+.1f percentage points (95%% CI %+.1f to %+.1f)"
          % (100 * d, 100 * lo, 100 * hi))
    print()


bp = K["by_phase"]["PRIMARY"]
block("P01", "No eQTL data (lower bound), primary set, by highest phase",
      {ph: (bp[ph]["no_eqtl_data"], bp[ph]["n"]) for ph in PH})
block("P02", "Uninstrumentable, upper bound (no eQTL data + no qualifying instrument)",
      {ph: (bp[ph]["upper_bound"], bp[ph]["n"]) for ph in PH})
block("P03", "Autosomal-only lower bound",
      {ph: (bp[ph]["autosomal_lower_bound"], bp[ph]["n_not_on_x_y_mt"]) for ph in PH})
M = list(csv.DictReader(open(os.path.join(MP, "manuscript", "recovery_map.csv"),
                             encoding="utf-8")))
block("P04", "Locked, among targets with no eQTLGen whole-blood instrument (recovery map)",
      {ph: (sum(1 for x in M if x["phase_primary"] == ph
                and x["reporting_tier"].startswith("locked")),
            sum(1 for x in M if x["phase_primary"] == ph)) for ph in PH})
