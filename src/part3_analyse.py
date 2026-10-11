#!/usr/bin/env python
"""PART 3 analysis (PREREGISTRATION R14.3-R14.5).

R14.3  does the abf-vs-susie PP4 discrepancy grow with eQTL lead z and with r2 < 1?
R14.4  head-to-head: coloc.abf PP4 with GTEx blood (n=853) vs eQTLGen (n=31,684),
       paired per target. THE PRE-REGISTERED REFUTATION TEST.
R14.5  distribution of the DENTIST residual-z RMSE across loci.
"""
import csv
import math

import numpy as np
from scipy import stats

OUT = r"D:\causal2drug\deliverables\d1h"
FLOOR = 1e-300


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def l10(x):
    return math.log10(max(x, FLOOR))


def main():
    rows = list(csv.DictReader(open(OUT + r"\part3_results.tsv", encoding="utf-8"),
                               delimiter="\t"))
    print("loci attempted: %d" % len(rows))
    notes = {}
    for r in rows:
        if r["note"]:
            notes[r["note"].split(":")[0]] = notes.get(r["note"].split(":")[0], 0) + 1
    if notes:
        print("notes:", notes)

    ok = [r for r in rows if r["abf_PP4_eqtlgen"] != ""]
    print("loci with a coloc result (eQTLGen exposure): %d" % len(ok))

    # ---------------- R14.5 DENTIST residual-z RMSE distribution ----------------
    print("\n" + "=" * 78)
    print("R14.5  DENTIST residual-z RMSE across loci  (Chen et al. 2021, Nat Commun 12:7117)")
    print("=" * 78)
    for lab, key in (("GWAS side", "dentist_rmse_gwas"),
                     ("eQTLGen side", "dentist_rmse_eqtlgen")):
        v = np.array([f(r[key]) for r in ok])
        v = v[~np.isnan(v)]
        if len(v) == 0:
            continue
        q = np.percentile(v, [25, 50, 75])
        print("  %-13s n=%-4d median=%-7.2f IQR=[%.2f, %.2f]  max=%-7.2f  %% > 2 : %.0f%%"
              % (lab, len(v), q[1], q[0], q[2], v.max(), 100 * (v > 2).mean()))
    print("  (~1 expected if the LD reference matched the study; >2 indicates misspecification)")

    # ---------------- R14.3 degradation: abf vs susie ----------------
    print("\n" + "=" * 78)
    print("R14.3  does the abf-vs-susie PP4 discrepancy grow with eQTL lead z, and with r2<1?")
    print("=" * 78)
    d = [r for r in ok if r["susie_PP4_eqtlgen"] != ""]
    print("  loci with BOTH abf and susie PP4: %d" % len(d))
    if len(d) >= 8:
        disc = np.array([l10(f(r["susie_PP4_eqtlgen"])) - l10(f(r["abf_PP4_eqtlgen"]))
                         for r in d])
        z = np.array([abs(f(r["eqtl_lead_z"])) for r in d])
        r2 = np.array([f(r["lead_r2"]) for r in d])
        gap = np.array([f(r["abf_gap"]) for r in d])
        for lab, x in (("eQTL lead |z|", z), ("lead-to-lead r2", r2), ("abf lBF gap", gap)):
            rho, p = stats.spearmanr(x, disc)
            print("  discrepancy log10(PP4_susie) - log10(PP4_abf)  vs  %-16s "
                  "Spearman rho=%+.3f  p=%.3g" % (lab, rho, p))
        print("  discrepancy: median=%+.2f  IQR=[%+.2f, %+.2f]"
              % (np.median(disc), *np.percentile(disc, [25, 75])))
    else:
        print("  too few loci with both methods; test not run")

    # ---------------- R14.4 HEAD-TO-HEAD, the refutation test ----------------
    print("\n" + "=" * 78)
    print("R14.4  HEAD-TO-HEAD: coloc.abf PP4, GTEx blood (n=853) vs eQTLGen (n=31,684)")
    print("       PRE-REGISTERED: if PP4 INCREASES with power, Class C is REFUTED.")
    print("=" * 78)
    h = [r for r in ok if r["abf_PP4_gtexblood"] != ""]
    print("  targets with BOTH exposures: %d" % len(h))
    if len(h) >= 6:
        a = np.array([f(r["abf_PP4_gtexblood"]) for r in h])   # low power
        b = np.array([f(r["abf_PP4_eqtlgen"]) for r in h])     # high power
        la = np.array([l10(x) for x in a])
        lb = np.array([l10(x) for x in b])
        lower = int((b < a).sum())
        higher = int((b > a).sum())
        print("  PP4 LOWER in the more powerful resource (eQTLGen) : %d / %d" % (lower, len(h)))
        print("  PP4 HIGHER in the more powerful resource          : %d / %d" % (higher, len(h)))
        bt = stats.binomtest(lower, lower + higher, 0.5, alternative="two-sided")
        print("  exact binomial (lower vs higher)                  : p = %.4g" % bt.pvalue)
        try:
            w = stats.wilcoxon(lb, la)
            print("  Wilcoxon signed-rank on log10(PP4)                : W=%.1f p=%.4g"
                  % (w.statistic, w.pvalue))
        except ValueError as e:
            print("  Wilcoxon not computable: %s" % e)
        print("  median log10 PP4: GTEx blood %+.2f   eQTLGen %+.2f   difference %+.2f"
              % (np.median(la), np.median(lb), np.median(lb - la)))
        print("  median PP4      : GTEx blood %.3g    eQTLGen %.3g"
              % (np.median(a), np.median(b)))
        print()
        if higher > lower and bt.pvalue < 0.05:
            print("  VERDICT: PP4 INCREASES with exposure power ->")
            print("           **CLASS C POWER-DEPENDENCE CLAIM IS REFUTED** (pre-registered).")
        elif lower > higher and bt.pvalue < 0.05:
            print("  VERDICT: PP4 DECREASES with exposure power -> consistent with Class C.")
        else:
            print("  VERDICT: NULL - no significant directional difference. Reported as a null.")
    else:
        print("  too few paired targets; test not run")

    print("\n" + "=" * 78)
    print("R14.5b  intermediate power point: UNAVAILABLE (see PREREGISTRATION R14.5).")
    print("  Non-monotonicity is UNVERIFIABLE with available data. Two power points give a")
    print("  slope, not a shape. No curve is claimed.")
    print("=" * 78)


if __name__ == "__main__":
    main()
