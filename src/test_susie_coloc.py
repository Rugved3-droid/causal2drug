#!/usr/bin/env python
"""Validate SuSiE-RSS + coloc.susie against synthetic data with KNOWN truth.

IMPORTANT CORRECTION, recorded here because it changes what the method can be claimed to
do. The Part A report recommended coloc.susie as the remedy for the HMGCR failure mode
(lead-variant displacement inside one LD block). That recommendation was imprecise.

For any single-effect Bayes factor,  lbf_j = -0.5*log(1+v/s2) + 0.5*(z_j^2/s2)*(v/s2)/(1+v/s2),
so the lbf DIFFERENCE between two variants scales with (z_j^2 - z_k^2)/s2. The residual
variance s2 is what flattens it: with s2 estimated, displacement is partially rescued on
synthetic data (D1 below). With s2 fixed at 1 -- a defect in the first draft -- coloc.susie
inherits coloc.abf's full (1/2)z^2 peakiness and rescues nothing.

What coloc.susie does fix is MULTIPLE CAUSAL VARIANTS: where one trait carries two signals
and the other one, coloc.abf blends them and returns H3, while coloc.susie separates the
components and finds the shared one. Cases M1/M2 test exactly that; D1 documents the
displacement limitation honestly rather than hiding it.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from susie_coloc import coloc_susie, coloc_from_lbf, ld_condition_number  # noqa: E402

W = 0.15 ** 2
NAMES = ["PP0", "PP1", "PP2", "PP3", "PP4"]


def abf_coloc(z1, z2, se1, se2):
    def labf(z, se):
        v = se * se
        r = W / (W + v)
        return 0.5 * (math.log(1 - r) + r * z * z)
    return coloc_from_lbf(np.array([labf(z, se1) for z in z1]),
                          np.array([labf(z, se2) for z in z2]))


def ar1(p, rho):
    i = np.arange(p)
    return rho ** np.abs(i[:, None] - i[None, :])


def run(name, z1, z2, R, want, note):
    a = abf_coloc(z1, z2, 0.0085, 0.0014)
    # strict=False: these are GROUND-TRUTH synthetic cases where a shared causal variant
    # is the correct answer, and the FIX-1 assertion (condition (a): both sides top at the
    # same variant) would reject exactly that. See test_coloc_assertion.py, which measures
    # the cost of the assertion explicitly.
    s = coloc_susie(z1, z2, R, L=3, strict=False)
    b = s["best"]
    got = NAMES[int(np.argmax([b[k] for k in NAMES]))] if b else "NO_CS"
    ok = got == want
    print("%-46s %-7s %-11.3g %-11s %-7s %s"
          % (name[:46], want, a["PP4"], ("%.3g" % b["PP4"]) if b else "no CS",
             "%d/%d" % (s["n_cs1"], s["n_cs2"]), "PASS" if ok else "*** FAIL ***"))
    print("%-46s   %s" % ("", note))
    return ok


def main():
    p, c = 101, 50
    rng = np.random.default_rng(11)
    print("=" * 108)
    print("SuSiE-RSS + coloc.susie VALIDATION")
    print("=" * 108)
    print("%-46s %-7s %-11s %-11s %-7s %s"
          % ("case", "want", "abf PP4", "susie PP4", "nCS", "result"))
    print("-" * 108)
    res = []

    R = ar1(p, 0.90)
    z1 = 15 * R[c] + rng.normal(0, 1, p)
    z2 = 15 * R[c] + rng.normal(0, 1, p)
    res.append(run("B1 shared single causal variant", z1, z2, R, "PP4",
                   "both traits one signal at the same variant -> H4"))

    R2 = ar1(p, 0.70)
    z1 = 15 * R2[20] + rng.normal(0, 1, p)
    z2 = 15 * R2[80] + rng.normal(0, 1, p)
    res.append(run("B2 genuinely distinct causal variants", z1, z2, R2, "PP3",
                   "well-separated causal variants -> H3"))

    Rm = ar1(p, 0.80)
    z1 = 18 * Rm[20] + 18 * Rm[75] + rng.normal(0, 1, p)
    z2 = 30 * Rm[75] + rng.normal(0, 1, p)
    res.append(run("M1 two signals vs one shared  <-- susie's purpose", z1, z2, Rm,
                   "PP4", "abf blends trait1's two peaks; susie should separate them"))

    z1 = 18 * Rm[20] + 18 * Rm[75] + rng.normal(0, 1, p)
    z2 = 30 * Rm[45] + rng.normal(0, 1, p)
    res.append(run("M2 two signals, neither shared", z1, z2, Rm, "PP3",
                   "trait2's signal matches neither of trait1's -> H3"))

    print("-" * 108)
    print("RESULT: %d/%d" % (sum(res), len(res)))

    print()
    print("D1  DOCUMENTED LIMITATION: shared causal variant, observed peaks displaced")
    print("    (20 trials; both methods reported, no pass criterion attached)")
    Rd = ar1(p, 0.95)
    rng2 = np.random.default_rng(3)
    ab, su, ndisp = [], [], 0
    for _ in range(20):
        z1 = 15 * Rd[c] + rng2.normal(0, 1, p)
        z2 = 50 * Rd[c] + rng2.normal(0, 1, p)
        if int(np.argmax(np.abs(z1))) != int(np.argmax(np.abs(z2))):
            ndisp += 1
            a = abf_coloc(z1, z2, 0.0085, 0.0014)
            s = coloc_susie(z1, z2, Rd, L=3, strict=False)
            ab.append(a["PP4"])
            su.append(s["best"]["PP4"] if s["best"] else float("nan"))
    print("    trials with displaced peaks: %d/20" % ndisp)
    if ab:
        print("    median PP4  coloc.abf = %.3g   coloc.susie = %.3g"
              % (np.median(ab), np.nanmedian(su)))
    print("    -> with residual-variance estimation, credible-set coloc DOES partially")
    print("       rescue lead displacement on synthetic data (abf ~0.03 -> susie ~0.57).")
    print("       It does NOT rescue it on the real HMGCR locus; see deliverables/d1g.")
    print()
    print("LD conditioning (rho=0.95 block):", ld_condition_number(ar1(p, 0.95)))
    return 0 if all(res) else 1


if __name__ == "__main__":
    raise SystemExit(main())
