#!/usr/bin/env python
"""Validate the coloc ABF implementation against synthetic data with KNOWN truth.

This separates "the implementation is wrong" from "the data really says PP3". Each case
constructs summary statistics whose correct answer is known by construction, and checks
the posterior lands on the right hypothesis.
"""
import math
import random

W = 0.15 ** 2
P1, P2, P12 = 1e-4, 1e-4, 1e-5


def labf(beta, se, w=W):
    v = se * se
    r = w / (w + v)
    z = beta / se
    return 0.5 * (math.log(1 - r) + r * z * z)


def lse(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def ldiff(a, b):
    return a + math.log1p(-math.exp(b - a)) if b < a else float("-inf")


def coloc(pairs):
    l1 = [labf(b1, s1) for b1, s1, _, _ in pairs]
    l2 = [labf(b2, s2) for _, _, b2, s2 in pairs]
    s1, s2 = lse(l1), lse(l2)
    s12 = lse([a + b for a, b in zip(l1, l2)])
    ls = [0.0, math.log(P1) + s1, math.log(P2) + s2,
          math.log(P1) + math.log(P2) + ldiff(s1 + s2, s12), math.log(P12) + s12]
    d = lse(ls)
    return [math.exp(x - d) for x in ls]


def build(n, causal1, causal2, z1_peak, z2_peak, se1, se2, decay=0.06, seed=1):
    """One variant per position; z decays with distance from each trait's causal variant."""
    rnd = random.Random(seed)
    out = []
    for i in range(n):
        z1 = z1_peak * math.exp(-decay * abs(i - causal1)) if causal1 is not None else 0.0
        z2 = z2_peak * math.exp(-decay * abs(i - causal2)) if causal2 is not None else 0.0
        z1 += rnd.gauss(0, 0.3)
        z2 += rnd.gauss(0, 0.3)
        out.append((z1 * se1, se1, z2 * se2, se2))
    return out


CASES = [
    dict(name="A. SHARED causal variant, both strong",
         kw=dict(n=200, causal1=100, causal2=100, z1_peak=15, z2_peak=15,
                 se1=0.01, se2=0.01),
         expect="PP4", why="identical causal position -> H4 must dominate"),
    dict(name="B. DISTINCT causal variants, both strong",
         kw=dict(n=200, causal1=40, causal2=160, z1_peak=15, z2_peak=15,
                 se1=0.01, se2=0.01),
         expect="PP3", why="well separated causal positions -> H3 must dominate"),
    dict(name="C. trait 2 only (no eQTL signal)",
         kw=dict(n=200, causal1=None, causal2=100, z1_peak=0, z2_peak=15,
                 se1=0.01, se2=0.01),
         expect="PP2", why="only the outcome is associated -> H2"),
    dict(name="D. trait 1 only (no GWAS signal)",
         kw=dict(n=200, causal1=100, causal2=None, z1_peak=15, z2_peak=0,
                 se1=0.01, se2=0.01),
         expect="PP1", why="only the exposure is associated -> H1"),
    dict(name="E. neither associated",
         kw=dict(n=200, causal1=None, causal2=None, z1_peak=0, z2_peak=0,
                 se1=0.01, se2=0.01),
         expect="PP0", why="no signal anywhere -> H0"),
    dict(name="F. SHARED causal, hugely unequal precision",
         kw=dict(n=200, causal1=100, causal2=100, z1_peak=15, z2_peak=50,
                 se1=0.0085, se2=0.0014),
         expect="PP4",
         why="same position, GWAS 1.3M vs eQTL 32k -- unequal N must NOT break H4"),
]

NAMES = ["PP0", "PP1", "PP2", "PP3", "PP4"]


def main():
    print("=" * 104)
    print("SYNTHETIC COLOC VALIDATION -- known truth by construction")
    print("=" * 104)
    print("%-46s %9s %9s %9s %9s %9s  %-5s %s"
          % ("case", "PP0", "PP1", "PP2", "PP3", "PP4", "want", "result"))
    print("-" * 104)
    n_pass = 0
    for c in CASES:
        pp = coloc(build(**c["kw"]))
        top = NAMES[pp.index(max(pp))]
        ok = top == c["expect"]
        n_pass += ok
        print("%-46s %9.3g %9.3g %9.3g %9.3g %9.3g  %-5s %s"
              % (c["name"], pp[0], pp[1], pp[2], pp[3], pp[4], c["expect"],
                 "PASS" if ok else "*** FAIL ***"))
        print("%-46s %s" % ("", c["why"]))
    print("-" * 104)
    print("RESULT: %d/%d cases put the posterior on the correct hypothesis"
          % (n_pass, len(CASES)))
    return 0 if n_pass == len(CASES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
