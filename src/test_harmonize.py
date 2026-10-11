#!/usr/bin/env python
"""Deliverable Part B2: deliberate synthetic test of the harmonization branches that
NEVER executed on real data.

In Deliverable 1 every one of 1,838 merged variants resolved `aligned_same`, because the
GWAS Catalog harmonised files and eQTL Catalogue both express alleles against the GRCh38
reference. The allele-flip, strand-flip, and strand-flip+allele-flip branches were
therefore completely untested. This constructs records that force each branch and checks
the SIGN of beta in every case.

Run:  python src/test_harmonize.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonize import harmonize_one  # noqa: E402

# Each case: exposure ref/alt, outcome effect/other, beta_out, expected status,
# expected harmonized beta, and a plain-language statement of what is being tested.
CASES = [
    dict(name="(0) control: no flip needed",
         exp=("G", "A"), out=("A", "G"), beta=0.30,
         exp_status="aligned_same", exp_beta=+0.30,
         why="outcome effect allele == exposure ALT; sign must be preserved"),

    dict(name="(a) ALLELE FLIP",
         exp=("G", "A"), out=("G", "A"), beta=0.30,
         exp_status="aligned_flipped", exp_beta=-0.30,
         why="outcome effect allele == exposure REF; sign MUST invert"),

    dict(name="(b) STRAND FLIP",
         exp=("G", "A"), out=("T", "C"), beta=0.30,
         exp_status="strand_flipped", exp_beta=+0.30,
         why="rc(T)=A == exposure ALT, rc(C)=G == exposure REF; sign preserved"),

    dict(name="(c) STRAND FLIP + ALLELE FLIP",
         exp=("G", "A"), out=("C", "T"), beta=0.30,
         exp_status="strand_flipped_and_allele_flipped", exp_beta=-0.30,
         why="rc(C)=G == exposure REF, rc(T)=A == exposure ALT; sign MUST invert"),

    dict(name="(d) genuine mismatch",
         exp=("G", "A"), out=("C", "G"), beta=0.30,
         exp_status="allele_mismatch", exp_beta=None,
         why="alleles cannot be reconciled by flip or strand; must be dropped"),

    dict(name="(e) palindromic, resolvable by frequency",
         exp=("A", "T"), out=("T", "A"), beta=0.30, af_exp=0.20, af_out=0.22,
         exp_status="aligned_same", exp_beta=+0.30,
         why="A/T palindrome, both AF far from 0.5 and concordant -> keep"),

    dict(name="(f) palindromic, AF too close to 0.5",
         exp=("C", "G"), out=("G", "C"), beta=0.30, af_exp=0.49, af_out=0.51,
         exp_status="ambiguous_unresolved", exp_beta=None,
         why="C/G palindrome, AF uninformative -> must be dropped, not guessed"),

    dict(name="(g) palindromic, AF discordant",
         exp=("A", "T"), out=("T", "A"), beta=0.30, af_exp=0.20, af_out=0.80,
         exp_status="ambiguous_unresolved", exp_beta=None,
         why="A/T palindrome, AFs on opposite sides of 0.5 -> unresolvable, drop"),

    dict(name="(h) palindromic allele-flip, frequency-resolved",
         exp=("A", "T"), out=("A", "T"), beta=0.30, af_exp=0.20, af_out=0.80,
         exp_status="aligned_flipped", exp_beta=-0.30,
         why="outcome effect == exposure REF; after flip AF_alt_out=0.20 matches 0.20"),

    dict(name="(i) indel must be rejected",
         exp=("G", "GA"), out=("GA", "G"), beta=0.30,
         exp_status="not_snp", exp_beta=None,
         why="multi-base allele; harmonization not attempted"),
]


def main():
    print("=" * 100)
    print("SYNTHETIC HARMONIZATION TEST -- exercising branches that never ran on real data")
    print("=" * 100)
    print("%-42s %-9s %-9s %-11s %-11s %s"
          % ("case", "exp", "out", "beta_in", "beta_out", "status"))
    print("-" * 100)
    n_pass = 0
    failures = []
    for c in CASES:
        r = harmonize_one(c["exp"][0], c["exp"][1], c["out"][0], c["out"][1], c["beta"],
                          af_alt_exp=c.get("af_exp"), af_effect_out=c.get("af_out"))
        got_status = r["status"]
        got_beta = r["beta_out_harmonized"]
        ok_status = got_status == c["exp_status"]
        ok_beta = (got_beta is None and c["exp_beta"] is None) or (
            got_beta is not None and c["exp_beta"] is not None
            and abs(got_beta - c["exp_beta"]) < 1e-12)
        ok = ok_status and ok_beta
        n_pass += ok
        if not ok:
            failures.append((c, got_status, got_beta))
        print("%-42s %-9s %-9s %-11s %-11s %s  %s"
              % (c["name"], "%s/%s" % c["exp"], "%s/%s" % c["out"],
                 "%+.2f" % c["beta"],
                 ("%+.2f" % got_beta) if got_beta is not None else "dropped",
                 got_status, "PASS" if ok else "*** FAIL ***"))
        print("%-42s expected: %-28s beta %s"
              % ("", c["exp_status"],
                 ("%+.2f" % c["exp_beta"]) if c["exp_beta"] is not None else "dropped"))
        print("%-42s why: %s" % ("", c["why"]))
    print("-" * 100)
    print("RESULT: %d/%d cases pass" % (n_pass, len(CASES)))
    if failures:
        print("\nFAILURES:")
        for c, gs, gb in failures:
            print("  %s: expected (%s, %s) got (%s, %s)"
                  % (c["name"], c["exp_status"], c["exp_beta"], gs, gb))
        return 1
    print("\nAll flip branches verified: sign inverts exactly when the outcome effect "
          "allele is the exposure REF allele, and is preserved otherwise, including "
          "under strand flip.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
