#!/usr/bin/env python
"""Allele harmonization for cis-MR, extracted so it can be unit-tested directly.

Implements pre-registration R6. The exposure side is expressed as (ref, alt) and the
outcome side as (effect_allele, other_allele); the harmonized effect allele is always the
exposure ALT allele, and beta_outcome carries a sign flip when the outcome's effect allele
is the exposure REF allele.

The allele-flip and strand-flip branches never executed on real data in Deliverable 1
(every merge resolved aligned_same, because both sources express alleles against the same
GRCh38 reference). They are therefore exercised by an explicit synthetic test --
`python src/test_harmonize.py` -- rather than trusted untested.
"""
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def rc(allele):
    """Reverse complement."""
    return "".join(COMP.get(b, "N") for b in reversed(allele))


def is_snp(a):
    return len(a) == 1 and a in COMP


def harmonize_one(exp_ref, exp_alt, out_effect, out_other,
                  beta_out, af_alt_exp=None, af_effect_out=None):
    """Align one variant.

    Returns dict with:
      status : aligned_same | aligned_flipped | strand_flipped |
               strand_flipped_and_allele_flipped | allele_mismatch |
               ambiguous_unresolved | not_snp
      sign   : +1.0 / -1.0 / None
      beta_out_harmonized : beta_out * sign, expressed w.r.t. exp_alt
    """
    def fail(status):
        return dict(status=status, sign=None, beta_out_harmonized=None,
                    ambiguous=(exp_ref, exp_alt) in AMBIG)

    if not all(is_snp(x) for x in (exp_ref, exp_alt, out_effect, out_other)):
        return fail("not_snp")

    ambiguous = (exp_ref, exp_alt) in AMBIG

    if (out_effect, out_other) == (exp_alt, exp_ref):
        status, sign = "aligned_same", 1.0
    elif (out_effect, out_other) == (exp_ref, exp_alt):
        status, sign = "aligned_flipped", -1.0
    else:
        re_, ro_ = rc(out_effect), rc(out_other)
        if (re_, ro_) == (exp_alt, exp_ref):
            status, sign = "strand_flipped", 1.0
        elif (re_, ro_) == (exp_ref, exp_alt):
            status, sign = "strand_flipped_and_allele_flipped", -1.0
        else:
            return fail("allele_mismatch")

    # Strand-ambiguous SNPs are resolved by allele frequency, or dropped (R6).
    if ambiguous:
        if af_alt_exp is None or af_effect_out is None:
            return fail("ambiguous_unresolved")
        af_alt_out = af_effect_out if sign > 0 else 1.0 - af_effect_out
        informative = ((af_alt_exp < 0.42 or af_alt_exp > 0.58)
                       and (af_alt_out < 0.42 or af_alt_out > 0.58))
        concordant = (af_alt_exp < 0.5) == (af_alt_out < 0.5)
        if not (informative and concordant):
            return fail("ambiguous_unresolved")

    return dict(status=status, sign=sign, beta_out_harmonized=beta_out * sign,
                ambiguous=ambiguous)
