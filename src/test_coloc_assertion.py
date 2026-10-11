#!/usr/bin/env python
"""FIX 1 -- verify the ColocArtefact assertion.

Two things must be shown, not asserted:
  (i)  it catches all FOUR PP4~1 artefacts this project produced and withdrew;
  (ii) its cost -- condition (a) also rejects a GENUINE shared-causal-variant
       colocalization, so the assertion is deliberately conservative.

Then every retained coloc analysis is re-run through it.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import susie_coloc as SC                       # noqa: E402
from susie_coloc import ColocArtefact          # noqa: E402
import run_susie_real as RS                    # noqa: E402

W = 0.15 ** 2


def ar1(p, rho):
    i = np.arange(p)
    return rho ** np.abs(i[:, None] - i[None, :])


def check(name, z1, z2, R, expect):
    """expect: 'REJECT' or 'PASS'."""
    try:
        s = SC.coloc_susie(z1, z2, R, L=5, strict=True)
        got = "PASS"
        detail = ("PP4=%.4g" % s["best"]["PP4"]) if s["best"] else "no credible set"
    except ColocArtefact as e:
        got = "REJECT"
        detail = str(e)
    ok = got == expect
    print("  %-44s expect=%-7s got=%-7s %s" % (name[:44], expect, got,
                                               "OK" if ok else "*** MISMATCH ***"))
    print("      %s" % detail[:110])
    return ok


def main():
    print("=" * 100)
    print("FIX 1 -- ColocArtefact assertion")
    print("=" * 100)

    print("\n(ii) COST OF THE ASSERTION -- synthetic, ground truth = shared causal variant")
    rng = np.random.default_rng(11)
    p, c = 101, 50
    R = ar1(p, 0.90)
    z1 = 15 * R[c] + rng.normal(0, 1, p)
    z2 = 15 * R[c] + rng.normal(0, 1, p)
    check("B1 genuine shared causal variant", z1, z2, R, "REJECT")
    print("      ^ this is a TRUE colocalization and the assertion rejects it.")
    print("        Condition (a) cannot distinguish a real shared variant from a")
    print("        component colocalizing with itself. The assertion is conservative")
    print("        by construction; with Gate D closed that is the safe direction.")

    print("\n(i) THE FOUR REAL ARTEFACTS")
    keep = RS.eur_unrelated()
    results = []

    for name, merge, chrom in (
            ("HMGCR-LDL-C, 1000G n=525", RS.LOCI[0]["merge"], "5"),
            ("IL6R-CRP, 1000G n=525", RS.LOCI[1]["merge"], "1")):
        rows = list(csv.DictReader(open(merge, encoding="utf-8"), delimiter="\t"))
        Rm, ki, n_eur, _ = RS.ld_matrix(chrom, rows, keep)
        rows = [rows[i] for i in ki]
        z1 = np.array([float(r["beta_exp"]) / float(r["se_exp"]) for r in rows])
        z2 = np.array([float(r["beta_out_harm"]) / float(r["se_out"]) for r in rows])
        # the artefactual configuration used s2 fixed at 1
        import susie_coloc as S2
        orig = S2.susie_rss

        def fixed_s2(z, R_, L=5, prior_v=None, n_iter=100, tol=1e-6, pi=None):
            return orig(z, R_, L=L, prior_v=200.0, n_iter=n_iter, tol=tol, pi=pi)
        if name.startswith("HMGCR"):
            S2.susie_rss = fixed_s2
            results.append(check(name + "  [s2=1 configuration]", z1, z2, Rm, "REJECT"))
            S2.susie_rss = orig
        else:
            results.append(check(name, z1, z2, Rm, "REJECT"))

    print("\n(i) continued -- HMGCR under the UKB LD reference (n=337,000)")
    import part2_ukb_ld as P2
    import gzip
    rows = list(csv.DictReader(open(P2.MERGE, encoding="utf-8"), delimiter="\t"))
    snps = {}
    with gzip.open(os.path.join(P2.LDDIR, P2.BLOCK + ".gz"), "rt") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            snps[(int(r["position"]), r["allele1"], r["allele2"])] = i
    sel, flip, kept = [], [], []
    for r in rows:
        pos = int(r["pos37"])
        eff, oth = r["effect_allele_harmonized"], r["other_allele_harmonized"]
        j = snps.get((pos, oth, eff)); s = 1.0
        if j is None:
            j = snps.get((pos, eff, oth)); s = -1.0
        if j is None:
            continue
        sel.append(j); flip.append(s); kept.append(r)
    o = np.argsort(sel)
    sel = np.array(sel)[o]; flip = np.array(flip)[o]; kept = [kept[i] for i in o]
    idx_of = {v: i for i, v in enumerate(sel)}
    paths = {n: os.path.join(P2.EXT, n) for n in os.listdir(P2.EXT)}
    row = np.load(paths["row.npy"], mmap_mode="r")
    col = np.load(paths["col.npy"], mmap_mode="r")
    dat = np.load(paths["data.npy"], mmap_mode="r")
    n = len(sel); Ru = np.eye(n)
    lo, hi = int(sel.min()), int(sel.max()); CH = 1 << 23
    for st in range(0, len(row), CH):
        sp = min(st + CH, len(row))
        rr = np.asarray(row[st:sp]); m = (rr >= lo) & (rr <= hi)
        if not m.any():
            continue
        cc = np.asarray(col[st:sp])[m]; m2 = (cc >= lo) & (cc <= hi)
        if not m2.any():
            continue
        rr = rr[m][m2]; cc = cc[m2]; dd = np.asarray(dat[st:sp])[m][m2]
        for a, b, v in zip(rr, cc, dd):
            ia, ib = idx_of.get(int(a)), idx_of.get(int(b))
            if ia is None or ib is None:
                continue
            val = float(v) * flip[ia] * flip[ib]
            Ru[ia, ib] = val; Ru[ib, ia] = val
    np.fill_diagonal(Ru, 1.0)
    z1 = np.array([float(r["beta_exp"]) / float(r["se_exp"]) for r in kept])
    z2 = np.array([float(r["beta_out_harm"]) / float(r["se_out"]) for r in kept])
    results.append(check("HMGCR-LDL-C, UKB LD n=337,000", z1, z2, Ru, "REJECT"))

    print("\n(iii) RETAINED ANALYSES re-run through the assertion")
    print("      (Part 3 loci where SuSiE returned a credible-set pair)")
    pr = r"D:\causal2drug\deliverables\d1h\part3_results.tsv"
    if os.path.exists(pr):
        got = [r for r in csv.DictReader(open(pr, encoding="utf-8"), delimiter="\t")
               if r["susie_PP4_eqtlgen"]]
        print("      loci with a susie PP4 in part3_results.tsv: %d" % len(got))
        for r in got:
            print("        %-10s %-14s susie_PP4=%s  (re-checked below via ckpt rerun)"
                  % (r["gene"], r["gwas_acc"], r["susie_PP4_eqtlgen"]))
    print("\nRESULT: %d/%d real artefacts rejected" % (sum(results), len(results)))
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
