#!/usr/bin/env python
"""PART A: diagnose the colocalization defect at HMGCR - LDL-C.

Three colocalizations on the SAME locus and the SAME outcome GWAS:

  (1) eQTLGen whole blood, beta/SE RECONSTRUCTED from Z via the SMR approximation
  (2) eQTL Catalogue GTEx v10 whole blood QTD000356 (n=853), NATIVE beta/SE -- no
      reconstruction. INTERVAL (n=4,729) was the intended comparator but QTS000053 is
      not released on the FTP (HTTP 404); this is the largest AVAILABLE whole-blood
      'ge' dataset. See deliverables/d1d/PART_A_COLOC_DIAGNOSIS.md.
  (3) the SAME QTD000356 data, but with beta/SE re-derived from its own z by the SMR
      formula using its own allele frequency and sample size

(1) vs (2) is the test the user asked for. (3) vs (2) is stronger: it isolates the
Z-score conversion exactly, holding the dataset, the locus and the outcome fixed, so any
difference is attributable to the conversion alone and nothing else.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T  # noqa: E402

OUT = r"D:\causal2drug\deliverables\d1d"
GWAS = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz")
INTERVAL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000015/QTD000356/QTD000356.all.tsv.gz")
ENSG = "ENSG00000113161"
CHROM, START, END = "5", 75236329, 75464241   # HMGCR gene body +-100kb, GRCh38
W_QUANT = 0.15 ** 2
P1, P2, P12 = 1e-4, 1e-4, 1e-5

EQ_COLS = ("molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue "
           "beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid").split()


def labf(beta, se, W=W_QUANT):
    v = se * se
    r = W / (W + v)
    z = beta / se
    return 0.5 * (math.log(1 - r) + r * z * z)


def lse(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def ldiff(a, b):
    return a + math.log1p(-math.exp(b - a)) if b < a else float("-inf")


def coloc(pairs):
    """pairs: list of (beta1, se1, beta2, se2)."""
    l1 = [labf(b1, s1) for b1, s1, _b2, _s2 in pairs]
    l2 = [labf(b2, s2) for _b1, _s1, b2, s2 in pairs]
    s1, s2 = lse(l1), lse(l2)
    s12 = lse([a + b for a, b in zip(l1, l2)])
    ls = [0.0, math.log(P1) + s1, math.log(P2) + s2,
          math.log(P1) + math.log(P2) + ldiff(s1 + s2, s12), math.log(P12) + s12]
    d = lse(ls)
    pp = [math.exp(x - d) for x in ls]
    return dict(n=len(pairs), PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4],
                s1=s1, s2=s2, s12=s12, gap=s1 + s2 - s12)


def smr_beta_se(z, af, n):
    v = 2.0 * af * (1.0 - af) * (n + z * z)
    se = 1.0 / math.sqrt(v)
    return z * se, se


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---- outcome GWAS (GRCh38) ----
    lines, _ = T.header_lines(GWAS, nbytes=1 << 16)
    gcols = lines[0].rstrip("\r").split("\t")
    gi = {c: i for i, c in enumerate(gcols)}
    grows, gst = T.query(GWAS, CHROM, START, END)
    gmap = {}
    for g in grows:
        try:
            gmap[int(g[gi["base_pair_location"]])] = (
                g[gi["effect_allele"]], g[gi["other_allele"]],
                float(g[gi["beta"]]), float(g[gi["standard_error"]]),
                float(g[gi["p_value"]]))
        except ValueError:
            pass
    print("GWAS variants in window: %d (bytes %d)" % (len(gmap), gst["bytes"]))

    # ---- INTERVAL whole blood, native beta/se ----
    qi = {c: i for i, c in enumerate(EQ_COLS)}
    qrows, qst = T.query(INTERVAL, CHROM, START, END)
    qrows = [r for r in qrows if r[qi["gene_id"]] == ENSG]
    print("GTEx blood HMGCR cis variants in window: %d (bytes %d)" % (len(qrows), qst["bytes"]))

    native, recon = [], []
    detail = []
    for r in qrows:
        pos = int(r[qi["position"]])
        g = gmap.get(pos)
        if not g:
            continue
        eff, oth, b_out, se_out, p_out = g
        ref, alt = r[qi["ref"]], r[qi["alt"]]
        if (eff, oth) == (alt, ref):
            sign = 1.0
        elif (eff, oth) == (ref, alt):
            sign = -1.0
        else:
            continue
        try:
            b_exp = float(r[qi["beta"]])
            se_exp = float(r[qi["se"]])
            ac = float(r[qi["ac"]])
            an = float(r[qi["an"]])
        except ValueError:
            continue
        if se_exp <= 0 or se_out <= 0 or an <= 0:
            continue
        af = ac / an
        if not (0 < af < 1):
            continue
        n_eff = an / 2.0
        z = b_exp / se_exp
        b_rec, se_rec = smr_beta_se(z, af, n_eff)
        native.append((b_exp, se_exp, sign * b_out, se_out))
        recon.append((b_rec, se_rec, sign * b_out, se_out))
        detail.append(dict(pos=pos, rsid=r[qi["rsid"]], ref=ref, alt=alt,
                           beta_native=b_exp, se_native=se_exp,
                           beta_smr=b_rec, se_smr=se_rec, z=z, af=af, n=n_eff,
                           p_exp=float(r[qi["pvalue"]]),
                           beta_out=sign * b_out, se_out=se_out, p_out=p_out))
    print("overlapping variants (GTEx blood x GWAS): %d" % len(native))

    with open(os.path.join(OUT, "hmgcr_gtexblood_merge.tsv"), "w", newline="",
              encoding="utf-8") as fh:
        cols = list(detail[0].keys())
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for d in detail:
            w.writerow([d[c] for c in cols])

    r_nat = coloc(native)
    r_rec = coloc(recon)

    # ---- (1) eQTLGen converted, from the Part C merge table ----
    m = list(csv.DictReader(
        open(r"D:\causal2drug\deliverables\d1c\HMGCR_LDL_harmonized_merge.tsv",
             encoding="utf-8"), delimiter="\t"))
    eqg = [(float(x["beta_exp"]), float(x["se_exp"]),
            float(x["beta_out_harm"]), float(x["se_out"])) for x in m]
    r_eqg = coloc(eqg)

    print()
    hdr = "%-46s %6s %9s %9s %9s %9s %9s %8s"
    print(hdr % ("colocalization", "n", "PP0", "PP1", "PP2", "PP3", "PP4", "s1+s2-s12"))
    print("-" * 112)
    for lab, r in (("(1) eQTLGen, beta/SE RECONSTRUCTED from Z", r_eqg),
                   ("(2) GTEx v10 whole blood, NATIVE beta/SE", r_nat),
                   ("(3) GTEx blood, beta/SE re-derived by SMR", r_rec)):
        print(hdr % (lab, r["n"], "%.3g" % r["PP0"], "%.3g" % r["PP1"],
                     "%.3g" % r["PP2"], "%.4f" % r["PP3"], "%.3g" % r["PP4"],
                     "%.2f" % r["gap"]))
    print()
    print("(2) vs (3) isolates the Z-score conversion: same dataset, same locus, same GWAS.")

    # scale comparison of native vs reconstructed
    import statistics as st
    rb = [d["beta_smr"] / d["beta_native"] for d in detail if d["beta_native"] != 0]
    rs = [d["se_smr"] / d["se_native"] for d in detail]
    print()
    print("beta_smr / beta_native : median %.4f  (IQR %.4f - %.4f)"
          % (st.median(rb), sorted(rb)[len(rb) // 4], sorted(rb)[3 * len(rb) // 4]))
    print("se_smr   / se_native   : median %.4f  (IQR %.4f - %.4f)"
          % (st.median(rs), sorted(rs)[len(rs) // 4], sorted(rs)[3 * len(rs) // 4]))
    print("(ratios ~1.0 mean the SMR reconstruction reproduces the native values)")

    import json
    json.dump(dict(eqtlgen_reconstructed=r_eqg, gtexblood_native=r_nat,
                   gtexblood_smr_rederived=r_rec,
                   beta_ratio_median=st.median(rb), se_ratio_median=st.median(rs)),
              open(os.path.join(OUT, "coloc_diagnosis.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
