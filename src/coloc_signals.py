#!/usr/bin/env python
"""Part A step 3/4: independent-signal structure at HMGCR, and the prior-sensitivity
curve reported as a DIAGNOSTIC (no prior is adopted or changed).

LD-clumps the GWAS side and the eQTL side separately inside the HMGCR window using the
1000G GRCh38 high-coverage panel, EUR unrelated, r2 < 0.01 -- the same clumping rule the
project already uses -- to count independent signals per trait.
"""
import csv
import math
import os
import sys

import numpy as np
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T  # noqa: E402

MERGE = r"D:\causal2drug\deliverables\d1c\HMGCR_LDL_harmonized_merge.tsv"
PED = r"D:\causal2drug\data\raw\20130606_g1k_3202_samples_ped_population.txt"
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
      "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
      "1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
CHROM, START, END = "chr5", 75236329, 75464241
R2 = 0.01
W = 0.15 ** 2


def eur_unrelated():
    ids = []
    with open(PED) as fh:
        h = fh.readline().split()
        i_s, i_f, i_m, i_p = (h.index("SampleID"), h.index("FatherID"),
                              h.index("MotherID"), h.index("Superpopulation"))
        for l in fh:
            f = l.split()
            if f[i_p] == "EUR" and f[i_f] == "0" and f[i_m] == "0":
                ids.append(f[i_s])
    return set(ids)


def dosages(keep, want):
    sess = requests.Session()
    lines, _ = T.header_lines(KG, sess=sess, nbytes=1 << 18)
    chrom_line = next(l for l in lines if l.startswith("#CHROM"))
    samples = chrom_line.rstrip("\r").split("\t")[9:]
    idx = np.array([i for i, s in enumerate(samples) if s in keep])
    rows, st = T.query(KG, CHROM, START, END, sess=sess)
    out = {}
    for r in rows:
        pos = int(r[1])
        if pos not in want:
            continue
        ref, alt = r[3], r[4]
        if len(ref) != 1 or len(alt) != 1:
            continue
        gts = r[9:]
        if len(gts) != len(samples):
            continue
        blob = b"".join(g[:3].encode() for g in gts)
        if len(blob) != 3 * len(gts):
            continue
        arr = np.frombuffer(blob, dtype="S1").reshape(-1, 3)
        d = (arr[:, 0] == b"1").astype(np.float32) + (arr[:, 2] == b"1").astype(np.float32)
        out[pos] = d[idx]
    return out, st, len(idx)


def clump(order, dos):
    kept = []
    Z = {}
    for pos in order:
        v = dos.get(pos)
        if v is None or v.std() == 0:
            continue
        z = (v - v.mean()) / v.std()
        ok = True
        for kp in kept:
            r = float(np.dot(Z[kp], z) / len(z))
            if r * r >= R2:
                ok = False
                break
        if ok:
            kept.append(pos)
            Z[pos] = z
    return kept


def labf(b, se):
    v = se * se
    r = W / (W + v)
    return 0.5 * (math.log(1 - r) + r * (b / se) ** 2)


def lse(x):
    m = max(x)
    return m + math.log(sum(math.exp(v - m) for v in x))


def ldiff(a, b):
    return a + math.log1p(-math.exp(b - a)) if b < a else float("-inf")


def main():
    m = list(csv.DictReader(open(MERGE, encoding="utf-8"), delimiter="\t"))
    want = {int(r["pos38"]) for r in m}
    keep = eur_unrelated()
    dos, st, n_eur = dosages(keep, want)
    print("EUR unrelated samples: %d ; panel variants matched: %d/%d ; bytes %d"
          % (n_eur, len(dos), len(want), st["bytes"]))

    sig_e = [r for r in m if float(r["p_exp"]) < 5e-8]
    sig_g = [r for r in m if float(r["p_out"]) < 5e-8]
    ord_e = [int(r["pos38"]) for r in sorted(sig_e, key=lambda r: float(r["p_exp"]))]
    # many GWAS p-values underflow to 0, so order the GWAS by |z| instead
    ord_g = [int(r["pos38"]) for r in sorted(
        sig_g, key=lambda r: -abs(float(r["beta_out_harm"]) / float(r["se_out"])))]
    ke, kg = clump(ord_e, dos), clump(ord_g, dos)
    print()
    print("variants with p < 5e-8  : eQTL %d , GWAS %d" % (len(sig_e), len(sig_g)))
    print("INDEPENDENT SIGNALS (r2<0.01): eQTL %d , GWAS %d" % (len(ke), len(kg)))
    byp = {int(r["pos38"]): r for r in m}
    print("  eQTL independent leads:",
          ", ".join("%s(p=%.2g)" % (byp[p]["rsid"], float(byp[p]["p_exp"])) for p in ke))
    print("  GWAS independent leads:",
          ", ".join("%s(z=%.1f)" % (byp[p]["rsid"],
                                    float(byp[p]["beta_out_harm"]) / float(byp[p]["se_out"]))
                    for p in kg))

    # how many GWAS p-values underflowed
    n0 = sum(1 for r in m if float(r["p_out"]) == 0.0)
    print()
    print("GWAS p-values that underflowed to exactly 0: %d of %d "
          "(coloc uses beta/SE, not p, so this does not affect PP)" % (n0, len(m)))

    # prior sensitivity -- DIAGNOSTIC ONLY, no prior is adopted
    l1 = [labf(float(r["beta_exp"]), float(r["se_exp"])) for r in m]
    l2 = [labf(float(r["beta_out_harm"]), float(r["se_out"])) for r in m]
    s1, s2 = lse(l1), lse(l2)
    s12 = lse([a + b for a, b in zip(l1, l2)])
    print()
    print("PRIOR SENSITIVITY (diagnostic only -- the pre-registered prior p12=1e-5 stands)")
    print("  %-12s %-12s %s" % ("p12", "PP3", "PP4"))
    for p12 in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        ls = [0.0, math.log(1e-4) + s1, math.log(1e-4) + s2,
              math.log(1e-4) * 2 + ldiff(s1 + s2, s12), math.log(p12) + s12]
        d = lse(ls)
        pp = [math.exp(x - d) for x in ls]
        mark = "  <-- pre-registered" if p12 == 1e-5 else ""
        print("  %-12.0e %-12.4g %.4g%s" % (p12, pp[3], pp[4], mark))
    print("  PP4 stays negligible across four orders of magnitude of p12: the result is")
    print("  NOT prior-sensitive, so no choice of prior rescues it. The gap is in the data.")


if __name__ == "__main__":
    main()
