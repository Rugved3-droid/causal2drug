#!/usr/bin/env python
"""Run coloc.abf and coloc.susie side by side on the real control loci.

LD: 1000 Genomes GRCh38 high-coverage panel, EUR unrelated (n=525) -- the same panel used
for clumping throughout the project. Dosages are re-oriented to the harmonized effect
allele (the eQTLGen assessed allele) before the correlation matrix is formed, so LD signs
match the z-score signs. Getting that wrong silently breaks SuSiE.
"""
import csv
import json
import math
import os
import sys

import numpy as np
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T          # noqa: E402
import susie_coloc as SC        # noqa: E402

PED = r"D:\causal2drug\data\raw\20130606_g1k_3202_samples_ped_population.txt"
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
      "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
      "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
OUT = r"D:\causal2drug\deliverables\d1g"
W = 0.15 ** 2

LOCI = [
    dict(name="HMGCR_LDL", merge=r"D:\causal2drug\deliverables\d1c\HMGCR_LDL_harmonized_merge.tsv",
         chrom="5"),
    dict(name="IL6R_CRP", merge=r"D:\causal2drug\deliverables\d1c\IL6R_CRP_harmonized_merge.tsv",
         chrom="1"),
]


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


def ld_matrix(chrom, rows, keep):
    """Return (R, kept_row_indices, n_eur). Dosages oriented to the harmonized
    effect allele so that LD sign matches z sign."""
    sess = requests.Session()
    url = KG % chrom
    lines, _ = T.header_lines(url, sess=sess, nbytes=1 << 18)
    chrom_line = next(l for l in lines if l.startswith("#CHROM"))
    samples = chrom_line.rstrip("\r").split("\t")[9:]
    idx = np.array([i for i, s in enumerate(samples) if s in keep])
    pos = [int(r["pos38"]) for r in rows]
    recs, st = T.query(url, "chr" + chrom, min(pos) - 1000, max(pos) + 1000, sess=sess)
    panel = {}
    for r in recs:
        p = int(r[1])
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
        d = (arr[:, 0] == b"1").astype(np.float64) + (arr[:, 2] == b"1").astype(np.float64)
        panel[(p, ref, alt)] = d[idx]

    vecs, keep_i = [], []
    for i, r in enumerate(rows):
        p = int(r["pos38"])
        eff = r["effect_allele_harmonized"]
        oth = r["other_allele_harmonized"]
        v = panel.get((p, oth, eff))          # effect allele is VCF ALT
        if v is None:
            v = panel.get((p, eff, oth))      # effect allele is VCF REF -> flip dosage
            if v is None:
                continue
            v = 2.0 - v
        if v.std() == 0:
            continue
        vecs.append((v - v.mean()) / v.std())
        keep_i.append(i)
    if not vecs:
        return None, [], len(idx), st
    X = np.vstack(vecs)
    R = (X @ X.T) / X.shape[1]
    np.fill_diagonal(R, 1.0)
    return R, keep_i, len(idx), st


def abf(z1, se1, z2, se2):
    def labf(z, se):
        v = se * se
        r = W / (W + v)
        return 0.5 * (math.log(1 - r) + r * z * z)
    return SC.coloc_from_lbf(np.array([labf(a, b) for a, b in zip(z1, se1)]),
                             np.array([labf(a, b) for a, b in zip(z2, se2)]))


def main():
    os.makedirs(OUT, exist_ok=True)
    keep = eur_unrelated()
    summary = []
    for loc in LOCI:
        rows = list(csv.DictReader(open(loc["merge"], encoding="utf-8"), delimiter="\t"))
        print("\n" + "=" * 78)
        print("%s   merged variants: %d" % (loc["name"], len(rows)))
        R, ki, n_eur, st = ld_matrix(loc["chrom"], rows, keep)
        rows = [rows[i] for i in ki]
        print("  LD: 1000G GRCh38 EUR unrelated n=%d ; variants with LD: %d ; bytes %d"
              % (n_eur, len(rows), st["bytes"]))
        cond = SC.ld_condition_number(R)
        print("  LD conditioning: min eig %.3g, max eig %.3g, condition number %.3g"
              % (cond["min_eigenvalue"], cond["max_eigenvalue"], cond["condition_number"]))
        z1 = np.array([float(r["beta_exp"]) / float(r["se_exp"]) for r in rows])
        z2 = np.array([float(r["beta_out_harm"]) / float(r["se_out"]) for r in rows])
        se1 = np.array([float(r["se_exp"]) for r in rows])
        se2 = np.array([float(r["se_out"]) for r in rows])

        a = abf(z1, se1, z2, se2)
        s = SC.coloc_susie(z1, z2, R, L=5)
        i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
        r2 = float(R[i1, i2] ** 2)
        d = abs(int(rows[i1]["pos38"]) - int(rows[i2]["pos38"]))
        print("  eQTL lead %s z=%.2f ; GWAS lead %s z=%.2f ; r2=%.3f ; distance %d bp"
              % (rows[i1]["rsid"], z1[i1], rows[i2]["rsid"], z2[i2], r2, d))
        print("  single-variant coloc.abf : PP3=%.4g PP4=%.4g  (gap s1+s2-s12 = %.2f)"
              % (a["PP3"], a["PP4"], a["gap"]))
        print("  credible sets: exposure %d, outcome %d" % (s["n_cs1"], s["n_cs2"]))
        for pr in sorted(s["pairs"], key=lambda x: -x["PP4"]):
            print("    CS%d(size %d, purity %.2f, top %s) x CS%d(size %d, purity %.2f, top %s)"
                  "  PP3=%.4g PP4=%.4g"
                  % (pr["cs1_component"], pr["cs1_size"], pr["cs1_purity"],
                     rows[pr["cs1_top"]]["rsid"],
                     pr["cs2_component"], pr["cs2_size"], pr["cs2_purity"],
                     rows[pr["cs2_top"]]["rsid"], pr["PP3"], pr["PP4"]))
        best = s["best"]
        print("  coloc.susie BEST PP4 = %s"
              % (("%.4g" % best["PP4"]) if best else "no credible set"))
        summary.append(dict(locus=loc["name"], n_variants=len(rows), n_eur=n_eur,
                            ld_condition_number=cond["condition_number"],
                            ld_min_eigenvalue=cond["min_eigenvalue"],
                            eqtl_lead=rows[i1]["rsid"], eqtl_lead_z=float(z1[i1]),
                            gwas_lead=rows[i2]["rsid"], gwas_lead_z=float(z2[i2]),
                            lead_r2=r2, lead_distance_bp=d,
                            abf_PP3=a["PP3"], abf_PP4=a["PP4"], abf_gap=a["gap"],
                            n_cs_exposure=s["n_cs1"], n_cs_outcome=s["n_cs2"],
                            susie_best_PP4=(best["PP4"] if best else None),
                            cs_pairs=[{k: v for k, v in pr.items()
                                       if k in ("cs1_component", "cs2_component", "PP3",
                                                "PP4", "cs1_size", "cs2_size",
                                                "cs1_purity", "cs2_purity")}
                                      for pr in s["pairs"]]))
    json.dump(summary, open(os.path.join(OUT, "susie_controls.json"), "w"), indent=2)
    print("\nwrote", os.path.join(OUT, "susie_controls.json"))


if __name__ == "__main__":
    main()
