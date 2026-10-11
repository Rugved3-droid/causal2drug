#!/usr/bin/env python
"""PART 2: single-shot test of a LARGER LD reference at the HMGCR locus.

Reference: UK Biobank LD matrices released with PolyFun (Weissbrod et al.), bucket
broad-alkesgroup-ukbb-ld, block chr5_74000001_77000001 -- n = 337,000 unrelated British
ancestry individuals, GRCh37, 3 Mb blocks.

Pre-declared purpose: answer the reviewer objection that 1000G EUR n=525 is simply too
small, with one number. One reference, one result, whatever it returns (Amendment 3).

The .npz is a 218-million-entry sparse COO matrix (~2.6 GB when materialised) and this
machine has ~1.9 GB free, so the members are extracted to disk and memory-mapped, then
filtered in chunks.

The residual-z RMSE reported here is the DENTIST statistic (Chen et al. 2021,
Nat Commun 12:7117), reimplemented and cited as theirs.
"""
import csv
import gzip
import json
import math
import os
import sys
import zipfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import susie_coloc as SC  # noqa: E402

LDDIR = r"D:\causal2drug\data\raw\ukbld"
BLOCK = "chr5_74000001_77000001"
EXT = os.path.join(LDDIR, "extracted")
MERGE = r"D:\causal2drug\deliverables\d1c\HMGCR_LDL_harmonized_merge.tsv"
OUT = r"D:\causal2drug\deliverables\d1i"
W = 0.15 ** 2


def extract():
    os.makedirs(EXT, exist_ok=True)
    z = zipfile.ZipFile(os.path.join(LDDIR, BLOCK + ".npz"))
    for name in z.namelist():
        dst = os.path.join(EXT, name)
        if not os.path.exists(dst):
            print("  extracting %s ..." % name, flush=True)
            with z.open(name) as src, open(dst, "wb") as fh:
                while True:
                    b = src.read(1 << 24)
                    if not b:
                        break
                    fh.write(b)
    return {n: os.path.join(EXT, n) for n in z.namelist()}


def labf(z, se):
    v = se * se
    r = W / (W + v)
    return 0.5 * (math.log(1 - r) + r * z * z)


def dentist_rmse(z, R):
    i = int(np.argmax(np.abs(z)))
    pred = z[i] * R[i]
    m = np.abs(R[i]) > 0.2
    return float(np.sqrt(np.mean((z[m] - pred[m]) ** 2))), int(m.sum()), i


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(MERGE, encoding="utf-8"), delimiter="\t"))
    print("HMGCR merged variants (eQTLGen x LDL GWAS): %d" % len(rows))

    # --- UKB snp list, GRCh37 ---
    snps = {}
    with gzip.open(os.path.join(LDDIR, BLOCK + ".gz"), "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for i, r in enumerate(rd):
            snps[(int(r["position"]), r["allele1"], r["allele2"])] = i
    print("UKB LD block %s: %d variants (GRCh37)" % (BLOCK, len(snps)))

    # match our variants by GRCh37 position + allele pair (either orientation)
    sel, flip, kept = [], [], []
    for r in rows:
        p = int(r["pos37"])
        eff, oth = r["effect_allele_harmonized"], r["other_allele_harmonized"]
        j = snps.get((p, oth, eff))
        s = 1.0
        if j is None:
            j = snps.get((p, eff, oth))
            s = -1.0
        if j is None:
            continue
        sel.append(j)
        flip.append(s)
        kept.append(r)
    print("matched into the UKB LD block: %d / %d" % (len(kept), len(rows)))
    if len(kept) < 30:
        print("RESULT: too few variants matched; test cannot proceed.")
        return

    order = np.argsort(sel)
    sel = np.array(sel)[order]
    flip = np.array(flip)[order]
    kept = [kept[i] for i in order]
    idx_of = {v: i for i, v in enumerate(sel)}
    selset = set(int(v) for v in sel)

    paths = extract()
    row = np.load(paths["row.npy"], mmap_mode="r")
    col = np.load(paths["col.npy"], mmap_mode="r")
    dat = np.load(paths["data.npy"], mmap_mode="r")
    shp = np.load(paths["shape.npy"])
    print("sparse LD: %d stored entries, matrix shape %s" % (len(row), tuple(shp)))

    n = len(sel)
    R = np.eye(n)
    lo, hi = int(sel.min()), int(sel.max())
    CH = 1 << 23
    nfill = 0
    for start in range(0, len(row), CH):
        stop = min(start + CH, len(row))
        rr = np.asarray(row[start:stop])
        m = (rr >= lo) & (rr <= hi)
        if not m.any():
            continue
        cc = np.asarray(col[start:stop])[m]
        m2 = (cc >= lo) & (cc <= hi)
        if not m2.any():
            continue
        rr = rr[m][m2]
        cc = cc[m2]
        dd = np.asarray(dat[start:stop])[m][m2]
        for a, b, v in zip(rr, cc, dd):
            ia, ib = idx_of.get(int(a)), idx_of.get(int(b))
            if ia is None or ib is None:
                continue
            val = float(v) * flip[ia] * flip[ib]
            R[ia, ib] = val
            R[ib, ia] = val
            nfill += 1
    np.fill_diagonal(R, 1.0)
    print("filled %d off-diagonal LD entries for the %dx%d submatrix" % (nfill, n, n))

    z1 = np.array([float(r["beta_exp"]) / float(r["se_exp"]) for r in kept])
    se1 = np.array([float(r["se_exp"]) for r in kept])
    z2 = np.array([float(r["beta_out_harm"]) / float(r["se_out"]) for r in kept])
    se2 = np.array([float(r["se_out"]) for r in kept])

    r_e, n_e, i1 = dentist_rmse(z1, R)
    r_g, n_g, i2 = dentist_rmse(z2, R)
    a = SC.coloc_from_lbf(np.array([labf(x, s) for x, s in zip(z1, se1)]),
                          np.array([labf(x, s) for x, s in zip(z2, se2)]))
    s = SC.coloc_susie(z1, z2, R, L=5)
    cond = SC.ld_condition_number(R)

    print()
    print("=" * 74)
    print("PART 2 RESULT - UKB LD reference (n=337,000 British), GRCh37")
    print("=" * 74)
    print("  variants used                     : %d" % n)
    print("  LD condition number               : %.4g (min eig %.3g)"
          % (cond["condition_number"], cond["min_eigenvalue"]))
    print("  eQTL lead                         : %s z=%.2f"
          % (kept[i1]["rsid"], z1[i1]))
    print("  GWAS lead                         : %s z=%.2f"
          % (kept[i2]["rsid"], z2[i2]))
    print("  DENTIST residual-z RMSE, GWAS     : %.3f   (1000G n=525 gave 7.80)" % r_g)
    print("  DENTIST residual-z RMSE, eQTL     : %.3f   (1000G n=525 gave 1.33)" % r_e)
    print("  coloc.abf   PP3 / PP4             : %.4g / %.4g" % (a["PP3"], a["PP4"]))
    print("  coloc.susie best PP4              : %s"
          % (("%.4g" % s["best"]["PP4"]) if s["best"] else "no credible set"))
    print("  credible sets exposure / outcome  : %d / %d" % (s["n_cs1"], s["n_cs2"]))
    tops = [kept[c["top"]]["rsid"] for c in s["cs2"]]
    print("  outcome CS top variants           : %s" % ", ".join(tops))
    print("  any CS tops at rs12916 (z=49.59)  : %s"
          % ("YES" if "rs12916" in tops else "NO"))
    json.dump(dict(reference="UKB LD (PolyFun/Weissbrod), broad-alkesgroup-ukbb-ld",
                   block=BLOCK, n_reference=337000, ancestry="British (UK Biobank)",
                   build="GRCh37", n_variants=n,
                   ld_condition_number=cond["condition_number"],
                   dentist_rmse_gwas=r_g, dentist_rmse_eqtl=r_e,
                   abf_PP3=a["PP3"], abf_PP4=a["PP4"],
                   susie_PP4=(s["best"]["PP4"] if s["best"] else None),
                   n_cs_exposure=s["n_cs1"], n_cs_outcome=s["n_cs2"],
                   outcome_cs_tops=tops, cs_tops_rs12916=("rs12916" in tops)),
              open(os.path.join(OUT, "part2_ukb_ld.json"), "w"), indent=2)
    print("\nwrote", os.path.join(OUT, "part2_ukb_ld.json"))


if __name__ == "__main__":
    main()
