#!/usr/bin/env python
"""PART 3: degradation measurement (PREREGISTRATION R14).

Per locus, using the SAME outcome GWAS, window and harmonization throughout:
  - coloc.abf with eQTLGen (n=31,684) as exposure
  - coloc.abf with GTEx whole blood (QTD000356, n=853) as exposure   <-- head-to-head
  - coloc.susie with eQTLGen as exposure
  - residual-z RMSE for every trait -- the DENTIST statistic (Chen et al. 2021,
    Nat Commun 12:7117), reimplemented here and cited as theirs
  - eQTL/GWAS lead z, r2 between leads, distance, lBF gap s1+s2-s12

Checkpoints per locus so an interruption costs at most one locus.
"""
import collections
import csv
import gzip
import math
import os
import sys
import time

import numpy as np
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T                 # noqa: E402
import susie_coloc as SC               # noqa: E402
from harmonize import harmonize_one    # noqa: E402
from pyliftover import LiftOver        # noqa: E402

D1B = r"D:\causal2drug\deliverables\d1b"
OUT = r"D:\causal2drug\deliverables\d1h"
CKPT = OUT + r"\ckpt"
RAW = r"D:\causal2drug\data\raw"
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
GTEX_BLOOD = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
              "QTS000015/QTD000356/QTD000356.all.tsv.gz")
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
      "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
      "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
PED = RAW + r"\20130606_g1k_3202_samples_ped_population.txt"
FLANK = 100_000
W = 0.15 ** 2
EQ_COLS = ("molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue "
           "beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid").split()

COLS = ["gene", "gwas_acc", "mesh_term", "n_merged_eqtlgen", "n_merged_gtex",
        "eqtl_lead_rsid", "eqtl_lead_z", "gwas_lead_rsid", "gwas_lead_z",
        "lead_r2", "lead_dist_bp", "abf_gap",
        "abf_PP3_eqtlgen", "abf_PP4_eqtlgen", "susie_PP4_eqtlgen",
        "n_cs_exposure", "n_cs_outcome",
        "abf_PP4_gtexblood", "gtex_lead_z", "n_lead_shared_gtex_eqtlgen",
        "dentist_rmse_eqtlgen", "dentist_rmse_gwas", "dentist_rmse_gtex", "note"]


def eur():
    ids = []
    with open(PED) as fh:
        h = fh.readline().split()
        a, b, c, d = (h.index("SampleID"), h.index("FatherID"), h.index("MotherID"),
                      h.index("Superpopulation"))
        for l in fh:
            f = l.split()
            if f[d] == "EUR" and f[b] == "0" and f[c] == "0":
                ids.append(f[a])
    return set(ids)


def labf(z, se):
    v = se * se
    r = W / (W + v)
    return 0.5 * (math.log(1 - r) + r * z * z)


def abf(z1, se1, z2, se2):
    return SC.coloc_from_lbf(np.array([labf(a, b) for a, b in zip(z1, se1)]),
                             np.array([labf(a, b) for a, b in zip(z2, se2)]))


def dentist_rmse(z, R):
    """Chen et al. 2021 (DENTIST): observed vs LD-predicted test statistic.
    Predict every variant from the lead using the reference LD; RMSE ~1 if consistent."""
    i = int(np.argmax(np.abs(z)))
    pred = z[i] * R[i]
    m = np.abs(R[i]) > 0.2
    if m.sum() < 5:
        return float("nan")
    return float(np.sqrt(np.mean((z[m] - pred[m]) ** 2)))


def load_af(rsids):
    af = {}
    with gzip.open(RAW + r"\eqtlgen_SNP_AF.txt.gz", "rt", encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_s, i_f, i_a, i_b = (h.index("SNP"), h.index("AlleleB_all"),
                              h.index("AlleleA"), h.index("AlleleB"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_s] in rsids:
                try:
                    af[f[i_s]] = (float(f[i_f]), f[i_a], f[i_b])
                except ValueError:
                    pass
    return af


def smr(z, a, n):
    v = 2.0 * a * (1.0 - a) * (n + z * z)
    if v <= 0:
        return None, None
    se = 1.0 / math.sqrt(v)
    return z * se, se


def main():
    os.makedirs(CKPT, exist_ok=True)
    loci = list(csv.DictReader(open(OUT + r"\loci.tsv", encoding="utf-8"), delimiter="\t"))
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D1B + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    import re
    harm = {}
    for l in open(RAW + r"\harmonised_list.txt"):
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m:
            harm[m.group(1)] = l.strip().lstrip("./")

    genes = {r["gene"] for r in loci}
    hits = collections.defaultdict(list)
    with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = r["MinikelSymbol"]
            if s not in genes:
                continue
            c = coords[s]
            if r["SNPChr"] != c["chrom"]:
                continue
            if int(c["start"]) - FLANK <= int(r["SNPPos"]) <= int(c["end"]) + FLANK:
                hits[s].append(r)
    need = {r["SNP"] for v in hits.values() for r in v}
    print("loci %d ; exposure variants %d ; loading AF..." % (len(loci), len(need)), flush=True)
    af = load_af(need)
    lo = LiftOver(RAW + r"\hg19ToHg38.over.chain.gz")
    keep = eur()

    done = set()
    for f in os.listdir(CKPT):
        for r in csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"), delimiter="\t"):
            done.add((r["gene"], r["gwas_acc"]))
    print("resuming: %d loci already done" % len(done), flush=True)

    t0 = time.time()
    for loc in loci:
        g, acc = loc["gene"], loc["gwas_acc"]
        if (g, acc) in done:
            continue
        rec = {c: "" for c in COLS}
        rec.update(gene=g, gwas_acc=acc, mesh_term=loc["mesh_term"])
        path = os.path.join(CKPT, "%s__%s.tsv" % (g, acc))
        try:
            ev = []
            for r in hits[g]:
                a = af.get(r["SNP"])
                if not a:
                    continue
                afb, alA, alB = a
                fa = afb if r["AssessedAllele"] == alB else (
                    1.0 - afb if r["AssessedAllele"] == alA else None)
                if fa is None:
                    continue
                z = float(r["Zscore"])
                b, se = smr(z, fa, float(r["NrSamples"]))
                if b is None:
                    continue
                res = lo.convert_coordinate("chr" + r["SNPChr"], int(r["SNPPos"]) - 1)
                if not res:
                    continue
                ev.append(dict(rsid=r["SNP"], chrom=res[0][0].replace("chr", ""),
                               pos38=res[0][1] + 1, assessed=r["AssessedAllele"],
                               other=r["OtherAllele"], af=fa, beta=b, se=se,
                               p=float(r["Pvalue"])))
            if not ev:
                rec["note"] = "no_exposure"
                raise StopIteration
            chrom = ev[0]["chrom"]
            lo_p, hi_p = min(e["pos38"] for e in ev) - 1000, max(e["pos38"] for e in ev) + 1000

            sess = requests.Session()
            url = GW + harm[acc]
            lines, _ = T.header_lines(url, sess=sess, nbytes=1 << 16)
            gc = lines[0].rstrip("\r").split("\t")
            gi = {c: i for i, c in enumerate(gc)}
            use_or = ("beta" not in gi) and ("odds_ratio" in gi)
            grows, _ = T.query(url, chrom, lo_p, hi_p, sess=sess)
            gmap = collections.defaultdict(list)
            for x in grows:
                try:
                    gmap[int(x[gi["base_pair_location"]])].append(x)
                except (ValueError, IndexError):
                    pass

            merged = []
            for e in ev:
                for x in gmap.get(e["pos38"], []):
                    try:
                        so = float(x[gi["standard_error"]])
                        bo = (math.log(float(x[gi["odds_ratio"]])) if use_or
                              else float(x[gi["beta"]]))
                    except (ValueError, KeyError):
                        continue
                    if so <= 0:
                        continue
                    try:
                        afo = float(x[gi["effect_allele_frequency"]])
                    except (ValueError, KeyError):
                        afo = None
                    h = harmonize_one(e["other"], e["assessed"], x[gi["effect_allele"]],
                                      x[gi["other_allele"]], bo,
                                      af_alt_exp=e["af"], af_effect_out=afo)
                    if h["sign"] is None:
                        continue
                    merged.append(dict(e, bo=h["beta_out_harmonized"], so=so))
                    break
            if len(merged) < 30:
                rec["note"] = "too_few_merged:%d" % len(merged)
                raise StopIteration

            # ---- LD (1000G EUR), oriented to the exposure assessed allele ----
            lines2, _ = T.header_lines(KG % chrom, sess=sess, nbytes=1 << 18)
            cl = next(l for l in lines2 if l.startswith("#CHROM"))
            samples = cl.rstrip("\r").split("\t")[9:]
            idx = np.array([i for i, s in enumerate(samples) if s in keep])
            recs, _ = T.query(KG % chrom, "chr" + chrom,
                              min(m["pos38"] for m in merged) - 500,
                              max(m["pos38"] for m in merged) + 500, sess=sess)
            panel = {}
            for r in recs:
                ref, alt = r[3], r[4]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                gts = r[9:]
                if len(gts) != len(samples):
                    continue
                blob = b"".join(q[:3].encode() for q in gts)
                if len(blob) != 3 * len(gts):
                    continue
                arr = np.frombuffer(blob, dtype="S1").reshape(-1, 3)
                d = ((arr[:, 0] == b"1").astype(np.float64)
                     + (arr[:, 2] == b"1").astype(np.float64))
                panel[(int(r[1]), ref, alt)] = d[idx]
            vecs, kept = [], []
            for m in merged:
                v = panel.get((m["pos38"], m["other"], m["assessed"]))
                if v is None:
                    v = panel.get((m["pos38"], m["assessed"], m["other"]))
                    if v is None:
                        continue
                    v = 2.0 - v
                if v.std() == 0:
                    continue
                vecs.append((v - v.mean()) / v.std())
                kept.append(m)
            if len(kept) < 30:
                rec["note"] = "too_few_with_LD:%d" % len(kept)
                raise StopIteration
            X = np.vstack(vecs)
            R = (X @ X.T) / X.shape[1]
            np.fill_diagonal(R, 1.0)

            z1 = np.array([m["beta"] / m["se"] for m in kept])
            se1 = np.array([m["se"] for m in kept])
            z2 = np.array([m["bo"] / m["so"] for m in kept])
            se2 = np.array([m["so"] for m in kept])
            a = abf(z1, se1, z2, se2)
            try:
                s = SC.coloc_susie(z1, z2, R, L=5, strict=True)
            except SC.ColocArtefact as ce:
                s = dict(best=None, n_cs1=-1, n_cs2=-1)
                rec["note"] = (rec["note"] + ";" if rec["note"] else "") +                     "coloc_artefact_rejected"
                print("      REJECTED by assertion: %s" % ce, flush=True)
            i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
            rec.update(n_merged_eqtlgen=len(kept),
                       eqtl_lead_rsid=kept[i1]["rsid"], eqtl_lead_z="%.4g" % z1[i1],
                       gwas_lead_rsid=kept[i2]["rsid"], gwas_lead_z="%.4g" % z2[i2],
                       lead_r2="%.4f" % (R[i1, i2] ** 2),
                       lead_dist_bp=abs(kept[i1]["pos38"] - kept[i2]["pos38"]),
                       abf_gap="%.3f" % a["gap"],
                       abf_PP3_eqtlgen="%.6g" % a["PP3"], abf_PP4_eqtlgen="%.6g" % a["PP4"],
                       susie_PP4_eqtlgen=("%.6g" % s["best"]["PP4"]) if s["best"] else "",
                       n_cs_exposure=s["n_cs1"], n_cs_outcome=s["n_cs2"],
                       dentist_rmse_eqtlgen="%.4g" % dentist_rmse(z1, R),
                       dentist_rmse_gwas="%.4g" % dentist_rmse(z2, R))

            # ---- head-to-head: GTEx whole blood as exposure, same GWAS ----
            try:
                qi = {c: i for i, c in enumerate(EQ_COLS)}
                qrows, _ = T.query(GTEX_BLOOD, chrom, lo_p, hi_p, sess=sess)
                ensg = coords[g]["ensembl_id"]
                qrows = [r for r in qrows if r[qi["gene_id"]] == ensg]
                gm = []
                for r in qrows:
                    pos = int(r[qi["position"]])
                    for x in gmap.get(pos, []):
                        try:
                            be, se = float(r[qi["beta"]]), float(r[qi["se"]])
                            so = float(x[gi["standard_error"]])
                            bo = (math.log(float(x[gi["odds_ratio"]])) if use_or
                                  else float(x[gi["beta"]]))
                        except (ValueError, KeyError):
                            continue
                        if se <= 0 or so <= 0:
                            continue
                        h = harmonize_one(r[qi["ref"]], r[qi["alt"]], x[gi["effect_allele"]],
                                          x[gi["other_allele"]], bo)
                        if h["sign"] is None:
                            continue
                        gm.append((be, se, h["beta_out_harmonized"], so, r[qi["rsid"]]))
                        break
                if len(gm) >= 30:
                    zg = np.array([m[0] / m[1] for m in gm])
                    ag = abf([m[0] / m[1] for m in gm], [m[1] for m in gm],
                             [m[2] / m[3] for m in gm], [m[3] for m in gm])
                    j = int(np.argmax(np.abs(zg)))
                    rec.update(n_merged_gtex=len(gm),
                               abf_PP4_gtexblood="%.6g" % ag["PP4"],
                               gtex_lead_z="%.4g" % zg[j],
                               n_lead_shared_gtex_eqtlgen=int(gm[j][4] == kept[i1]["rsid"]))
                else:
                    rec["note"] = (rec["note"] + ";" if rec["note"] else "") + \
                        "gtex_too_few:%d" % len(gm)
            except Exception as e:
                rec["note"] = (rec["note"] + ";" if rec["note"] else "") + \
                    "gtex_fail:%s" % type(e).__name__
        except StopIteration:
            pass
        except Exception as e:
            rec["note"] = (rec["note"] + ";" if rec["note"] else "") + \
                "fail:%s" % type(e).__name__
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerow(rec)
        done.add((g, acc))
        print("  %-10s %-14s abf=%-10s susie=%-10s gtex=%-10s  %.0fs"
              % (g, acc, rec["abf_PP4_eqtlgen"] or "-", rec["susie_PP4_eqtlgen"] or "-",
                 rec["abf_PP4_gtexblood"] or "-", time.time() - t0), flush=True)

    rows = []
    for f in sorted(os.listdir(CKPT)):
        rows += list(csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                    delimiter="\t"))
    with open(OUT + r"\part3_results.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("wrote %d rows" % len(rows))


if __name__ == "__main__":
    main()
