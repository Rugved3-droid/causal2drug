#!/usr/bin/env python
"""PART 3 R14.4 -- the head-to-head, run WITHOUT LD.

coloc.abf uses only beta and SE, so the pre-registered refutation test does not need the
1000 Genomes LD fetch that dominates the runtime of part3_degradation.py. Same loci, same
outcome GWAS, same window, same harmonization; only the exposure resource differs:

    GTEx whole blood  QTD000356  n = 853        (low power)
    eQTLGen phase I              n = 31,684     (high power)

PRE-REGISTERED (PROTOCOL.md, Amendment 1): if PP4 INCREASES with exposure power, the
Class C power-dependence claim is REFUTED and reported as such.
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
CKPT = OUT + r"\ckpt_h2h"
RAW = r"D:\causal2drug\data\raw"
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
GTEX = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
        "QTS000015/QTD000356/QTD000356.all.tsv.gz")
FLANK = 100_000
W = 0.15 ** 2
EQ = ("molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se "
      "type ac an r2 molecular_trait_object_id gene_id median_tpm rsid").split()
COLS = ["gene", "gwas_acc", "mesh_term", "n_eqtlgen", "n_gtex",
        "eqtlgen_lead_z", "gtex_lead_z", "PP4_eqtlgen", "PP3_eqtlgen",
        "PP4_gtexblood", "PP3_gtexblood", "gap_eqtlgen", "gap_gtex", "note"]


def labf(z, se):
    v = se * se
    r = W / (W + v)
    return 0.5 * (math.log(1 - r) + r * z * z)


def abf(pairs):
    return SC.coloc_from_lbf(np.array([labf(a / b, b) for a, b, _c, _d in pairs]),
                             np.array([labf(c / d, d) for _a, _b, c, d in pairs]))


def load_af(rs):
    af = {}
    with gzip.open(RAW + r"\eqtlgen_SNP_AF.txt.gz", "rt", encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_s, i_f, i_a, i_b = (h.index("SNP"), h.index("AlleleB_all"),
                              h.index("AlleleA"), h.index("AlleleB"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_s] in rs:
                try:
                    af[f[i_s]] = (float(f[i_f]), f[i_a], f[i_b])
                except ValueError:
                    pass
    return af


def smr(z, a, n):
    v = 2.0 * a * (1.0 - a) * (n + z * z)
    return (z / math.sqrt(v), 1.0 / math.sqrt(v)) if v > 0 else (None, None)


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
            if r["SNPChr"] == c["chrom"] and \
               int(c["start"]) - FLANK <= int(r["SNPPos"]) <= int(c["end"]) + FLANK:
                hits[s].append(r)
    af = load_af({r["SNP"] for v in hits.values() for r in v})
    lo = LiftOver(RAW + r"\hg19ToHg38.over.chain.gz")
    print("loci %d ; AF loaded %d" % (len(loci), len(af)), flush=True)

    done = set()
    for f in os.listdir(CKPT):
        for r in csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"), delimiter="\t"):
            done.add((r["gene"], r["gwas_acc"]))
    print("resuming: %d done" % len(done), flush=True)

    t0 = time.time()
    for loc in loci:
        g, acc = loc["gene"], loc["gwas_acc"]
        if (g, acc) in done:
            continue
        rec = {c: "" for c in COLS}
        rec.update(gene=g, gwas_acc=acc, mesh_term=loc["mesh_term"])
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
                b, se = smr(float(r["Zscore"]), fa, float(r["NrSamples"]))
                if b is None:
                    continue
                res = lo.convert_coordinate("chr" + r["SNPChr"], int(r["SNPPos"]) - 1)
                if not res:
                    continue
                ev.append(dict(pos=res[0][1] + 1, chrom=res[0][0].replace("chr", ""),
                               assessed=r["AssessedAllele"], other=r["OtherAllele"],
                               af=fa, beta=b, se=se))
            if not ev:
                raise StopIteration("no_exposure")
            chrom = ev[0]["chrom"]
            lo_p, hi_p = min(e["pos"] for e in ev) - 1000, max(e["pos"] for e in ev) + 1000

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

            def outcome(x):
                so = float(x[gi["standard_error"]])
                bo = (math.log(float(x[gi["odds_ratio"]])) if use_or
                      else float(x[gi["beta"]]))
                return bo, so

            m1 = []
            for e in ev:
                for x in gmap.get(e["pos"], []):
                    try:
                        bo, so = outcome(x)
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
                    m1.append((e["beta"], e["se"], h["beta_out_harmonized"], so))
                    break
            if len(m1) < 30:
                raise StopIteration("eqtlgen_too_few:%d" % len(m1))

            qi = {c: i for i, c in enumerate(EQ)}
            qrows, _ = T.query(GTEX, chrom, lo_p, hi_p, sess=sess)
            ensg = coords[g]["ensembl_id"]
            qrows = [r for r in qrows if r[qi["gene_id"]] == ensg]
            m2 = []
            for r in qrows:
                for x in gmap.get(int(r[qi["position"]]), []):
                    try:
                        be, se = float(r[qi["beta"]]), float(r[qi["se"]])
                        bo, so = outcome(x)
                    except (ValueError, KeyError):
                        continue
                    if se <= 0 or so <= 0:
                        continue
                    h = harmonize_one(r[qi["ref"]], r[qi["alt"]], x[gi["effect_allele"]],
                                      x[gi["other_allele"]], bo)
                    if h["sign"] is None:
                        continue
                    m2.append((be, se, h["beta_out_harmonized"], so))
                    break
            if len(m2) < 30:
                raise StopIteration("gtex_too_few:%d" % len(m2))

            a1, a2 = abf(m1), abf(m2)
            rec.update(n_eqtlgen=len(m1), n_gtex=len(m2),
                       eqtlgen_lead_z="%.4g" % max(abs(p[0] / p[1]) for p in m1),
                       gtex_lead_z="%.4g" % max(abs(p[0] / p[1]) for p in m2),
                       PP4_eqtlgen="%.6g" % a1["PP4"], PP3_eqtlgen="%.6g" % a1["PP3"],
                       PP4_gtexblood="%.6g" % a2["PP4"], PP3_gtexblood="%.6g" % a2["PP3"],
                       gap_eqtlgen="%.3f" % a1["gap"], gap_gtex="%.3f" % a2["gap"])
        except StopIteration as e:
            rec["note"] = str(e)
        except Exception as e:
            rec["note"] = "fail:%s" % type(e).__name__
        with open(os.path.join(CKPT, "%s__%s.tsv" % (g, acc)), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerow(rec)
        print("  %-10s %-14s eqtlgen=%-11s gtex=%-11s %s  %.0fs"
              % (g, acc, rec["PP4_eqtlgen"] or "-", rec["PP4_gtexblood"] or "-",
                 rec["note"], time.time() - t0), flush=True)

    rows = []
    for f in sorted(os.listdir(CKPT)):
        rows += list(csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                    delimiter="\t"))
    with open(OUT + r"\part3_headtohead.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("wrote %d rows" % len(rows))


if __name__ == "__main__":
    main()
