#!/usr/bin/env python
"""FIX 2 + FIX 3 -- re-specified descriptive measurement under PREREGISTRATION R15.

Stage 1 (cheap): screen every one-per-target candidate for an outcome GWAS signal at
                 p < 5e-8 inside the cis-window -- the criterion R14.1 omitted.
Stage 2 (costly): for qualifying loci only, fetch LD and compute the DENTIST residual-z
                 RMSE (Chen et al. 2021) and credible-set availability, with the FIX-1
                 ColocArtefact assertion active.

No cap: every qualifying locus is run to completion, so nothing is left UNVERIFIED.
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
import part3_degradation as P3         # noqa: E402

D1B = r"D:\causal2drug\deliverables\d1b"
OUT = r"D:\causal2drug\deliverables\d1j"
CKPT = OUT + r"\ckpt"
RAW = r"D:\causal2drug\data\raw"
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
FLANK = 100_000
ZSIG = 5.4513          # |z| for p < 5e-8
COLS = ["gene", "gwas_acc", "mesh_term", "n_merged", "eqtl_lead_rsid", "eqtl_lead_z",
        "gwas_lead_rsid", "gwas_lead_z", "gwas_genome_wide_sig", "lead_r2", "lead_dist_bp",
        "abf_gap", "abf_PP3", "abf_PP4", "susie_PP4", "n_cs_exposure", "n_cs_outcome",
        "dentist_rmse_eqtlgen", "dentist_rmse_gwas", "note"]


def main():
    os.makedirs(CKPT, exist_ok=True)
    part_c = [r for r in csv.DictReader(
        open(r"D:\causal2drug\deliverables\d1f\partc_results.tsv", encoding="utf-8"),
        delimiter="\t") if r["concordant"] in ("YES", "NO")]
    cand = {}
    for r in part_c:
        if int(r["n_overlap"]) < 100:
            continue
        g = r["gene"]
        if g not in cand or int(r["n_overlap"]) > int(cand[g]["n_overlap"]):
            cand[g] = r
    cands = sorted(cand.values(), key=lambda r: r["gene"])
    print("R15 candidates (n_overlap>=100, one per target): %d" % len(cands), flush=True)

    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D1B + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    import re
    harm = {}
    for l in open(RAW + r"\harmonised_list.txt"):
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m:
            harm[m.group(1)] = l.strip().lstrip("./")

    genes = {r["gene"] for r in cands}
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
    af = P3.load_af({r["SNP"] for v in hits.values() for r in v})
    lo = LiftOver(RAW + r"\hg19ToHg38.over.chain.gz")
    keep = P3.eur()
    print("AF loaded for %d rsIDs" % len(af), flush=True)

    done = {}
    for f in os.listdir(CKPT):
        for r in csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"), delimiter="\t"):
            done[(r["gene"], r["gwas_acc"])] = r
    print("resuming: %d loci done" % len(done), flush=True)

    t0 = time.time()
    for loc in cands:
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
                b, se = P3.smr(float(r["Zscore"]), fa, float(r["NrSamples"]))
                if b is None:
                    continue
                res = lo.convert_coordinate("chr" + r["SNPChr"], int(r["SNPPos"]) - 1)
                if not res:
                    continue
                ev.append(dict(rsid=r["SNP"], chrom=res[0][0].replace("chr", ""),
                               pos38=res[0][1] + 1, assessed=r["AssessedAllele"],
                               other=r["OtherAllele"], af=fa, beta=b, se=se))
            if not ev:
                raise StopIteration("no_exposure")
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
                raise StopIteration("too_few_merged:%d" % len(merged))
            z2all = np.array([m["bo"] / m["so"] for m in merged])
            z1all = np.array([m["beta"] / m["se"] for m in merged])
            i2 = int(np.argmax(np.abs(z2all)))
            i1 = int(np.argmax(np.abs(z1all)))
            sig = bool(abs(z2all[i2]) > ZSIG)
            rec.update(n_merged=len(merged), gwas_genome_wide_sig=str(sig),
                       eqtl_lead_rsid=merged[i1]["rsid"], eqtl_lead_z="%.4g" % z1all[i1],
                       gwas_lead_rsid=merged[i2]["rsid"], gwas_lead_z="%.4g" % z2all[i2])
            if not sig:
                rec["note"] = "excluded_R15_gwas_not_genome_wide_significant"
                raise StopIteration(rec["note"])

            # ---- stage 2: LD, only for qualifying loci ----
            lines2, _ = T.header_lines(P3.KG % chrom, sess=sess, nbytes=1 << 18)
            cl = next(l for l in lines2 if l.startswith("#CHROM"))
            samples = cl.rstrip("\r").split("\t")[9:]
            idx = np.array([i for i, s in enumerate(samples) if s in keep])
            recs, _ = T.query(P3.KG % chrom, "chr" + chrom,
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
                raise StopIteration("too_few_with_LD:%d" % len(kept))
            X = np.vstack(vecs)
            R = (X @ X.T) / X.shape[1]
            np.fill_diagonal(R, 1.0)
            z1 = np.array([m["beta"] / m["se"] for m in kept])
            se1 = np.array([m["se"] for m in kept])
            z2 = np.array([m["bo"] / m["so"] for m in kept])
            se2 = np.array([m["so"] for m in kept])
            a = P3.abf(z1, se1, z2, se2)
            j1, j2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
            rec.update(n_merged=len(kept), lead_r2="%.4f" % (R[j1, j2] ** 2),
                       lead_dist_bp=abs(kept[j1]["pos38"] - kept[j2]["pos38"]),
                       abf_gap="%.3f" % a["gap"], abf_PP3="%.6g" % a["PP3"],
                       abf_PP4="%.6g" % a["PP4"],
                       dentist_rmse_eqtlgen="%.4g" % P3.dentist_rmse(z1, R),
                       dentist_rmse_gwas="%.4g" % P3.dentist_rmse(z2, R))
            try:
                s = SC.coloc_susie(z1, z2, R, L=5, strict=True)
                rec.update(susie_PP4=("%.6g" % s["best"]["PP4"]) if s["best"] else "",
                           n_cs_exposure=s["n_cs1"], n_cs_outcome=s["n_cs2"])
                if not s["best"]:
                    rec["note"] = "no_credible_set_pair"
            except SC.ColocArtefact as ce:
                rec["note"] = "coloc_artefact_rejected"
                print("      REJECTED: %s" % ce, flush=True)
        except StopIteration as e:
            if not rec["note"]:
                rec["note"] = str(e)
        except Exception as e:
            rec["note"] = "fail:%s" % type(e).__name__
        with open(os.path.join(CKPT, "%s__%s.tsv" % (g, acc)), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerow(rec)
        print("  %-10s %-14s sig=%-6s abf=%-11s susie=%-10s %s  %.0fs"
              % (g, acc, rec["gwas_genome_wide_sig"], rec["abf_PP4"] or "-",
                 rec["susie_PP4"] or "-", rec["note"][:34], time.time() - t0), flush=True)

    rows = []
    for f in sorted(os.listdir(CKPT)):
        rows += list(csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                    delimiter="\t"))
    with open(OUT + r"\r15_results.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("wrote %d rows" % len(rows))


if __name__ == "__main__":
    main()
