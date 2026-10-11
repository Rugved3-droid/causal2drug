#!/usr/bin/env python
"""Part C: re-run the two surviving positive controls under the frozen rules.

Chain per pair:
  eQTLGen cis hits (GRCh37, p<5e-8, gene body +-100 kb)
    -> beta/se from Z via the SMR approximation, using eQTLGen's own allele frequencies
    -> LD clump vs 1000G phase 3 EUR (GRCh37), r2 < 0.01
    -> liftover GRCh37 -> GRCh38 (UCSC chain), failure rate reported
    -> GWAS Catalog harmonised outcome region (GRCh38) by tabix range query
    -> allele harmonization incl. flip / strand-flip branches (src/harmonize.py, unit-tested)
    -> MR: lead-SNP Wald ratio (PREREGISTRATION R12.2)
    -> colocalization (R8)
    -> provenance

Rules are frozen: threshold p<5e-8 (DECISION 6), window +-100 kb (R1), clumping r2<0.01
(R5), estimator lead-SNP Wald ratio (R12.2). Nothing here may be tuned to a result.
"""
import collections
import csv
import datetime
import gzip
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T          # noqa: E402
import d1b_clump as C           # noqa: E402
from harmonize import harmonize_one  # noqa: E402
from pyliftover import LiftOver      # noqa: E402

D1B = r"D:\causal2drug\deliverables\d1b"
OUT = r"D:\causal2drug\deliverables\d1c"
RAW = r"D:\causal2drug\data\raw"
CHAIN = RAW + r"\hg19ToHg38.over.chain.gz"
AF_FILE = RAW + r"\eqtlgen_SNP_AF.txt.gz"
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

FLANK = 100_000
P_INSTRUMENT = 5e-8
R2_CLUMP = 0.01
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_QUANT = 0.15 ** 2

PAIRS = [
    dict(name="HMGCR_LDL", gene="HMGCR", ensg="ENSG00000113161",
         trait="LDL cholesterol", efo="EFO_0004611", acc="GCST90239658",
         rel="GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         expect="POSITIVE"),
    dict(name="IL6R_CRP", gene="IL6R", ensg="ENSG00000160712",
         trait="C-reactive protein", efo="EFO_0004458", acc="GCST90179146",
         rel="GCST90179001-GCST90180000/GCST90179146/harmonised/GCST90179146.h.tsv.gz",
         expect="NEGATIVE"),
]


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def z_to_beta_se(z, af, n):
    """SMR approximation (Zhu et al. 2016): beta and se from Z, allele freq and N.

    Note beta/se == z exactly under this parameterisation, so F = z^2 and the instrument
    strength does not depend on the allele frequency.
    """
    v = 2.0 * af * (1.0 - af) * (n + z * z)
    if v <= 0:
        return None, None
    se = 1.0 / math.sqrt(v)
    return z * se, se


def load_af(rsids):
    af = {}
    with gzip.open(AF_FILE, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_snp = hdr.index("SNP")
        i_af = hdr.index("AlleleB_all")
        i_a = hdr.index("AlleleA")
        i_b = hdr.index("AlleleB")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_snp] in rsids:
                try:
                    af[f[i_snp]] = (float(f[i_af]), f[i_a], f[i_b])
                except ValueError:
                    pass
    return af


def hits_for(sym, coords):
    c = coords[sym]
    lo, hi = int(c["start"]) - FLANK, int(c["end"]) + FLANK
    rows = []
    with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["MinikelSymbol"] != sym or r["SNPChr"] != c["chrom"]:
                continue
            if lo <= int(r["SNPPos"]) <= hi:
                rows.append(r)
    return rows, c, lo, hi


# ---------------- coloc (R8, unchanged) ----------------
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


def coloc(rows):
    l1 = [labf(r["beta_exp"], r["se_exp"]) for r in rows]
    l2 = [labf(r["beta_out_harm"], r["se_out"]) for r in rows]
    s1, s2 = lse(l1), lse(l2)
    s12 = lse([a + b for a, b in zip(l1, l2)])
    ls = [0.0, math.log(P1) + s1, math.log(P2) + s2,
          math.log(P1) + math.log(P2) + ldiff(s1 + s2, s12), math.log(P12) + s12]
    d = lse(ls)
    pp = [math.exp(x - d) for x in ls]
    return dict(n_variants=len(rows), PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3],
                PP4=pp[4])


def norm_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def run(p, coords, keep, lo_chain, log_all):
    log = {"pair": p["name"], "started_utc": now(),
           "rules": "R1 window +-100kb; DECISION 6 p<5e-8; R5 clump r2<0.01; "
                    "R12.2 lead-SNP Wald ratio"}
    rows, c, lo, hi = hits_for(p["gene"], coords)
    log["gene_grch37"] = {"symbol": p["gene"], "ensembl_id": c["ensembl_id"],
                          "chrom": c["chrom"], "start": int(c["start"]),
                          "end": int(c["end"]), "window": [lo, hi],
                          "source": "Ensembl GRCh37 release 87 GTF"}
    log["n_eqtlgen_hits_p_lt_5e-8"] = len(rows)

    # ---- beta/se from Z ----
    afmap = log_all["af"]
    exp_rows = []
    n_no_af = 0
    for r in rows:
        a = afmap.get(r["SNP"])
        if a is None:
            n_no_af += 1
            continue
        af_b, allele_a, allele_b = a
        # AlleleB is the assessed allele in the AF file; align to eQTLGen AssessedAllele
        if r["AssessedAllele"] == allele_b:
            af_assessed = af_b
        elif r["AssessedAllele"] == allele_a:
            af_assessed = 1.0 - af_b
        else:
            n_no_af += 1
            continue
        z = float(r["Zscore"])
        n = float(r["NrSamples"])
        b, se = z_to_beta_se(z, af_assessed, n)
        if b is None:
            n_no_af += 1
            continue
        exp_rows.append(dict(rsid=r["SNP"], chrom37=r["SNPChr"], pos37=int(r["SNPPos"]),
                             assessed=r["AssessedAllele"], other=r["OtherAllele"],
                             z=z, n=n, af_assessed=af_assessed, beta_exp=b, se_exp=se,
                             p_exp=float(r["Pvalue"]), F=z * z))
    log["n_dropped_no_allele_frequency"] = n_no_af
    log["n_exposure_variants"] = len(exp_rows)

    # ---- clumping (GRCh37) ----
    ctx = C.ChromCtx(c["chrom"], keep)
    raw_for_clump = sorted(rows, key=lambda r: float(r["Pvalue"]))
    cl = C.clump_one(ctx, p["gene"], c, raw_for_clump)
    log["clumping"] = {"n_hits": cl["n_hits"], "n_in_ld_panel": cl["n_in_ld_panel"],
                       "n_independent": cl["n_independent"],
                       "independent_snps": cl["independent_snps"].split(";")}
    lead = min(exp_rows, key=lambda r: r["p_exp"])
    log["lead_snp"] = {k: lead[k] for k in
                       ("rsid", "pos37", "assessed", "other", "z", "n", "af_assessed",
                        "beta_exp", "se_exp", "p_exp", "F")}

    # ---- liftover GRCh37 -> GRCh38 ----
    n_lift_ok = n_lift_fail = n_lift_multi = 0
    fails = []
    for r in exp_rows:
        res = lo_chain.convert_coordinate("chr" + r["chrom37"], r["pos37"] - 1)
        if not res:
            n_lift_fail += 1
            r["pos38"] = None
            fails.append(r["rsid"])
            continue
        if len(res) > 1:
            n_lift_multi += 1
        ch38, p38, strand, _ = res[0]
        r["pos38"] = p38 + 1
        r["chrom38"] = ch38.replace("chr", "")
        r["lift_strand"] = strand
        n_lift_ok += 1
    log["liftover"] = {
        "chain": "UCSC hg19ToHg38.over.chain.gz",
        "sha256": "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d",
        "n_in": len(exp_rows), "n_lifted": n_lift_ok, "n_failed": n_lift_fail,
        "n_multi_mapping": n_lift_multi,
        "failure_rate": round(n_lift_fail / len(exp_rows), 6) if exp_rows else None,
        "failed_rsids": fails[:50],
        "handling": "variants that fail to lift are DROPPED and counted; never "
                    "carried forward on their GRCh37 coordinate"}
    exp_rows = [r for r in exp_rows if r["pos38"] is not None]
    neg = [r for r in exp_rows if r.get("lift_strand") == "-"]
    log["liftover"]["n_negative_strand"] = len(neg)

    # ---- outcome GWAS region (GRCh38) ----
    p38 = [r["pos38"] for r in exp_rows]
    w_lo, w_hi = min(p38) - 1000, max(p38) + 1000
    url = GW + p["rel"]
    gcols, hst = None, None
    lines, hst = T.header_lines(url, nbytes=1 << 16)
    gcols = lines[0].rstrip("\r").split("\t")
    t0 = time.time()
    grows, gstats = T.query(url, exp_rows[0]["chrom38"], w_lo, w_hi)
    gstats["seconds"] = round(time.time() - t0, 2)
    gstats["retrieved_utc"] = now()
    gstats["accession"] = p["acc"]
    gstats["columns"] = gcols
    log["gwas_retrieval"] = gstats
    log["gwas_window_grch38"] = [w_lo, w_hi]

    gi = {c2: i for i, c2 in enumerate(gcols)}
    gmap = collections.defaultdict(list)
    for g in grows:
        gmap[int(g[gi["base_pair_location"]])].append(g)

    # ---- harmonize ----
    counts = collections.Counter()
    merged = []
    for r in exp_rows:
        cands = gmap.get(r["pos38"])
        if not cands:
            counts["no_position_match"] += 1
            continue
        best = None
        for g in cands:
            try:
                b_out = float(g[gi["beta"]])
                se_out = float(g[gi["standard_error"]])
                p_out = float(g[gi["p_value"]])
            except ValueError:
                continue
            try:
                af_out = float(g[gi["effect_allele_frequency"]])
            except (ValueError, KeyError):
                af_out = None
            h = harmonize_one(r["other"], r["assessed"],
                              g[gi["effect_allele"]], g[gi["other_allele"]], b_out,
                              af_alt_exp=r["af_assessed"], af_effect_out=af_out)
            if h["sign"] is not None:
                best = (g, h, b_out, se_out, p_out, af_out)
                break
            best = best or (g, h, b_out, se_out, p_out, af_out)
        if best is None:
            counts["no_parsable_gwas_record"] += 1
            continue
        g, h, b_out, se_out, p_out, af_out = best
        counts[h["status"]] += 1
        if h["sign"] is None:
            continue
        if se_out <= 0 or r["se_exp"] <= 0:
            counts["bad_se"] += 1
            continue
        merged.append(dict(
            rsid=r["rsid"], chrom=r["chrom38"], pos37=r["pos37"], pos38=r["pos38"],
            exp_other=r["other"], exp_assessed=r["assessed"],
            gwas_effect=g[gi["effect_allele"]], gwas_other=g[gi["other_allele"]],
            harmonization=h["status"], sign_applied=h["sign"],
            ambiguous="TRUE" if h["ambiguous"] else "FALSE",
            effect_allele_harmonized=r["assessed"],
            other_allele_harmonized=r["other"],
            beta_exp=r["beta_exp"], se_exp=r["se_exp"], p_exp=r["p_exp"],
            z_exp=r["z"], n_exp=r["n"], af_assessed_exp=r["af_assessed"], F=r["F"],
            beta_out_raw=b_out, beta_out_harm=h["beta_out_harmonized"],
            se_out=se_out, p_out=p_out, af_effect_out=af_out))
    log["harmonization"] = dict(counts)
    log["n_overlapping_after_harmonization"] = len(merged)

    cols = list(merged[0].keys())
    with open(os.path.join(OUT, "%s_harmonized_merge.tsv" % p["name"]), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for m in merged:
            w.writerow([m[c2] for c2 in cols])

    # ---- MR: lead-SNP Wald ratio (R12.2) ----
    lead_m = min(merged, key=lambda m: m["p_exp"])
    b = lead_m["beta_out_harm"] / lead_m["beta_exp"]
    se = abs(lead_m["se_out"] / lead_m["beta_exp"])
    z = b / se
    log["mr"] = {
        "estimator": "lead-SNP Wald ratio (R12.2)",
        "lead_snp": lead_m["rsid"], "lead_p_exposure": lead_m["p_exp"],
        "beta_exposure": lead_m["beta_exp"], "se_exposure": lead_m["se_exp"],
        "beta_outcome_harmonized": lead_m["beta_out_harm"], "se_outcome": lead_m["se_out"],
        "harmonization_of_lead": lead_m["harmonization"],
        "beta": b, "se": se, "z": z, "p": norm_p(z),
        "direction": "POSITIVE" if b > 0 else "NEGATIVE",
        "expected_direction": p["expect"],
        "F": lead_m["F"], "n_instruments_after_clumping": cl["n_independent"]}
    log["coloc"] = coloc(merged)
    log["pass"] = {
        "instrument_p_lt_5e-8": lead_m["p_exp"] < P_INSTRUMENT,
        "F_ge_10": lead_m["F"] >= 10.0,
        "mr_p_lt_0.05": log["mr"]["p"] < 0.05,
        "direction_matches": log["mr"]["direction"] == p["expect"]}
    log["VERDICT"] = "PASS" if all(log["pass"].values()) else "FAIL"
    log["finished_utc"] = now()
    json.dump(log, open(os.path.join(OUT, "%s_provenance.json" % p["name"]), "w"),
              indent=2)
    return log


def main():
    os.makedirs(OUT, exist_ok=True)
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D1B + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    keep = C.eur_samples()
    lo_chain = LiftOver(CHAIN)

    need = set()
    for p in PAIRS:
        rows, _c, _lo, _hi = hits_for(p["gene"], coords)
        need |= {r["SNP"] for r in rows}
    print("looking up allele frequencies for %d rsIDs ..." % len(need), flush=True)
    t0 = time.time()
    af = load_af(need)
    print("  found %d/%d in eQTLGen AF file (%.0fs)" % (len(af), len(need),
                                                        time.time() - t0))
    log_all = {"af": af}

    summary = []
    for p in PAIRS:
        print("\n" + "=" * 72 + "\n" + p["name"] + "\n" + "=" * 72, flush=True)
        lg = run(p, coords, keep, lo_chain, log_all)
        print(json.dumps({k: lg[k] for k in
                          ("n_eqtlgen_hits_p_lt_5e-8", "clumping", "liftover",
                           "harmonization", "n_overlapping_after_harmonization",
                           "mr", "coloc", "pass", "VERDICT")}, indent=2))
        summary.append(lg)
    json.dump(summary, open(os.path.join(OUT, "d1c_summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
