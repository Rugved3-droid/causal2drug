#!/usr/bin/env python
"""Positive controls with candidate instruments from the recovery map.
Implements manuscript/RESCUE_CONTROLS_PROTOCOL.md (written before this was run).

    python manuscript_prep/manuscript/rescue_controls/run_rescue_controls.py

Uses, unchanged and by import: src/select_gwas.py (R2 + R2a), src/tabix_http.py,
src/harmonize.py (R6). Estimator: single-variant Wald ratio (R12.2).
Writes  rescue_controls/gwas_selection_<EFO>.csv, exposure_rows.json, outcome_rows.json,
        rescue_controls_results.json, rescue_controls_results.csv, run_log.txt
"""
import contextlib
import csv
import datetime
import gzip
import io
import json
import math
import os
import re
import sys
import tarfile

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
ROOT = os.path.dirname(MP)
sys.path.insert(0, os.path.join(ROOT, "src"))
import select_gwas as SG            # noqa: E402
import tabix_http as T              # noqa: E402
from harmonize import harmonize_one  # noqa: E402
from pyliftover import LiftOver     # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
DATA = os.path.join(MS, "recovery_data")
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
CHAIN = os.path.join(ROOT, "data", "raw", "hg19ToHg38.over.chain.gz")
P_INSTRUMENT = 5e-8
CONTROLS = [
    dict(name="PCSK9_LDL", gene="PCSK9", outcome="LDL cholesterol", efo="EFO_0004611",
         expect="POSITIVE"),
    dict(name="ANGPTL3_TG", gene="ANGPTL3", outcome="triglycerides", efo="EFO_0004530",
         expect="POSITIVE"),
    dict(name="APOC3_TG", gene="APOC3", outcome="triglycerides", efo="EFO_0004530",
         expect="POSITIVE"),
    dict(name="AGTR1_SBP", gene="AGTR1", outcome="systolic blood pressure", efo="EFO_0006335",
         expect="POSITIVE"),
    dict(name="REN_SBP", gene="REN", outcome="systolic blood pressure", efo="EFO_0006335",
         expect="POSITIVE"),
]
LOG = []


def log(s):
    LOG.append(s)
    print(s)


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def clean(v):
    return "" if v is None else str(v).strip().lstrip("'").strip()


def sheet(path, name, header_has):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    it = wb[name].iter_rows(values_only=True)
    for row in it:
        h = [clean(c).replace("\n", " ") for c in row]
        if header_has in h:
            break
    out = [dict(zip(h, row)) for row in it if any(c is not None for c in row)]
    wb.close()
    return out


def norm_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


M = {r["gene"]: r for r in csv.DictReader(open(os.path.join(MS, "recovery_map.csv"),
                                               encoding="utf-8"))}

# ------------------------------------------------------------------ exposure, per protocol


def exposure(c):
    x = M[c["gene"]]
    if x["pqtl_studies_with_instrument"]:
        if x["ukbppp_instrument"] == "True":
            sh = x["ukbppp_note"]
            rs, pos = re.match(r"(\S+) \(chr(\S+?):(\d+)\)", x["ukbppp_lead_variant"]).group(1, 3)
            pcol = "log10(p)" if sh == "ST10" else "log10(p) (discovery)"
            bcol = "BETA" if sh == "ST10" else "BETA (discovery, wrt. A1)"
            scol = "SE" if sh == "ST10" else "SE (discovery)"
            fcol = "A1FREQ" if sh == "ST10" else "A1FREQ (discovery)"
            hit = [r for r in sheet(os.path.join(DATA, "41586_2023_6592_MOESM3_ESM.xlsx"), sh,
                                    "UKBPPP ProteinID")
                   if clean(r["cis/trans"]).lower() == "cis"
                   and clean(r["UKBPPP ProteinID"]).split(":")[0] == c["gene"]
                   and str(int(float(r["GENPOS (hg38)"]))) == pos]
            assert len(hit) == 1, (c["gene"], len(hit))
            r = hit[0]
            vid = clean(r[[k for k in r if k.startswith("Variant ID")][0]]).split(":")
            return dict(resource="UKB-PPP plasma cis-pQTL (%s)" % sh, rsid=rs,
                        chrom38=clean(r["CHROM"]), pos38=int(pos), other=vid[2], effect=vid[3],
                        beta=float(r[bcol]), se=float(r[scol]), af_effect=float(r[fcol]),
                        p="10^-%s" % clean(r[pcol]), log10p=float(r[pcol]),
                        unit="SD of plasma protein level (as scaled in the UKB-PPP table)",
                        bp_from_gene_body=x["ukbppp_bp_from_gene_body"],
                        variant_id_in_table=":".join(vid))
        if x["fenland_instrument"] == "True":
            rs, ch, pos = re.match(r"(\S+) \(chr(\S+?):(\d+)\)",
                                   x["fenland_lead_variant"]).group(1, 2, 3)
            hit = [r for r in sheet(os.path.join(DATA, "pietzner2021_supplementary_tables.xlsx"),
                                    "ST2", "SomaScan.id")
                   if clean(r["cis.trans"]).lower() == "cis"
                   and c["gene"] in re.split(r"[|;,\s]+", clean(r["HGNC.symbol.protein"]))
                   and str(int(float(r["Position"]))) == pos]
            hit.sort(key=lambda r: float(r["P.value"]) if not clean(r["P.value"]).startswith("<")
                     else 0.0)
            r = hit[0]
            lo = LiftOver(CHAIN).convert_coordinate("chr" + ch, int(pos) - 1)
            if not lo:
                return dict(resource="Fenland plasma cis-pQTL (ST2)", rsid=rs,
                            not_harmonizable="liftover GRCh37 -> GRCh38 failed")
            se = float(r["SE"])
            b = float(r["Effect"])
            return dict(resource="Fenland plasma cis-pQTL (ST2)", rsid=rs,
                        chrom38=lo[0][0].replace("chr", ""), pos38=lo[0][1] + 1,
                        pos37=int(pos), lift_strand=lo[0][2], other=clean(r["NEA"]).upper(),
                        effect=clean(r["EA"]).upper(), beta=b, se=se,
                        af_effect=float(r["EAF"]), p=clean(r["P.value"]),
                        log10p=-math.log10(float(r["P.value"])) if float(r["P.value"]) > 0
                        else 300.0,
                        unit="SD of plasma protein level (as scaled in the Fenland table)",
                        bp_from_gene_body=x["fenland_bp_from_gene_body"],
                        aptamer=clean(r["SomaScan.id"]),
                        n_matching_rows=len(hit))
        return dict(resource="deCODE", not_harmonizable="deCODE table gives no standard error "
                    "or stated effect allele; protocol allows deCODE only if established")
    tissue, var = x["gtex_best_cv_tissue"], x["gtex_best_cv_variant_b38"]
    gid = x["gtex_ensembl_id"]
    ch, pos, ref, alt = var.split("_")[:4]
    with tarfile.open(os.path.join(DATA, "GTEx_Analysis_v8_eQTL.tar")) as tf:
        m = tf.getmember("GTEx_Analysis_v8_eQTL/%s.v8.signif_variant_gene_pairs.txt.gz" % tissue)
        with gzip.open(tf.extractfile(m), "rt", encoding="utf-8") as fh:
            h = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if line.startswith(var + "\t"):
                    f = dict(zip(h, line.rstrip("\n").split("\t")))
                    if f["gene_id"].split(".")[0] == gid:
                        return dict(resource="GTEx v8 %s cis-eQTL" % tissue, rsid=var,
                                    chrom38=ch.replace("chr", ""), pos38=int(pos), other=ref,
                                    effect=alt, beta=float(f["slope"]), se=float(f["slope_se"]),
                                    af_effect=None, maf=float(f["maf"]), p=f["pval_nominal"],
                                    log10p=-math.log10(float(f["pval_nominal"])),
                                    unit="normalised expression unit in GTEx %s (slope, ALT "
                                         "allele)" % tissue,
                                    bp_from_gene_body="within +-100 kb (recovery map)")
    return dict(resource="GTEx v8 %s" % tissue, rsid=var,
                not_harmonizable="variant row not found in the tissue file")


# ------------------------------------------------------------------ outcome selection (R2+R2a)
SEL = {}


def select(efo):
    if efo in SEL:
        return SEL[efo]
    out = os.path.join(HERE, "gwas_selection_%s.csv" % efo)
    if not os.path.exists(out):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            SG.main(efo, out)
        open(os.path.join(HERE, "gwas_selection_%s.log" % efo), "w",
             encoding="utf-8").write(buf.getvalue())
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    elig = [r for r in rows if r["tbi_present"] == "True" and r["european_only"] == "True"]
    SEL[efo] = (elig[0] if elig else None, len(rows), len(elig))
    return SEL[efo]


# ------------------------------------------------------------------ run
results, exp_dump, out_dump = [], {}, {}
log("run started %s" % now())
for c in CONTROLS:
    res = dict(control=c["name"], target=c["gene"], outcome=c["outcome"], efo=c["efo"],
               expected_sign=c["expect"])
    e = exposure(c)
    exp_dump[c["name"]] = e
    res.update(resource=e["resource"], variant=e.get("rsid", ""))
    sel, n_all, n_el = select(c["efo"])
    if sel is None:
        res.update(verdict="OUTCOME NOT RETRIEVABLE", reason="no eligible study under R2+R2a")
        results.append(res)
        continue
    res.update(outcome_study=sel["accession"], outcome_first_author=sel["first_author"],
               outcome_n=int(sel["max_parsed_N"]), outcome_sample=sel["initial_sample_size"],
               n_studies_with_harmonised_file=n_all, n_eligible_studies=n_el)
    if "not_harmonizable" in e:
        res.update(verdict="NOT HARMONIZABLE", reason=e["not_harmonizable"])
        results.append(res)
        continue
    res.update(chrom38=e["chrom38"], pos38=e["pos38"], effect_allele=e["effect"],
               other_allele=e["other"], beta_exposure=e["beta"], se_exposure=e["se"],
               p_exposure=e["p"], F=(e["beta"] / e["se"]) ** 2, exposure_unit=e["unit"],
               af_effect_exposure=e.get("af_effect"),
               bp_from_gene_body=e["bp_from_gene_body"])
    url = GW + sel["rel_path"].lstrip("./")
    try:
        lines, _ = T.header_lines(url, nbytes=1 << 16)
        gcols = lines[0].rstrip("\r").split("\t")
        grows, gstats = T.query(url, e["chrom38"], e["pos38"] - 1000, e["pos38"] + 1000)
    except Exception as ex:                                       # noqa: BLE001
        res.update(verdict="OUTCOME NOT RETRIEVABLE", reason=repr(ex)[:300])
        results.append(res)
        continue
    gi = {k: i for i, k in enumerate(gcols)}
    cands = [g for g in grows if int(g[gi["base_pair_location"]]) == e["pos38"]]
    out_dump[c["name"]] = dict(url=url, columns=gcols, rows_at_position=cands,
                               n_rows_in_window=len(grows), retrieved_utc=now())
    if not cands:
        res.update(verdict="NOT HARMONIZABLE",
                   reason="variant position absent from the outcome file (%d rows within "
                          "+-1 kb)" % len(grows))
        results.append(res)
        continue
    best = None
    for g in cands:
        try:
            b_out, se_out = float(g[gi["beta"]]), float(g[gi["standard_error"]])
            p_out = g[gi["p_value"]]
        except ValueError:
            continue
        try:
            af_out = float(g[gi["effect_allele_frequency"]])
        except (ValueError, KeyError):
            af_out = None
        h = harmonize_one(e["other"], e["effect"], g[gi["effect_allele"]],
                          g[gi["other_allele"]], b_out, af_alt_exp=e.get("af_effect"),
                          af_effect_out=af_out)
        cand = (g, h, b_out, se_out, p_out, af_out)
        if h["sign"] is not None:
            best = cand
            break
        best = best or cand
    if best is None:
        res.update(verdict="NOT HARMONIZABLE", reason="no parsable outcome record")
        results.append(res)
        continue
    g, h, b_out, se_out, p_out, af_out = best
    res.update(harmonization=h["status"], outcome_effect_allele=g[gi["effect_allele"]],
               outcome_other_allele=g[gi["other_allele"]], beta_outcome_raw=b_out,
               se_outcome=se_out, p_outcome=p_out, af_effect_outcome=af_out)
    if h["sign"] is None:
        res.update(verdict="NOT HARMONIZABLE", reason="harmonize_one status: " + h["status"])
        results.append(res)
        continue
    b = h["beta_out_harmonized"] / e["beta"]
    se = abs(se_out / e["beta"])
    z = b / se
    obs = "POSITIVE" if b > 0 else "NEGATIVE"
    crit = {"instrument_p_lt_5e-8": e["log10p"] > -math.log10(P_INSTRUMENT),
            "F_ge_10": res["F"] >= 10.0, "mr_p_lt_0.05": norm_p(z) < 0.05,
            "direction_matches": obs == c["expect"]}
    res.update(beta_outcome_harmonized=h["beta_out_harmonized"], mr_beta=b, mr_se=se,
               mr_ci_low=b - 1.96 * se, mr_ci_high=b + 1.96 * se, mr_z=z, mr_p=norm_p(z),
               observed_sign=obs, criteria=crit,
               verdict="PASS" if all(crit.values()) else "FAIL")
    results.append(res)

for r in results:
    log("%-11s %-34s %-18s -> %s" % (r["control"], r["resource"][:34], r["variant"][:18],
                                     r["verdict"]))
    if "mr_beta" in r:
        log("            exposure beta %+.4f (SE %.4f, F %.0f) for allele %s; outcome %s N=%s "
            "beta %+.5f (SE %.5f), harmonization %s" % (
                r["beta_exposure"], r["se_exposure"], r["F"], r["effect_allele"],
                r["outcome_study"], r["outcome_n"], r["beta_outcome_harmonized"],
                r["se_outcome"], r["harmonization"]))
        log("            MR beta %+.4f (95%% CI %+.4f to %+.4f), p %.3g; expected %s, observed "
            "%s; criteria %s" % (r["mr_beta"], r["mr_ci_low"], r["mr_ci_high"], r["mr_p"],
                                 r["expected_sign"], r["observed_sign"], r["criteria"]))
    else:
        log("            %s | outcome study %s" % (r.get("reason"), r.get("outcome_study")))
log("run finished %s" % now())
json.dump(results, open(os.path.join(HERE, "rescue_controls_results.json"), "w"), indent=1)
json.dump(exp_dump, open(os.path.join(HERE, "exposure_rows.json"), "w"), indent=1)
json.dump(out_dump, open(os.path.join(HERE, "outcome_rows.json"), "w"), indent=1)
keys = []
for r in results:
    for k in r:
        if k not in keys:
            keys.append(k)
with open(os.path.join(HERE, "rescue_controls_results.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(results)
open(os.path.join(HERE, "run_log.txt"), "w", encoding="utf-8").write("\n".join(LOG) + "\n")
