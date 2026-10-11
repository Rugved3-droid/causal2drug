#!/usr/bin/env python
"""Display items built from the recovery map.

    python -I manuscript_prep/manuscript/recovery_scripts/build_recovery_map.py > ...summary.txt
    python manuscript_prep/manuscript/recovery_scripts/make_recovery_outputs.py

Reads   manuscript/recovery_map.csv, manuscript/recovery_out/*.csv,
        analysis/out/per_target_primary.csv
Writes  manuscript/display/table2_named_targets.{csv,md}
        manuscript/display/table3_recovery_tiers.{csv,md}
        manuscript/display/tableS13..S17 (csv, md), manuscript/display/recovery_tables.xlsx
        manuscript/display/figure3_recovery_map.{pdf,svg,png}
"""
import collections
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.patches import Patch                             # noqa: E402
import openpyxl                                                  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
OUT = os.path.join(MS, "display")
os.makedirs(OUT, exist_ok=True)

M = list(csv.DictReader(open(os.path.join(MS, "recovery_map.csv"), encoding="utf-8")))
P = {r["gene"]: r for r in csv.DictReader(
    open(os.path.join(MP, "analysis", "out", "per_target_primary.csv"), encoding="utf-8"))}
NAMED = ["AGTR1", "ANGPTL3", "APOB", "APOC3", "AVPR2", "F2", "PCSK9", "PLG", "REN", "SLC5A1",
         "IL1B", "ITGA2B", "NPC1L1"]
PHASES = ["Preclinical", "Phase I", "Phase II", "Phase III", "Launched"]
AREAS = [k[len("area_phase: "):] for k in next(iter(P.values())) if k.startswith("area_phase: ")]
TIERS = ["candidate: CV-relevant tissue or plasma protein", "candidate: other tissue only", "locked: never assayed",
         "locked: assayed, no instrument"]
# Display labels (2026-10-11 revision). The tier KEYS in recovery_map.csv are unchanged; only the
# wording shown to readers changed, so that a tier is not read as a validated instrument or
# as drug-target rescue. "candidate: CV-relevant tissue or plasma protein" -> candidate in a CV-relevant tissue or plasma
# protein; "candidate: other tissue only" -> candidate in another tissue only.
TLAB = ["Candidate: CV-relevant tissue or plasma protein", "Candidate: other tissue only",
        "Locked: never assayed", "Locked: assayed, no instrument"]
TSHORT = ["CV-relevant or plasma-protein candidate", "Other-tissue candidate only",
          "Locked: never assayed", "Locked: assayed, no instrument"]
REV = os.path.join(MS, "revision_sources")
DRUG = list(csv.DictReader(open(os.path.join(REV, "drug_status.csv"), encoding="utf-8")))
for _d in DRUG:
    # Final reconciliation (2026-10-11): the FDA Purple Book lists ReoPro (abciximab) as
    # discontinued with the licence voluntarily revoked. That is the more specific source for
    # a biologic and supersedes the "Prescription" status in the openFDA Drugs@FDA feed.
    if _d["drug"] == "Abciximab":
        _d["fda_class"] = "discontinued"
        _d["fda_marketing_status"] += ("; FDA Purple Book (September 2026): discontinued, "
                                       "licence voluntarily revoked")
GTEXN = {r["tissueSiteDetailId"]: int(r["eqtl_samples"]) for r in csv.DictReader(
    open(os.path.join(REV, "gtex_v8_tissue_sample_sizes.csv"), encoding="utf-8"))}
# Site of drug action, taken from the mechanism-of-action text of a regulatory label or
# assessment. Source per row is in Supplementary Table S18 (column label_source) or named here.
SITE = {
    "AGTR1": "Angiotensin II type 1 receptor: vasoconstriction and aldosterone secretion by the "
             "adrenal cortex (FDA label, losartan)",
    "ANGPTL3": "ANGPTL3, expressed primarily in the liver; the antibody binds the circulating "
               "protein (FDA label, evinacumab)",
    "APOB": "Cells that produce apolipoprotein B (antisense blocks its production); organ not "
            "stated (EMA summary, Kynamro)",
    "APOC3": "Production of apolipoprotein C-III (antisense); organ not stated (EMA summary, "
             "Waylivra)",
    "AVPR2": "Vasopressin V2 receptor; effect is renal water excretion (FDA label, tolvaptan)",
    "F2": "Circulating thrombin, at its active site (FDA label, argatroban)",
    "PCSK9": "Circulating PCSK9; LDL receptors on hepatocytes (FDA label, alirocumab)",
    "PLG": "Not available from the sources checked",
    "REN": "Renin, secreted by the kidney, acting on circulating angiotensinogen (FDA label, "
           "aliskiren)",
    "SLC5A1": "SGLT1 in the intestine (and SGLT2 in the kidney); label states the mechanism of "
              "cardiovascular benefit is not established (FDA label, sotagliflozin)",
    "IL1B": "Soluble IL-1\u03b1 and IL-1\u03b2, bound by a decoy receptor (FDA label, rilonacept)",
    "ITGA2B": "Platelet glycoprotein IIb/IIIa (FDA label, eptifibatide)",
    "NPC1L1": "Brush border of the small intestine (FDA label, ezetimibe)",
}


# Angiotensin II (Giapreza) is linked to AGTR1 in ChEMBL as an AGONIST approved for shock. It
# acts in the opposite direction to the receptor blockers and is not listed as an approved
# drug for the target; AGTR1's verification rests on the receptor blockers.
NOT_LISTED = {("AGTR1", "angiotensin ii")}


def nice(name):
    n = name.lower().strip()
    return {"inclisiran sodium": "inclisiran"}.get(n, n)


def reg_status(g):
    ds = [d for d in DRUG if d["target"] == g and (g, d["drug"].lower()) not in NOT_LISTED]
    ds = [dict(d, drug=nice(d["drug"])) for d in ds]
    mk, dc, other = [], [], []
    for d in ds:
        ema_auth = ": Authorised" in d["ema_status"]
        if d["fda_class"] == "marketed":
            mk.append(d["drug"])
        elif ema_auth:
            mk.append(d["drug"] + " (EMA only)")
        elif d["fda_class"] == "discontinued":
            dc.append(d["drug"] + (" (FDA discontinued; EMA refused)" if "Refused" in
                                   d["ema_status"] else " (FDA discontinued)"))
        else:
            est = sorted({e.split(": ", 1)[1] for e in d["ema_status"].split("; ")
                          if ": " in e})
            other.append(d["drug"] + (" (no FDA record; EMA: %s)" % ", ".join(est).lower()
                                      if est else " (no FDA or EMA central record)"))
    head = ("Marketed: %d of %d drugs" % (len(mk), len(ds)) if mk else
            "No currently marketed drug found" if dc else
            "Current marketing not confirmed in either source")
    parts = [head]
    if dc:
        parts.append("discontinued: " + ", ".join(dc))
    if other:
        parts.append("not confirmed: " + ", ".join(other))
    return "; ".join(parts), bool(mk)


CATMAP = {"candidate-tissue": "candidate: tissue", "candidate-protein": "candidate: protein",
          "candidate-both": "candidate: tissue and protein", "locked": "locked"}


def disp(k, v):
    """Reader-facing wording for tier and category values; the CSV keeps the original keys."""
    if k == "reporting_tier":
        return TLAB[TIERS.index(v)]
    if k in ("category", "category_before_id_amendment"):
        return CATMAP.get(v, v)
    return v


MARK = {"IL1B": "\u2020", "APOB": "\u2021", "PLG": "\u00a7"}
CV = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Aorta", "Artery_Coronary",
      "Artery_Tibial", "Liver", "Kidney_Cortex", "Adrenal_Gland"]
STATUS = {"NO_EQTL_DATA": "No eQTL data", "NO_QUALIFYING_INSTRUMENT": "No qualifying instrument"}
TABLES = []


def tis(t):
    return t.replace("_", " ").replace("Cells Cultured fibroblasts", "Cultured fibroblasts")


def npct(n, d):
    return "%d (%.1f%%)" % (n, 100.0 * n / d) if d else "0"


def pq(x, k):
    c = x[k + "_cell"]
    if c == "instrument":
        return "%s, p = %s" % (x[k + "_lead_variant"].split(" (")[0], x[k + "_p"])
    return {"not assayed": "Not assayed", "assayed; no cis-pQTL in table": "Assayed, none listed",
            "cis-pQTL reported; outside +-100 kb": "Listed, outside ±100 kb",
            "cis-pQTL reported; distance UNVERIFIED": "Listed, distance UNVERIFIED"}[c]


# ------------------------------------------------------------------ Table 2: the 13 named
rows = []
for g in NAMED:
    x = next(q for q in M if q["gene"] == g)
    ts = [t.rsplit(" (", 1) for t in
          x["gtex_tissues_with_instrument_excl_blood_with_p"].split("; ") if t]
    tl = "; ".join("%s%s (%s" % (tis(t), "*" if t in CV else "", p) for t, p in ts) or "None"
    sub = ("X-linked, not analysed" if x["chromosome"] == "X" else
           "autosomal, not tested" if x["eqtlgen_status"] == "NO_EQTL_DATA" else
           x["eqtlgen_subtype"].split(" ", 1)[1])
    rs_txt, _ = reg_status(g)
    rows.append([g + (" (X)" if x["chromosome"] == "X" else "") + MARK.get(g, ""),
                 "; ".join(nice(d_) for d_ in P[g]["chembl_strict_drugs"].split("; ")
                           if (g, d_.lower()) not in NOT_LISTED), SITE[g], rs_txt,
                 "%s: %s" % (STATUS[x["eqtlgen_status"]], sub),
                 TSHORT[TIERS.index(x["reporting_tier"])], len(ts), tl,
                 "Yes" if x["gtex_whole_blood_instrument"] == "True" else "No",
                 pq(x, "ukbppp"), pq(x, "decode"), pq(x, "fenland")])
TABLES.append(("table2_named_targets", "Table 2",
               "Table 2. Candidate instrument availability for the 13 launched targets with "
               "an approved cardiovascular drug in ChEMBL and no eQTLGen whole-blood instrument",
               ["Targets are launched in the primary set and pass the strict drug-verification "
                "rule; drug links and approval records were verified (Table S18). "
                "Angiotensin II, an AGTR1 agonist approved for shock, is linked to AGTR1 in "
                "ChEMBL and is not listed. "
                "Instrument: a variant at p < 5e-8 within the gene body +/- 100 kb in the "
                "resource's own build. GTEx v8 tissues other than whole blood are listed in "
                "order of lead p-value (in brackets); * marks the eight pre-specified "
                "cardiovascular-relevant tissues. pQTL columns give the lead qualifying cis "
                "variant in the published table of each study. (X) X-linked gene.",
                "Tiers describe candidate instrument availability, not validated instruments. "
                "CV-relevant or plasma-protein candidate: a qualifying variant in a "
                "cardiovascular-relevant tissue or a plasma cis-pQTL. Other-tissue candidate "
                "only: a qualifying variant only in other tissues. Locked: none in any checked "
                "resource; this is relative to the resources checked and is not a statement "
                "that no instrument exists.",
                "Site of drug action is taken from the mechanism-of-action text of an FDA label "
                "or an EMA public assessment summary; it was not available for PLG. Regulatory "
                "status: FDA marketing status from Drugs@FDA (openFDA) and status in the EMA "
                "medicines table, both accessed October 2026; drug-by-drug records and sources "
                "are in Supplementary Table S18. 'Marketed' means at least one product with "
                "FDA status Prescription, or EMA status Authorised; nationally authorised "
                "European products are not in the EMA table.",
                "\u2020 IL1B: the approved drug, rilonacept, binds both IL-1\u03b1 and "
                "IL-1\u03b2 (FDA label), and its approved cardiovascular indication in ChEMBL "
                "is pericarditis, whereas the source data label the target launched for "
                "atherosclerosis. IL1B is kept in this table and is not used as a named "
                "example elsewhere.",
                "\u2021 APOB: the only linked drug, mipomersen (Kynamro), is listed as "
                "discontinued by the FDA and was refused marketing authorisation by the EMA "
                "in 2013. It was therefore approved at one time (by the FDA) and is no longer "
                "marketed; APOB remains a named example on that basis.",
                "\u00a7 PLG: neither linked drug (anistreplase, streptokinase) has a record in "
                "Drugs@FDA or the EMA medicines table, so approval by either agency could not "
                "be confirmed from these sources. PLG is kept in this table and is not used as "
                "a named example elsewhere.",
                "P values below 1\u00d710\u207b\u00b3\u2070\u2070 are derived from the \u2212log10 P "
                "reported in the source table, rounded to the nearest power of 10.",
                "'Approved' in this study means ever approved by the FDA or the EMA for a "
                "cardiovascular indication; the cardiovascular indication is taken from "
                "ChEMBL and the approval record from the two agencies' databases."],
               ["Target", "Approved cardiovascular drug(s) in ChEMBL", "Site of drug action",
                "Regulatory status (FDA / EMA)", "eQTLGen status",
                "Candidate-instrument tier", "GTEx tissues with instrument, n",
                "GTEx tissues with instrument (lead P)", "GTEx whole-blood instrument",
                "UKB-PPP cis-pQTL", "deCODE cis-pQTL", "Fenland cis-pQTL"], rows))

# ------------------------------------------------------------------ Table 3: tiers


def trow(label, rs):
    c = collections.Counter(x["reporting_tier"] for x in rs)
    n = len(rs)
    lk = c[TIERS[2]] + c[TIERS[3]]
    return [label, n] + [npct(c[t], n) for t in TIERS[:2]] + [npct(lk, n)] + \
           [npct(c[t], n) for t in TIERS[2:]]


def groups():
    g = [("All uninstrumentable targets", M), ("By eQTLGen status", None),
         ("  No eQTL data: autosomal, not tested",
          [x for x in M if x["eqtlgen_status"] == "NO_EQTL_DATA"
           and x["chromosome"] not in ("X", "MT")]),
         ("  No eQTL data: X-linked, not analysed", [x for x in M if x["chromosome"] == "X"]),
         ("  No eQTL data: mitochondrial, not analysed",
          [x for x in M if x["chromosome"] == "MT"]),
         ("  No qualifying instrument",
          [x for x in M if x["eqtlgen_status"] == "NO_QUALIFYING_INSTRUMENT"]),
         ("By highest clinical phase", None)]
    g += [("  " + p, [x for x in M if x["phase_primary"] == p]) for p in PHASES]
    g += [("By indication area", None)]
    g += [("  " + a, [x for x in M if P[x["gene"]]["area_phase: " + a]]) for a in AREAS]
    return g


r3 = [trow(l, rs) if rs is not None else [l] + [""] * 6 for l, rs in groups()]
TABLES.append(("table3_recovery_tiers", "Table 3",
               "Table 3. Recovery map of candidate instrument availability for the %d "
               "cardiovascular targets with no eQTLGen whole-blood instrument" % len(M),
               ["Resources checked: GTEx v8 cis-eQTLs in 48 tissues other than whole blood, and "
                "the published cis-pQTL tables of UKB-PPP, deCODE and Fenland. A candidate "
                "instrument is a variant at p < 5e-8 within the gene body +/- 100 kb; its "
                "validity as an instrument was not tested. Candidate, CV-relevant tissue or "
                "plasma protein: a candidate in one of eight pre-specified cardiovascular-relevant "
                "tissues (heart left ventricle, atrial appendage, aorta, coronary artery, tibial "
                "artery, liver, kidney cortex, adrenal gland) or a plasma cis-pQTL. "
                "Candidate, other tissue only: a candidate only in other tissues. Locked: no "
                "candidate in any checked resource, split by whether any of the three "
                "proteomic studies measured the protein. The tiers were defined after the "
                "first results were seen.",
                "A target is counted in every indication area in which it has a programme, so "
                "area rows do not sum to the total. Percentages are of the row.",
                "Never-assayed counts include microRNA and mitochondrial genes, which no plasma "
                "protein panel can measure (9 and 2 overall; 4 and 1 among heart failure and "
                "cardiomyopathy targets), so they exceed the protein-coding counts given in the "
                "text."],
               ["", "Targets", TLAB[0], TLAB[1], "Locked, total",
                "Locked: never assayed", "Locked: assayed, no instrument"], r3))

# ------------------------------------------------------------------ supplementary tables
keep = ["gene", "ensembl_id", "gtex_ensembl_id", "chromosome", "eqtlgen_status",
        "eqtlgen_subtype", "phase_primary", "indication_areas", "chembl_verified_strict",
        "reporting_tier", "category", "category_before_id_amendment", "flags",
        "gtex_n_tissues_tested", "gtex_n_tissues_with_instrument_excl_blood",
        "gtex_tissues_with_instrument_excl_blood_with_p",
        "gtex_cv_relevant_tissues_with_instrument", "gtex_whole_blood_instrument",
        "ukbppp_cell", "ukbppp_lead_variant", "ukbppp_p", "ukbppp_bp_from_gene_body",
        "decode_cell", "decode_lead_variant", "decode_p", "decode_bp_from_gene_body",
        "fenland_cell", "fenland_lead_variant", "fenland_p", "fenland_bp_from_gene_body"]
TABLES.append(("tableS13_recovery_map_per_target", "S13 recovery map",
               "Supplementary Table S13. Recovery map, one row per target (%d targets)" % len(M),
               ["Full file with every column: recovery_map.csv. P values in brackets are "
                "for the lead qualifying variant in that tissue."],
               keep, [[disp(k, x[k]) for k in keep] for x in M]))
ot = [x for x in M if x["reporting_tier"] == TIERS[1]]
TABLES.append(("tableS14_other_tissue_only_candidates", "S14 off-target rescues",
               "Supplementary Table S14. Targets whose only candidate instrument is in a tissue "
               "outside the eight cardiovascular-relevant tissues "
               "(%d targets)" % len(ot), [],
               ["Target", "Highest phase", "Indication areas", "eQTLGen status", "Tissues, n",
                "Tissues with a candidate instrument (lead p; GTEx v8 eQTL sample size)",
                "Proteomic studies that assayed it"],
               [[x["gene"], x["phase_primary"], x["indication_areas"],
                 STATUS[x["eqtlgen_status"]], x["gtex_n_tissues_with_instrument_excl_blood"],
                 "; ".join("%s; n=%d)" % (tis(t)[:-1], GTEXN[t.rsplit(" (", 1)[0]]) for t in
                           x["gtex_tissues_with_instrument_excl_blood_with_p"].split("; ")),
                 x["pqtl_studies_assayed"] or "none"] for x in ot]))
lk = [x for x in M if x["reporting_tier"].startswith("locked")]
TABLES.append(("tableS15_locked_targets", "S15 locked",
               "Supplementary Table S15. Locked targets: no instrument in GTEx v8 non-blood "
               "tissues or in the three published plasma cis-pQTL tables (%d targets)" % len(lk),
               ["Locked is relative to the resources checked."],
               ["Target", "Chromosome", "Highest phase", "Indication areas", "eQTLGen status",
                "Tier", "GTEx tissues in which the gene was tested, n",
                "Proteomic studies that assayed it", "Flags"],
               [[x["gene"], x["chromosome"], x["phase_primary"], x["indication_areas"],
                 STATUS[x["eqtlgen_status"]], disp("reporting_tier", x["reporting_tier"]),
                 x["gtex_n_tissues_tested"],
                 x["pqtl_studies_assayed"] or "none", x["flags"]] for x in lk]))
r = []
ALLT = {t.rsplit(" (", 1)[0] for x in M
        for t in x["gtex_tissues_with_instrument_excl_blood_with_p"].split("; ") if t} | set(CV)
PQN = {"ukbppp": "54,219", "decode": "35,559", "fenland": "10,708"}
for t in CV:
    r.append(["GTEx v8: " + tis(t), GTEXN[t],
              sum(1 for x in M if t in x.get("gtex_cv_tissues_tested", "").split("; ")),
              sum(1 for x in M
                  if t in x["gtex_cv_relevant_tissues_with_instrument"].split("; "))])
r.append(["GTEx v8: any of the eight cardiovascular-relevant tissues", "", "",
          sum(1 for x in M if x["gtex_cv_tissue_fixable"] == "True")])
r.append(["GTEx v8: any tissue other than whole blood",
          "%d to %d" % (min(v for k, v in GTEXN.items() if k != "Whole_Blood" and k in ALLT),
                        max(v for k, v in GTEXN.items() if k != "Whole_Blood" and k in ALLT)),
          sum(1 for x in M if int(x["gtex_n_tissues_tested"]) > 0),
          sum(1 for x in M if int(x["gtex_n_tissues_with_instrument_excl_blood"]) > 0)])
r.append(["GTEx v8: whole blood (does not count toward a tier)", GTEXN["Whole_Blood"],
          sum(1 for x in M if x["gtex_whole_blood_tested"] == "True"),
          sum(1 for x in M if x["gtex_whole_blood_instrument"] == "True")])
for k, nm in (("ukbppp", "UKB-PPP (Olink)"), ("decode", "deCODE (SomaScan)"),
              ("fenland", "Fenland (SomaScan)")):
    r.append(["Plasma cis-pQTL: " + nm, PQN[k], sum(1 for x in M if x[k + "_assayed"] == "True"),
              sum(1 for x in M if x[k + "_instrument"] == "True")])
r.append(["Plasma cis-pQTL: any of the three studies", "",
          sum(1 for x in M if x["pqtl_studies_assayed"]),
          sum(1 for x in M if x["pqtl_studies_with_instrument"])])
TABLES.append(("tableS16_recovery_by_resource", "S16 by resource",
               "Supplementary Table S16. Recovery by resource (%d targets)" % len(M),
               ["Tested or assayed: the gene is in the tissue's GTEx egenes file, or the "
                "protein is in the study's list of measured proteins. GTEx sample size is "
                "the number of samples in the v8 eQTL analysis, from the GTEx Portal "
                "(tissueSiteDetail, dataset gtex_v8, accessed 2026-10-11); for 'any tissue' the "
                "range is over the tissues that supply a candidate for at least one target. "
                "Proteomic sample sizes are from the abstracts of the cited papers and the "
                "Fenland workbook."],
               ["Resource", "Sample size", "Targets tested or assayed",
                "Targets with a candidate instrument"], r))
idc = list(csv.DictReader(open(os.path.join(MS, "recovery_out", "gtex_id_change_resolution.csv"),
                               encoding="utf-8")))
bm = {x["gene"]: x for x in M}
TABLES.append(("tableS17_ensembl_id_changes", "S17 ID changes",
               "Supplementary Table S17. Targets absent from GTEx v8 under their GRCh37 Ensembl "
               "gene ID, and the post-hoc symbol match", [
                   "A symbol match was allowed only where HGNC records a different current "
                   "Ensembl ID for the same symbol, GTEx carries that ID under the same symbol, "
                   "the chromosome agrees, and at least 95% of the GTEx gene body lies inside "
                   "the current Ensembl GRCh38 record."],
               ["Target", "GRCh37 Ensembl ID", "HGNC current Ensembl ID", "New ID in GTEx v8",
                "GRCh37 coordinates", "GTEx (GENCODE v26) coordinates",
                "Ensembl GRCh38 coordinates", "Share of GTEx body inside GRCh38 record",
                "Same strand as GRCh37", "Match applied", "Category before", "Category after"],
               [[x["gene"], x["grch37_ensembl_id"], x["hgnc_current_ensembl_id"] or "unchanged",
                 x["new_id_in_gtex_v8"], x["grch37_coordinates"], x["gtex_v26_coordinates"],
                 x["ensembl_grch38_coordinates"],
                 x["share_of_gtex_gene_body_inside_ensembl_grch38_record"], x["same_strand"],
                 x["symbol_match_applied"],
                 CATMAP[bm[x["gene"]]["category_before_id_amendment"]],
                 CATMAP[bm[x["gene"]]["category"]]] for x in idc]))

TABLES.append(("tableS18_regulatory_status", "S18 regulatory status",
               "Supplementary Table S18. Regulatory status and labelled mechanism of each drug "
               "linked to the 13 targets of Table 2",
               ["FDA: marketing status of products in Drugs@FDA, through the openFDA endpoint. "
                "EMA: status in the EMA medicines table (centrally authorised products only), "
                "matched on the substance fields. Accessed 2026-10-11. 'no record' means the "
                "drug is absent from that source, not that it was never approved. Label source: "
                "the DailyMed set identifier of the FDA label from which the mechanism text was "
                "read. Angiotensin II is an AGTR1 agonist and is shown for completeness only."],
               ["Target", "Drug", "FDA marketing status (products)", "FDA class",
                "EMA medicines table", "FDA label (brand)", "Label source"],
               [[d["target"], nice(d["drug"]), d["fda_marketing_status"], d["fda_class"],
                 d["ema_status"], d["label_brand"], d["label_source"]] for d in DRUG]))
import json                                                        # noqa: E402
RC = json.load(open(os.path.join(MS, "rescue_controls", "rescue_controls_results.json")))


def f_(v, fmt):
    return "" if v is None or v == "" else fmt % v


TABLES.append(("tableS19_candidate_instrument_controls", "S19 controls",
               "Supplementary Table S19. Positive controls using candidate instruments from "
               "the recovery map",
               ["Protocol written before the outcome data were retrieved "
                "(RESCUE_CONTROLS_PROTOCOL.md). Single-variant Wald ratio; expected sign "
                "positive for all five (higher target level, higher outcome). Pass: instrument "
                "p < 5e-8 within +/- 100 kb, F >= 10, MR p < 0.05 and the expected sign. Not "
                "harmonizable: the recorded variant could not be aligned under the pipeline's "
                "rule; no other variant was substituted."],
               ["Control", "Resource", "Variant", "Effect allele", "Exposure beta (SE)", "F",
                "Outcome study (N)", "Outcome beta for the effect allele (SE)", "Harmonization",
                "MR beta (95% CI)", "MR p", "Observed sign", "Verdict", "Reason"],
               [[c["control"].replace("_", " \u2192 "), c["resource"], c["variant"],
                 c.get("effect_allele", ""),
                 ("%+.4f (%.4f)" % (c["beta_exposure"], c["se_exposure"]))
                 if "beta_exposure" in c else "", f_(c.get("F"), "%.0f"),
                 "%s (%s)" % (c.get("outcome_study", ""), "{:,}".format(c["outcome_n"])
                              if "outcome_n" in c else ""),
                 ("%+.5f (%.5f)" % (c["beta_outcome_harmonized"], c["se_outcome"]))
                 if "beta_outcome_harmonized" in c else "", c.get("harmonization", ""),
                 ("%+.3f (%+.3f to %+.3f)" % (c["mr_beta"], c["mr_ci_low"], c["mr_ci_high"]))
                 if "mr_beta" in c else "",
                 ("<1e-300" if c.get("mr_p") == 0 else f_(c.get("mr_p"), "%.2g")),
                 c.get("observed_sign", ""), c["verdict"], c.get("reason", "")] for c in RC]))
TABLES.append(("tableS20_gtex_v8_sample_sizes", "S20 GTEx n",
               "Supplementary Table S20. GTEx v8 eQTL sample size by tissue, and the number of "
               "targets with a candidate instrument in each",
               ["Sample size: GTEx Portal, tissueSiteDetail, dataset gtex_v8, eQTL sample "
                "summary, accessed 2026-10-11. * cardiovascular-relevant tissue (pre-specified)."],
               ["Tissue", "eQTL samples", "Targets with a candidate instrument (of %d)" % len(M)],
               [[tis(t) + ("*" if t in CV else ""), n,
                 sum(1 for x in M if t in [q.rsplit(" (", 1)[0] for q in
                     (x["gtex_tissues_with_instrument_excl_blood_with_p"].split("; ")
                      + ([("Whole_Blood (x")] if x["gtex_whole_blood_instrument"] == "True"
                         else []))])]
                for t, n in sorted(GTEXN.items(), key=lambda kv: -kv[1])
                if n > 0]))

import re as _re                                                   # noqa: E402

_SUP = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def style(v):
    """House style for anything shown to a reader: one notation for P values and windows,
    no ISO dates. Applied to titles, headers, cells and footnotes of every table."""
    if not isinstance(v, str):
        return v
    v = _re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)e-(\d+)\b",
                lambda m: "%s\u00d710%s" % (m.group(1),
                                           ("-" + str(int(m.group(2)))).translate(_SUP)), v)
    v = v.replace("p < 5\u00d710", "P<5\u00d710").replace("p<5\u00d710", "P<5\u00d710")
    v = v.replace("+/- 100 kb", "\u00b1100 kb").replace("+/-100 kb", "\u00b1100 kb")
    v = v.replace("+-100 kb", "\u00b1100 kb").replace("+- 100 kb", "\u00b1100 kb")
    v = v.replace("F >= 10", "F\u226510").replace("MR p < 0.05", "MR P<0.05")
    v = v.replace("FDR < 0.05", "FDR<0.05").replace(", p = ", ", P=").replace("lead p ", "lead P ")
    v = v.replace("accessed 2026-10-11", "accessed October 2026")
    v = v.replace("Accessed 2026-10-11", "Accessed October 2026")
    v = v.replace("on 2026-10-11", "in October 2026").replace("2026-10-11", "October 2026")
    return v

TABLES = [(a_, b_, style(c_), [style(x_) for x_ in d_], [style(x_) for x_ in e_],
           [[style(x_) for x_ in r_] for r_ in f_]) for a_, b_, c_, d_, e_, f_ in TABLES]
wb = openpyxl.Workbook()
wb.remove(wb.active)
for stem, sheet, title, foot, header, rows in TABLES:
    with open(os.path.join(OUT, stem + ".csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    md = ["**%s**" % title, "", "| " + " | ".join(str(h) for h in header) + " |",
          "|" + "---|" * len(header)]
    md += ["| " + " | ".join(str(c).replace("|", "/") for c in r_) + " |" for r_ in rows]
    md += [""] + foot
    open(os.path.join(OUT, stem + ".md"), "w", encoding="utf-8", newline="\n").write(
        "\n".join(md) + "\n")
    ws = wb.create_sheet(sheet[:31])
    ws.append([title])
    ws.append(header)
    for r_ in rows:
        ws.append(list(r_))
    for f in foot:
        ws.append([f])
    print("wrote %s (%d rows)" % (stem, len(rows)))
wb.save(os.path.join(OUT, "recovery_tables.xlsx"))

# ------------------------------------------------------------------ Figure 3
C = ["#2a78d6", "#1baf7a", "#eb6834", "#f6b79d"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "axes.labelcolor": INK, "svg.fonttype": "none",
    "pdf.fonttype": 42, "figure.facecolor": SURF, "savefig.facecolor": SURF,
    "hatch.linewidth": 0.6, "hatch.color": "#a8502b"})


def panel(ax, grp, title):
    n_ = len(grp)
    for i, (lab, rs) in enumerate(grp):
        y = n_ - 1 - i
        c = collections.Counter(x["reporting_tier"] for x in rs)
        n = float(len(rs))
        left = 0.0
        for j, t in enumerate(TIERS):
            v = 100.0 * c[t] / n
            if v > 0:
                ax.barh(y, v, left=left, height=0.64, color=C[j], edgecolor=SURF, linewidth=1.4,
                        hatch="/////" if j == 3 else None)
                if j == 3:
                    ax.barh(y, v, left=left, height=0.64, color="none", edgecolor=SURF,
                            linewidth=1.4)
                if v >= 5.5:
                    ax.text(left + v / 2, y, "%d" % c[t], ha="center", va="center",
                            fontsize=6.9, color=INK if j == 3 else "white")
            left += v
        ax.text(101.5, y, "n = %d" % n, ha="left", va="center", fontsize=7.2, color=INK2)
    ax.set_yticks(range(n_))
    ax.set_yticklabels([g[0] for g in grp][::-1], fontsize=7.6, color=INK)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, n_ - 0.4)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title(title, fontsize=8.4, fontweight="bold", loc="left", pad=5)


ga = [("All", M)] + [(p, [x for x in M if x["phase_primary"] == p]) for p in PHASES]
gb = sorted([(a, [x for x in M if P[x["gene"]]["area_phase: " + a]]) for a in AREAS],
            key=lambda g: -sum(1 for x in g[1] if x["reporting_tier"].startswith("locked"))
            / float(len(g[1])))
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 6.3),
                             gridspec_kw=dict(height_ratios=[6, 10], hspace=0.32))
panel(a1, ga, "a  By highest clinical phase")
panel(a2, gb, "b  By indication area (ordered by share locked)")
a1.text(0, 1.02, "", transform=a1.transAxes)
a2.set_xlabel("Share of targets with no eQTLGen whole-blood instrument (%)")
hnd = [Patch(fc=C[0]), Patch(fc=C[1]), Patch(fc=C[2]),
       Patch(fc=C[3], hatch="/////", ec="#a8502b", lw=0)]
fig.legend(hnd, TLAB, loc="lower center", bbox_to_anchor=(0.5, 0.905), ncol=2, frameon=False,
           fontsize=7.2, handlelength=1.3, columnspacing=1.6, handletextpad=0.5)
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, "figure3_recovery_map." + ext), dpi=300, bbox_inches="tight")
print("wrote figure3_recovery_map; area order: %s" % [g[0] for g in gb])
