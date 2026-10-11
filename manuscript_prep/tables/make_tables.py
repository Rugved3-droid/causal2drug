#!/usr/bin/env python
"""Manuscript tables from the corrected classification.

    python manuscript_prep/analysis/classify.py     # first
    python manuscript_prep/tables/make_tables.py

Reads only manuscript_prep/analysis/out/ and the frozen pair table. Writes, for each
table, a CSV and a formatted Markdown file, plus one formatted workbook
(manuscript_tables.xlsx) with every table on its own sheet.
"""
import collections
import csv
import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
AOUT = os.path.join(MP, "analysis", "out")
sys.path.insert(0, os.path.join(MP, "analysis"))

K = json.load(open(os.path.join(AOUT, "key_numbers.json")))
prim = list(csv.DictReader(open(os.path.join(AOUT, "per_target_primary.csv"), encoding="utf-8")))
PHASES = ["Preclinical", "Phase I", "Phase II", "Phase III", "Launched"]
SETS = ["PRIMARY", "SENS1", "SENS2", "C14ONLY"]
SETNAME = {"PRIMARY": "Primary set", "SENS1": "Sensitivity 1 (original 630)",
           "SENS2": "Sensitivity 2", "C14ONLY": "Sensitivity 3 (MeSH C14 only)"}
AREAS = ["Lipids / atherosclerosis", "Ischaemic heart disease", "Heart failure / cardiomyopathy",
         "Arrhythmia", "Hypertension", "Pulmonary hypertension", "Thrombosis",
         "Cerebrovascular", "Other", "Unspecified CV"]
LAB = {"INSTRUMENTED": "Instrumented", "NO_EQTL_DATA": "No eQTL data",
       "NO_QUALIFYING_INSTRUMENT": "No qualifying instrument", "UNRESOLVED": "Unresolved"}
TABLES = []      # (file stem, sheet name, title, footnotes, header, rows)


def npct(n, d):
    return "%d (%.1f%%)" % (n, 100.0 * n / d) if d else "0"


def trow(label, t):
    return [label, t["n"], npct(t["instrumented"], t["n"]), npct(t["no_eqtl_data"], t["n"]),
            t["no_eqtl_data_autosomal"],
            "%d / %d" % (t["no_eqtl_data_chrX"], t["no_eqtl_data_chrMT"]),
            "%d / %d (%.1f%%)" % (t["autosomal_lower_bound"], t["n_not_on_x_y_mt"],
                                  100.0 * t["autosomal_lower_bound"] / t["n_not_on_x_y_mt"])
            if t["n_not_on_x_y_mt"] else "",
            npct(t["no_qualifying_instrument"], t["n"]), npct(t["upper_bound"], t["n"]),
            t["unresolved"]]


BH = ["", "Targets", "Instrumented", "No eQTL data (lower bound)",
      "...autosomal gene not tested", "...X-linked / mitochondrial, not analysed",
      "Autosomal-only lower bound", "No qualifying instrument",
      "Uninstrumentable, upper bound", "Unresolved gene record"]
PROV = ("The evidence rows on which the two ChEMBL rules disagree were reviewed and the "
        "strict rule was retained. Ever approved by FDA or EMA: the strict rule, and at "
        "least one of the target's strict-rule drugs has an approval record in Drugs@FDA or "
        "the EMA medicines table; this is the definition of 'approved' "
        "used in the text. Two of the 33 rest on an indication that is cardiovascular only in a "
        "broad sense: NR3C1 (dexamethasone, prednisolone and prednisone; ChEMBL approved "
        "indication myocarditis) and MS4A1 (rituximab; systemic vasculitis). Both are "
        "instrumented, so they enter the denominator only. Strict: an approved drug linked directly to the gene "
        "(single protein, specific complex or mRNA) has an approved indication that sits only "
        "in cardiovascular MeSH branches. Permissive: any ChEMBL link, any approved "
        "indication under MeSH C14 or a lipid term.")
FOOT_DEF = ("Instrumented: at least one eQTLGen whole-blood cis-eQTL at p < 5e-8 within the "
            "gene body +/- 100 kb. No eQTL data: the gene is not among the 19,250 genes eQTLGen "
            "tested; its sub-types are autosomal genes not tested, and X-linked and "
            "mitochondrial genes, which eQTLGen did not analyse (no target is Y-linked). The autosomal-only lower bound is autosomal "
            "no-eQTL-data targets over targets not on X, Y or MT. No qualifying instrument: "
            "tested, but no variant meets the criterion. The upper bound is no eQTL data plus "
            "no qualifying instrument and is reported only as a bound.")

# ------------------------------------------------------------------ Table 1
f, s = K["flow"], K["sets"]
r1 = [["Rows in the source target-indication table", f["raw_rows"], "", "", ""],
      ["Distinct target-indication pairs", f["distinct_ti_uid"], "", "", ""],
      ["Pairs after removing 4 programmes with conflicting phase labels", f["frozen_pairs"],
       "", "", ""],
      ["Distinct targets / indications in those pairs",
       "%d / %d" % (f["targets"], f["indications"]), "", "", ""],
      ["", SETNAME["PRIMARY"], SETNAME["SENS1"], SETNAME["SENS2"], SETNAME["C14ONLY"]],
      ["Cardiovascular indications (MeSH terms)"] + [s[k]["terms"] for k in SETS],
      ["Cardiovascular target-indication pairs"] + [s[k]["pairs"] for k in SETS],
      ["Cardiovascular targets"] + [s[k]["targets"] for k in SETS]]
for key, lab in (("instrumented", "Instrumented"),
                 ("no_eqtl_data", "No eQTL data (lower bound)"),
                 ("no_eqtl_data_autosomal", "  autosomal gene not tested"),
                 ("no_eqtl_data_chrX", "  X-linked, not analysed by eQTLGen"),
                 ("no_eqtl_data_chrMT", "  mitochondrial, not analysed by eQTLGen"),
                 ("no_qualifying_instrument", "No qualifying instrument"),
                 ("upper_bound", "Uninstrumentable, upper bound"),
                 ("unresolved", "Unresolved gene record")):
    r1.append(["  " + lab] + [npct(K["overall"][k][key], K["overall"][k]["n"]) for k in SETS])
r1.append(["  Autosomal-only lower bound"] + [
    "%d / %d (%.1f%%)" % (K["overall"][k]["autosomal_lower_bound"],
                          K["overall"][k]["n_not_on_x_y_mt"],
                          K["overall"][k]["autosomal_lower_pct"]) for k in SETS])
r1.append(["Highest phase reached in the set"] + [""] * 4)
for ph in PHASES:
    r1.append(["  " + ph] + [K["by_phase"][k][ph]["n"] for k in SETS])
r1.append(["Launched, approved CV drug verified in ChEMBL: strict rule"]
          + [K["launched"][k]["strict"]["tally"]["n"] for k in SETS])
r1.append(["Launched, drug ever approved by FDA or EMA for a CV indication"]
          + [K["launched"][k]["approved"]["tally"]["n"] for k in SETS])
r1.append(["Launched, approved CV drug verified in ChEMBL: permissive rule"]
          + [K["launched"][k]["permissive"]["tally"]["n"] for k in SETS])
TABLES.append(("table1_flow_and_cohort", "Table 1", "Table 1. Derivation of the cardiovascular "
               "target set and its composition",
               ["Primary set: every indication that Minikel et al. tag as cardiovascular, plus "
                "indications matching eight lipid/atherosclerosis tokens. Sensitivity 1: the "
                "originally analysed set (MeSH branch C14 plus the lipid tokens). Sensitivity 2: "
                "sensitivity 1 plus cardiovascular-tagged terms filed under signs and symptoms. "
                "Sensitivity 3: MeSH branch C14 only.",
                "Phase is the highest Minikel phase over the target's indications in the same "
                "set. Percentages are of all targets in the set.", FOOT_DEF, PROV],
               ["", "N", "", "", ""], r1))

# ------------------------------------------------------------------ Table 2
r2 = [trow(ph, K["by_phase"]["PRIMARY"][ph]) for ph in PHASES]
r2.append(trow("All targets", K["overall"]["PRIMARY"]))
r2.append(trow("Launched, drug ever approved by FDA or EMA for a CV indication "
               "", K["launched"]["PRIMARY"]["approved"]["tally"]))
r2.append(trow("Launched, approved CV drug verified: strict",
               K["launched"]["PRIMARY"]["strict"]["tally"]))
r2.append(trow("Launched, approved CV drug verified: permissive",
               K["launched"]["PRIMARY"]["permissive"]["tally"]))
r2.append(trow("Phase III or launched, drug ever approved by FDA or EMA",
               K["chembl_verified_any_label_primary"]["approved"]))
r2.append(trow("Phase III or launched, verified: strict",
               K["chembl_verified_any_label_primary"]["strict"]))
r2.append(trow("Phase III or launched, verified: permissive",
               K["chembl_verified_any_label_primary"]["permissive"]))
TABLES.append(("table2_status_by_phase", "Table 2", "Table 2. Instrument status by highest "
               "clinical phase, primary cardiovascular set (%d targets)" % s["PRIMARY"]["targets"],
               [FOOT_DEF, PROV, "Minikel phase labels were not altered; the last two rows add "
                "Phase III targets for which ChEMBL lists such a drug."], BH, r2))

# ------------------------------------------------------------------ Table 3
PHN = {p: i + 1 for i, p in enumerate(PHASES)}


def tal(rows):
    c = collections.Counter(r["status"] for r in rows)
    xl = sum(1 for r in rows if r["no_eqtl_data_subtype"]
             and not r["no_eqtl_data_subtype"].startswith("autosomal"))
    sx = sum(1 for r in rows if r["no_eqtl_data_subtype"].startswith("X-linked"))
    sm = sum(1 for r in rows if r["no_eqtl_data_subtype"].startswith("mitochondrial"))
    nx = sum(1 for r in rows if r["chromosome"] in ("X", "Y", "MT"))
    return dict(n=len(rows), instrumented=c["INSTRUMENTED"], no_eqtl_data=c["NO_EQTL_DATA"],
                no_eqtl_data_x_linked=xl, no_eqtl_data_autosomal=c["NO_EQTL_DATA"] - xl,
                no_eqtl_data_chrX=sx, no_eqtl_data_chrMT=sm,
                autosomal_lower_bound=c["NO_EQTL_DATA"] - xl, n_not_on_x_y_mt=len(rows) - nx,
                no_qualifying_instrument=c["NO_QUALIFYING_INSTRUMENT"],
                upper_bound=c["NO_EQTL_DATA"] + c["NO_QUALIFYING_INSTRUMENT"],
                unresolved=c["UNRESOLVED"])


r3 = []
for a in AREAS:
    col = "area_phase: " + a
    inn = [r for r in prim if r[col]]
    L = [r for r in inn if r[col] == "Launched"]
    LS = [r for r in L if r["chembl_verified_strict"] == "True"]
    t, tl, ts = tal(inn), tal(L), tal(LS)
    r3.append(trow(a, t) + [tl["n"], tl["no_eqtl_data"], tl["upper_bound"],
                            ts["n"], ts["no_eqtl_data"], ts["upper_bound"]])
t = tal(prim)
LL = [r for r in prim if r["phase_PRIMARY"] == "Launched"]
LS = [r for r in LL if r["chembl_verified_strict"] == "True"]
r3.append(trow("Any area (each target once)", t)
          + [len(LL), tal(LL)["no_eqtl_data"], tal(LL)["upper_bound"],
             len(LS), tal(LS)["no_eqtl_data"], tal(LS)["upper_bound"]])
TABLES.append(("table3_status_by_indication_area", "Table 3", "Table 3. Instrument status by "
               "cardiovascular indication area, primary set",
               ["A target is counted once in every area in which it has an indication, so rows "
                "do not sum to the total. 'Unspecified CV' holds the single generic MeSH term "
                "'Cardiovascular Diseases'. Launched in area: the target carries the Launched "
                "label for an indication in that area.", FOOT_DEF,
                "Assignment of MeSH terms to areas is in Supplementary Table S2."],
               ["Indication area"] + BH[1:] + ["Launched in area", "...no eQTL data",
                                              "...upper bound",
                                              "Launched in area, strict-verified",
                                              "...no eQTL data", "...upper bound"], r3))

# ------------------------------------------------------------------ S1 per-target
keep = [c for c in prim[0] if not c.startswith("area_phase: ")]
TABLES.append(("tableS1_per_target_classification", "S1 per-target",
               "Supplementary Table S1. Classification of every target in the primary "
               "cardiovascular set (%d targets)" % len(prim),
               ["One row per target. status uses the +/-100 kb window; status_window_500kb and "
                "status_window_1Mb are the window sensitivity analyses. Columns ending "
                "_upper_bound use the eQTL Catalogue top cis p-value over +/-1 Mb, unclumped.",
                "phase_<SET> is blank when the target is not in that set."],
               keep, [[r[c] for c in keep] for r in sorted(prim, key=lambda r: r["gene"])]))

# ------------------------------------------------------------------ S2 CV set terms
ct = list(csv.DictReader(open(os.path.join(AOUT, "cv_terms_primary.csv"), encoding="utf-8")))
TABLES.append(("tableS2_cv_indication_terms", "S2 CV terms",
               "Supplementary Table S2. The %d MeSH indications in the primary cardiovascular "
               "set, with the rule that admits each and its indication area" % len(ct),
               ["R1: Minikel top-level branch is C14. R2: term matches a lipid/atherosclerosis "
                "token. R3: Minikel tags the term cardiovascular and files it under signs and "
                "symptoms (C23). R4: tagged cardiovascular, filed under another branch.",
                "flag marks terms that are in the set by rule but are not primarily cardiac; "
                "none was removed. Sensitivity 2 contains no R4 term."],
               list(ct[0].keys()), [list(x.values()) for x in ct]))

# ------------------------------------------------------------------ S3 sensitivity sets
r = []
for k in SETS:
    r.append([SETNAME[k] + " - all targets"] + trow("", K["overall"][k])[1:])
    for ph in PHASES:
        r.append([SETNAME[k] + " - " + ph] + trow("", K["by_phase"][k][ph])[1:])
    r.append([SETNAME[k] + " - launched, drug ever approved by FDA or EMA"]
             + trow("", K["launched"][k]["approved"]["tally"])[1:])
    r.append([SETNAME[k] + " - launched, verified strict"]
             + trow("", K["launched"][k]["strict"]["tally"])[1:])
    r.append([SETNAME[k] + " - launched, verified permissive"]
             + trow("", K["launched"][k]["permissive"]["tally"])[1:])
r.append(["All 2,517 Minikel targets (any therapeutic area)"]
         + trow("", K["overall"]["ALL_2517"])[1:])
TABLES.append(("tableS3_sensitivity_cv_set", "S3 set sensitivity",
               "Supplementary Table S3. Instrument status under each definition of the "
               "cardiovascular set", [FOOT_DEF, PROV], BH, r))

# ------------------------------------------------------------------ S4 window
r = []
for k in ("PRIMARY", "SENS1", "SENS2"):
    for w in ("100kb", "500kb", "1Mb"):
        r.append(["%s, +/-%s, all targets" % (SETNAME[k], w)]
                 + trow("", K["window"][k][w]["overall"])[1:])
        r.append(["%s, +/-%s, launched" % (SETNAME[k], w)]
                 + trow("", K["window"][k][w]["launched"])[1:])
g = K["window"]["primary_targets_gaining_instrument"]
TABLES.append(("tableS4_sensitivity_cis_window", "S4 window sensitivity",
               "Supplementary Table S4. Sensitivity of instrument status to the cis-window",
               ["+/-100 kb is the pre-registered primary window. 'No eQTL data' cannot change "
                "with the window. Targets gaining an instrument in a wider window (primary "
                "set): %s." % ("; ".join("%s (%s, nearest variant %s bp from the gene body)"
                                         % (x["gene"], x["phase"], x["nearest_bp"]) for x in g)
                               or "none"),
                "eQTLGen tested variants within 1 Mb of the gene centre, so a +/-1 Mb window "
                "measured from the gene body is not completely covered; the largest distance "
                "observed for any target is %s bp." % K["window"]["max_distance_observed_bp"]],
               BH, r))

# ------------------------------------------------------------------ S5 sub-types
sub = [r for r in prim if r["status"] == "NO_QUALIFYING_INSTRUMENT"]
h5 = ["gene", "phase_PRIMARY", "no_instrument_subtype", "gtex_blood_tpm",
      "n_cis_eqtl_p5e8_any_distance", "nearest_p5e8_variant_bp_from_gene_body",
      "best_p_any_distance", "tissue_restricted", "chembl_verified_strict",
      "chembl_verified_permissive"]
c5 = collections.Counter(r["no_instrument_subtype"] for r in sub)
TABLES.append(("tableS5_no_qualifying_instrument_subtypes", "S5 sub-types",
               "Supplementary Table S5. Targets with no qualifying instrument, by sub-type "
               "(primary set, %d targets)" % len(sub),
               ["Counts: " + "; ".join("%s: %d" % kv for kv in sorted(c5.items())) + ".",
                "T1: a variant at p < 5e-8 exists but lies outside the +/-100 kb window. T2: "
                "eQTLGen reports a cis-eQTL at FDR < 0.05 for the gene, none at p < 5e-8. T3: "
                "no cis-eQTL at FDR < 0.05."], h5,
               [[r[c] for c in h5] for r in sorted(sub, key=lambda r: (r["no_instrument_subtype"],
                                                                       r["gene"]))]))

# ------------------------------------------------------------------ S6 mismatches
mm = list(csv.DictReader(open(os.path.join(AOUT, "label_chembl_mismatches.csv"),
                              encoding="utf-8")))
TABLES.append(("tableS6_phase_label_vs_chembl", "S6 label vs ChEMBL",
               "Supplementary Table S6. Targets whose Minikel phase label and ChEMBL approval "
               "record disagree (primary set, %d targets)" % len(mm),
               ["Minikel labels were used as published and were not edited. ChEMBL 37. Launched "
                "targets are listed when they are not strict-verified; Phase III targets when "
                "ChEMBL lists an approved drug with an approved cardiovascular indication under "
                "the permissive rule.", PROV,
                "'No ChEMBL mechanism' means ChEMBL links no drug to the gene, which is a limit "
                "of the link and not evidence that no drug exists."],
               list(mm[0].keys()), [list(m.values()) for m in mm]))

# ------------------------------------------------------------------ S7 launched detail
L = [r for r in prim if r["phase_PRIMARY"] == "Launched"]
h7 = ["gene", "status", "no_eqtl_data_subtype", "no_instrument_subtype",
      "chembl_verified_strict", "chembl_verified_permissive", "example_eligible",
      "chembl_strict_fail_reason", "cv_indications_at_highest_phase", "chembl_strict_drugs",
      "chembl_approved_cv_drugs",
      "gtex_blood_tpm", "gtex_liver_tpm", "gtex_heart_LV_tpm", "gtex_max_tissue",
      "tissue_restricted", "secreted", "in_SENS1", "in_SENS2", "phase_SENS1", "phase_SENS2"]
order = {"NO_EQTL_DATA": 0, "NO_QUALIFYING_INSTRUMENT": 1, "INSTRUMENTED": 2, "UNRESOLVED": 3}
TABLES.append(("tableS7_launched_targets", "S7 launched targets",
               "Supplementary Table S7. All %d launched cardiovascular targets in the primary "
               "set, with instrument status and ChEMBL drugs" % len(L),
               ["Sorted by status. Drug lists show up to 15 names. example_eligible: "
                "strict-verified and at the same phase in sensitivity set 2; only such targets "
                "are named as examples in the text and figures.", PROV], h7,
               [[r[c] for c in h7] for r in sorted(L, key=lambda r: (order[r["status"]],
                                                                     r["gene"]))]))

# ------------------------------------------------------------------ S8 gene records
gr = list(csv.DictReader(open(os.path.join(AOUT, "gene_record_resolution.csv"),
                              encoding="utf-8")))
TABLES.append(("tableS8_gene_record_resolution", "S8 gene records",
               "Supplementary Table S8. Gene symbols without a unique GRCh37 record, and their "
               "resolution (%d symbols; %d unresolved)"
               % (len(gr), sum(1 for x in gr if x["result"] == "UNRESOLVED")),
               ["Rule: take the Ensembl gene ID of the approved HGNC record; accept it if it is "
                "a primary-chromosome gene in Ensembl GRCh37 release 87."],
               list(gr[0].keys()), [list(x.values()) for x in gr]))

# ------------------------------------------------------------------ S9 enrichment
r = []
for k in ("PRIMARY", "SENS1", "SENS2", "ALL_2517"):
    for key, e in K["enrichment"][k].items():
        grp, fld = key.split("_vs_instrumented__")
        (a, b), (c, d) = e["table"]
        r.append([SETNAME.get(k, "All 2,517 targets"), LAB[grp], fld.replace("_", " "),
                  "%d / %d" % (a, a + b), "%d / %d" % (c, c + d), "%.2f" % e["odds_ratio"],
                  "%.2f-%.2f" % (e["ci_lo"], e["ci_hi"]), "%.3g" % e["p_fisher"]])
TABLES.append(("tableS9_enrichment", "S9 enrichment",
               "Supplementary Table S9. Tissue-restricted expression and secretion among "
               "uninstrumentable versus instrumented targets",
               ["Human Protein Atlas annotation. Tissue restricted: RNA tissue specificity "
                "'Tissue enriched', 'Group enriched' or 'Tissue enhanced' versus 'Low tissue "
                "specificity'. Odds ratio with Woolf 95% CI; Fisher exact p."],
               ["Set", "Group", "Feature", "Group with feature", "Instrumented with feature",
                "Odds ratio", "95% CI", "p"], r))


# ------------------------------------------------------------------ S10 chromosome
r = []
for k in ("PRIMARY", "SENS1", "SENS2", "C14ONLY", "ALL_2517"):
    c = K["no_eqtl_data_by_chromosome"][k]
    r.append([SETNAME.get(k, "All 2,517 targets"), c["n"], c["no_eqtl_data"], c["autosomal"],
              c["on_X"], c["on_MT"], c["on_Y"], c["targets_on_X_Y_MT_in_set"],
              c["X_with_blood_tpm_ge1"], c["MT_with_blood_tpm_ge1"],
              "%d / %d" % (c["autosomal_blood_tpm_below_1"], c["autosomal_with_tpm"]),
              "%.1f%%" % c["autosomal_lower_bound_pct"]])
TABLES.append(("tableS10_no_eqtl_data_by_chromosome", "S10 chromosome",
               "Supplementary Table S10. Targets with no eQTL data: autosomal genes versus "
               "genes on chromosomes eQTLGen did not analyse",
               ["None of the 19,250 genes eQTLGen tested lies on X, Y or the mitochondrial "
                "genome, so every such target has no eQTL data whatever its expression in "
                "blood. Autosomal lower bound: autosomal targets with no eQTL data as a share "
                "of autosomal and unresolved targets in the set."],
               ["Set", "Targets", "No eQTL data", "...autosomal", "...X-linked",
                "...mitochondrial", "...Y-linked", "All X/Y/MT targets in set",
                "X-linked with blood TPM >= 1", "Mitochondrial with blood TPM >= 1",
                "Autosomal no-eQTL-data with blood TPM < 1", "Autosomal lower bound"], r))


import re as _re                                                   # noqa: E402

_SUP = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def style(v):
    """House style for anything shown to a reader: one notation for P values and windows,
    no ISO dates. Applied to titles, headers, cells and footnotes of every table."""
    if not isinstance(v, str):
        return v
    v = _re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)e-(\d+)\b",
                lambda m: "%s\u00d710%s" % (m.group(1), ("-" + m.group(2)).translate(_SUP)), v)
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
# ================================================================== writers
def md(title, foot, header, rows, maxrows=None):
    out = ["**%s**" % title, "", "| " + " | ".join(str(h) for h in header) + " |",
           "|" + "---|" * len(header)]
    for r_ in (rows if maxrows is None else rows[:maxrows]):
        out.append("| " + " | ".join(str(x).replace("|", "/") for x in r_) + " |")
    if maxrows is not None and len(rows) > maxrows:
        out.append("")
        out.append("*First %d of %d rows shown; the CSV and workbook hold all rows.*"
                   % (maxrows, len(rows)))
    out.append("")
    for f_ in foot:
        out.append(f_)
        out.append("")
    return "\n".join(out)


from openpyxl import Workbook                                   # noqa: E402
from openpyxl.styles import Alignment, Border, Font, Side       # noqa: E402
from openpyxl.utils import get_column_letter                    # noqa: E402

wb = Workbook()
wb.remove(wb.active)
thin = Side(style="thin", color="000000")
for stem, sheet, title, foot, header, rows in TABLES:
    with open(os.path.join(HERE, stem + ".csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    big = len(rows) > 80
    open(os.path.join(HERE, stem + ".md"), "w", encoding="utf-8").write(
        md(title, foot, header, rows, 40 if big else None))
    ws = wb.create_sheet(sheet[:31])
    ws.append([title])
    ws["A1"].font = Font(bold=True, size=11)
    ws.append(list(header))
    for c in ws[2]:
        c.font = Font(bold=True)
        c.border = Border(top=thin, bottom=thin)
        c.alignment = Alignment(wrap_text=True, vertical="bottom")
    for r_ in rows:
        ws.append([("" if x is None else x) for x in r_])
    for c in ws[ws.max_row]:
        c.border = Border(bottom=thin)
    ws.append([])
    for f_ in foot:
        ws.append([f_])
        ws.cell(ws.max_row, 1).font = Font(italic=True, size=9)
    for i in range(1, len(header) + 1):
        width = max([len(str(header[i - 1]))] + [len(str(r_[i - 1])) for r_ in rows[:300]
                                                 if i - 1 < len(r_)])
        ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 46)
    ws.freeze_panes = "B3"
    print("%-46s %4d rows x %2d columns" % (stem, len(rows), len(header)))
wb.save(os.path.join(HERE, "manuscript_tables.xlsx"))
print("wrote manuscript_tables.xlsx (%d sheets)" % len(TABLES))
