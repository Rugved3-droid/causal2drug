#!/usr/bin/env python
"""Supplementary Table S21 (phase confidence intervals) and the STROBE-MR checklist with the
official item wording.

    python manuscript_prep/analysis/phase_ci.py > manuscript_prep/analysis/phase_ci_output.txt
    python manuscript_prep/manuscript/draft/make_revision_tables.py

The checklist wording is read from the official fillable Word checklist downloaded from
strobe-mr.org on 2026-10-11 (revision_sources/strobe_mr_checklist.docx, parsed to
strobe_mr_checklist_rows.json). The checklist is copyrighted by the Equator Network under
CC BY 3.0. Only the right-hand column ("where addressed") is ours.
"""
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
OUT = os.path.join(MS, "display")

# ---------------------------------------------------------------- S11, S12
# The former main-text Tables 2 and 3 are Supplementary Tables S11 and S12. They are re-copied
# from tables/ on every run so that they cannot go stale.
for src_name, dst_name, old, new in (
        ("table2_status_by_phase.md", "tableS11_status_by_phase.md", "**Table 2.",
         "**Supplementary Table S11."),
        ("table3_status_by_indication_area.md", "tableS12_status_by_indication_area.md",
         "**Table 3.", "**Supplementary Table S12.")):
    t = open(os.path.join(MP, "tables", src_name), encoding="utf-8").read()
    assert t.startswith(old), src_name
    open(os.path.join(OUT, dst_name), "w", encoding="utf-8", newline="\n").write(
        new + t[len(old):])
print("re-copied tables S11 and S12")

# ---------------------------------------------------------------- S21
txt = open(os.path.join(MP, "analysis", "phase_ci_output.txt"), encoding="utf-8").read()
md = ["**Supplementary Table S21. Share of targets by highest clinical phase, with 95% "
      "confidence intervals**", "",
      "| Measure | Phase | Targets with the status / all targets | Share | 95% CI |",
      "|---|---|---|---|---|"]
cur = ""
for line in txt.split("\n"):
    m = re.match(r"\[P\d+\] (.*)", line)
    if m:
        cur = m.group(1)
        continue
    m = re.match(r"\s+(\S.*?)\s+(\d+) /\s+(\d+) =\s+([\d.]+)% \(95% CI\s+([\d.]+)% to\s+([\d.]+)%\)",
                 line)
    if m:
        md.append("| %s | %s | %s / %s | %s%% | %s%% to %s%% |" % (cur, m.group(1), m.group(2),
                                                                  m.group(3), m.group(4),
                                                                  m.group(5), m.group(6)))
        continue
    m = re.match(r"\s+Launched minus Preclinical: (\S+) percentage points \(95% CI (\S+) to (\S+)\)",
                 line)
    if m:
        md.append("| %s | Launched minus preclinical | | %s points | %s to %s points |"
                  % (cur, m.group(1), m.group(2), m.group(3)))
md += ["", "Wilson score intervals for proportions; Newcombe hybrid score interval for the "
       "difference between two independent proportions. Descriptive; no hypothesis about phase "
       "was pre-specified. The last block is among the 280 targets with no eQTLGen whole-blood "
       "instrument."]
open(os.path.join(OUT, "tableS21_phase_confidence_intervals.md"), "w", encoding="utf-8",
     newline="\n").write("\n".join(md) + "\n")
print("wrote tableS21 (%d rows)" % (len(md) - 6))

# ---------------------------------------------------------------- STROBE-MR
rows = json.load(open(os.path.join(MS, "revision_sources", "strobe_mr_checklist_rows.json"),
                      encoding="utf-8"))
NA = ("Not applicable to the main result, which is whether a candidate instrument exists, "
      "not an MR estimate. ")
W = {
 "1": "Title; Abstract (Methods, Results). MR is the context; the main result is instrument availability, with MR positive controls",
 "2": "Introduction, paragraphs 1 to 3",
 "3": "Introduction, last paragraph. Pre-specified expected signs for the controls: Methods, Statistics and Controls; Supplemental Methods 6 and 7",
 "4": "Methods, Target Set; Instrument Status in eQTLGen; Recovery Map of Candidate Instruments. Supplemental Methods 1 (table of data sources)",
 "4a": "Public summary-level resources; cohorts, sample sizes and versions in Methods and Supplemental Methods 1; Table S20 for GTEx tissue sample sizes",
 "4b": "Sample sizes: Methods (eQTLGen 31,684; UKB-PPP 54,219; deCODE 35,559; Fenland 10,708); Tables S16, S19, S20. No power calculation was done",
 "4c": "Methods, Instrument Status in eQTLGen; Recovery Map (P<5×10⁻⁸, gene body ±100 kb, resource's own build and significance flag); Supplemental Methods 5",
 "4d": "Exposures: gene expression (eQTLGen, GTEx v8) and plasma protein level (3 pQTL studies). Outcomes of the controls: LDL cholesterol, triglycerides, systolic blood pressure, C-reactive protein from GWAS Catalog studies named in Supplemental Methods 6 and 7 and Table S19",
 "4e": "Methods, Data and Code Availability (public de-identified summary data; approval not required; institutional confirmation pending)",
 "5": NA + "For the controls the three core assumptions are not separately stated; Discussion (paragraph on validity) lists what an available cis-QTL does not establish",
 "6": "Methods, Statistics and Controls; Controls",
 "6a": "Counts and proportions. Controls: exposure per unit as scaled in the source table (SD of protein level; normalised expression), outcome in the GWAS's units (Table S19)",
 "6b": "One variant per control: the lead qualifying variant recorded in the recovery map; no weights (Supplemental Methods 7)",
 "6c": "Single-variant Wald ratio; SE = |SE_outcome / beta_exposure| (Methods, Statistics and Controls; Supplemental Methods 6 and 7). Covariate sets of the source studies were not compared",
 "6d": "Targets with no resolvable gene record are reported as unresolved (3). Controls that could not be harmonized are reported as such and not replaced",
 "6e": "No correction; all tests are descriptive (Methods, Statistics and Controls)",
 "7": NA + "Positive controls with pre-specified signs (Methods, Statistics and Controls); IL6R audit of allele harmonization and linkage disequilibrium (Supplemental Methods 8)",
 "8": "Alternative target sets and cis-windows (Tables S3, S4); permissive drug-verification rule (Table S7); identifier-change sensitivity (Table S17); controls of candidate instruments (Table S19)",
 "9a": "Supplemental Methods 10 (Python and package versions; commands)",
 "9b": "Not pre-registered. The recovery rule and the control protocol were written in advance as dated documents without an external timestamp (Methods; Limitations; Supplemental Methods 5 and 7)",
 "10a": "Figure 1; Table 1 (targets and pairs at each step). Individual-level counts are those of the source studies",
 "10b": "Not applicable: no phenotypic data were analysed",
 "10c": "Not assessed (eQTLGen is a meta-analysis; heterogeneity statistics were not used)",
 "10d": "Controls only: exposure and outcome samples are both of European ancestry by the outcome selection rule. Sample overlap (UK Biobank contributes to UKB-PPP and to the outcome GWAS) was not quantified",
 "11a": "Table S19 (variant-exposure and variant-outcome associations for the candidate-instrument controls); Results, Candidate Instruments as Positive Controls",
 "11b": "Results, Candidate Instruments as Positive Controls; Table S19 (Wald ratio, 95% CI, P)",
 "11c": "Not applicable",
 "11d": "Figures 1 to 3 show instrument availability; no MR plot (one variant per control)",
 "12a": "Results, Candidate Instruments as Positive Controls: 1 of 3 eQTLGen controls and 2 of 5 candidate-instrument controls passed; IL6R audit (Supplemental Methods 8)",
 "12b": "Not applicable (single-variant estimates)",
 "13a": "Not applicable to single-variant estimates; none performed",
 "13b": "Results, first section; Tables S3, S4, S7, S17",
 "13c": "Not assessed",
 "13d": "Expected signs are taken from the known effects of the approved drugs (Methods, Statistics and Controls)",
 "13e": "None",
 "14": "Discussion, paragraph 1",
 "15": "Discussion, Limitations; Discussion paragraph on what an available cis-QTL does not establish",
 "16a": "Discussion, paragraphs 2 to 7",
 "16b": "Discussion, paragraphs 2 to 4 (tissue of the candidate versus site of drug action; LD at IL6R; APOC3 locus). Causal language is avoided for candidates",
 "16c": "Discussion, last paragraph before Limitations; Conclusions (reporting recommendations). No effect sizes of interventions are inferred",
 "17": "Limitations (European ancestry; resources checked; one threshold)",
 "18": "Sources of Funding (none)",
 "19": "Methods, Data and Code Availability; Supplemental Methods 1 and 10",
 "20": "Disclosures (none)",
}
md = ["| Item | Section | Checklist item (official wording) | Where addressed in this manuscript |",
      "|---|---|---|---|"]
item, missing = "", []
for r in rows:
    c = [x.replace("|", "/").replace(" ", " ").strip() for x in r] + ["", "", ""]
    if c[0] == "Item No." or (not c[0] and not c[1] and not c[2]):
        continue
    if c[0].isdigit():
        item = c[0]
        key = item
        sec, text = c[1], c[2]
    elif re.fullmatch(r"[a-e]\)", c[1]):
        key = item + c[1][0]
        sec, text = c[1], c[2]
    elif c[1].isupper() or c[0].isupper():
        md.append("| | **%s** | | |" % (c[1] or c[0]))
        continue
    else:
        continue
    if key not in W and text:
        missing.append(key)
    md.append("| %s | %s | %s | %s |" % (c[0], sec, text, W.get(key, "See sub-items" if not text else "**NOT MAPPED**")))
md += ["", "Item wording: STROBE-MR checklist, fillable Word version, strobe-mr.org, downloaded "
       "2026-10-11; copyright Equator Network, CC BY 3.0. The \"Page No.\" and \"Relevant text\" "
       "columns of the official form are replaced here by section references and must be "
       "completed with page numbers at submission."]
open(os.path.join(OUT, "strobe_mr_checklist.md"), "w", encoding="utf-8", newline="\n").write(
    "\n".join(md) + "\n")
print("wrote STROBE-MR checklist: %d rows; unmapped keys: %s" % (len(md) - 4, missing))
