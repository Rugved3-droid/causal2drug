# RESCUE_CONTROLS_RESULTS — positive controls with candidate instruments from the recovery map

Run 2026-10-11 under `RESCUE_CONTROLS_PROTOCOL.md`, which was written first.

**Corrected later on 2026-10-11 after the second-implementation audit.** The protocol takes the lead variant and sheet recorded in `recovery_map.csv`. A P-value underflow in the map builder had recorded the UKB-PPP discovery sheet for PCSK9 and a non-lead variant for ANGPTL3. With the map corrected, PCSK9 uses the combined-cohort row (same variant; β changes from +0.376 to +0.387, F from 1,723 to 2,305) and ANGPTL3's lead is rs775677524 (still an indel). No verdict changed. The protocol was not edited. First-run numbers: β −1.153 (0.028), F 1,723, MR +0.376 (+0.367 to +0.385) for PCSK9; rs763015008 for ANGPTL3.

**Result: 2 of 5 passed, 1 failed, 2 could not be harmonized.** No variant, study or tissue was
substituted.

Reproduce: `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py`
Every number below is in `rescue_controls\run_log.txt` and
`rescue_controls\rescue_controls_results.csv` (one row per control), or in
`display\tableS19_candidate_instrument_controls.md`.

| control | candidate instrument | exposure β (SE), F | outcome study | MR β (95% CI), p | expected / observed | verdict |
|---|---|---|---|---|---|---|
| PCSK9 → LDL cholesterol | UKB-PPP plasma cis-pQTL, rs11591147, allele T, 0 bp from gene | −1.121 (0.023), F 2,305 | GCST90239658, N 1,320,016 | **+0.387 (+0.377 to +0.396)**, p < 1 × 10⁻³⁰⁰ | positive / positive | **PASS** |
| ANGPTL3 → triglycerides | UKB-PPP plasma cis-pQTL, rs775677524, a 9-bp deletion (AGTTAATGTG>A), 0 bp from gene | −0.282 (0.006), F 2,296 | GCST90239664, N 1,320,016 | none | positive / — | **NOT HARMONIZABLE**: the variant is an indel, and rule R6 drops indels |
| APOC3 → triglycerides | Fenland plasma cis-pQTL, rs964184, allele C, 51,505 bp from gene | −0.210 (0.020), F 115 | GCST90239664, N 1,320,016 | **+1.171 (+1.153 to +1.190)**, p < 1 × 10⁻³⁰⁰ | positive / positive | **PASS** |
| AGTR1 → systolic blood pressure | GTEx v8 tibial artery cis-eQTL, chr3:148624775 T>G (rs275692), allele G | +0.212 (0.027), F 61 | GCST90310294, N 1,028,980 | **−0.037 (−0.268 to +0.194)**, p = 0.75 | positive / negative | **FAIL** (not significant; point estimate has the wrong sign) |
| REN → systolic blood pressure | UKB-PPP plasma cis-pQTL, rs193280350, allele A, frequency 0.0076, 0 bp from gene | +0.835 (0.035), F 578 | GCST90310294, N 1,028,980 | none | positive / — | **NOT HARMONIZABLE**: no row at the variant's position in the outcome file |

Units: PCSK9, SD of LDL cholesterol per SD higher plasma PCSK9; APOC3, SD of log-transformed
triglycerides per SD higher plasma apolipoprotein C-III; AGTR1, mm Hg of systolic blood
pressure per unit of GTEx normalized expression. Sources of the scales and the sample-overlap
assessment for every control are in `display\tableS22_control_scales_and_overlap.md`.
The Wald-ratio standard error does not carry the exposure's uncertainty (pipeline rule R12.2).

## How each instrument and outcome was fixed

- Instruments are the lead qualifying variants already recorded in `recovery_map.csv`. PCSK9,
  ANGPTL3 and REN have a UKB-PPP instrument; APOC3 is not assayed in UKB-PPP, so Fenland was
  next under the protocol's order; AGTR1 is assayed by no proteomic study, so the
  cardiovascular-relevant tissue with the lowest lead p (tibial artery) was used.
- Outcome studies were chosen by the pipeline's own function (`src/select_gwas.py`, R2 + R2a):
  triglycerides EFO_0004530 → GCST90239664 (Graham; 1,007 European-only studies with a tabix
  index were eligible); systolic blood pressure EFO_0006335 → GCST90310294 (Keaton; 26
  eligible). Selection tables: `rescue_controls\gwas_selection_<EFO>.csv`.
- Harmonization used `src/harmonize.py` unchanged. All three estimates resolved
  `aligned_same`; allele frequencies agree between exposure and outcome (PCSK9 T: 0.016 and
  0.016; APOC3 C: 0.867 and 0.865).

## What to make of it

1. **The two passes are both lipid targets with a plasma pQTL.** PCSK9 is the control that
   could not be run at all in eQTLGen. It passes with its coding variant rs11591147.
2. **APOC3 passes on a variant that is not specific to APOC3.** rs964184 lies 51.5 kb from the
   gene, in the *APOA1/C3/A4/A5* cluster (the Fenland workbook, table ST5, lists it against
   APOA5 and its neighbours), and its triglyceride association is extremely strong
   (outcome p = 2.6 × 10⁻³³³⁶). The estimate (+1.17 per unit of protein) is large. The
   pass shows a direction, not that the instrument isolates APOC3. This is the "target
   specificity" caveat in a concrete case.
3. **The one tissue-eQTL control fails.** The tibial-artery AGTR1 variant has no detectable
   association with systolic blood pressure in over a million people (outcome β −0.008,
   SE 0.025, p = 0.75). An eQTL at p = 3 × 10⁻¹⁴ in a cardiovascular-relevant tissue did not
   reproduce the effect of a first-line drug class.
4. **Two candidates could not be used under the pipeline's own rules.** The ANGPTL3 lead is an
   indel and the REN lead is a rare variant absent from the outcome GWAS. Being listed in the
   recovery map did not make either usable for two-sample MR as specified.

Not an estimate under the protocol, reported only so it is not hidden: the outcome file does
contain the ANGPTL3 deletion (effect allele A, β −0.0709, SE 0.0016). With the exposure β of
−0.282 for the same allele, the ratio would be positive. It is **not counted as a pass**,
because the harmonization rule excludes indels and the protocol forbids changing the rule.

## Scales and sample overlap (added in the final round)

| control | exposure scale | outcome scale | overlap of exposure and outcome samples |
|---|---|---|---|
| PCSK9 | UKB-PPP combined cohort (52,363 UK Biobank participants, some of non-European ancestry), inverse-rank-normalized protein (SD) | Graham 2021, inverse-normal residuals (SD) | Yes: UK Biobank is in both. Up to all 52,363 exposure participants, at most 4.0% of the outcome sample; exact number not verified |
| ANGPTL3 | as PCSK9 | Graham 2021, log then inverse-normal (SD) | as PCSK9; no estimate |
| APOC3 | Fenland (10,708), per SD of normalized abundance | as ANGPTL3 | Not established; Fenland is not named in the outcome paper's main text. At most 0.8% of the outcome sample |
| AGTR1 | GTEx v8 tibial artery (584), regression slope on normalized expression | Keaton 2024, mm Hg | None known |
| REN | UKB-PPP combined cohort (52,363) | Keaton 2024, mm Hg | Yes had there been an estimate: at most 5.1% of the outcome sample |

With overlapping samples, weak-instrument bias is toward the confounded observational
association and shrinks with instrument strength (Burgess, Davies and Thompson, *Genet
Epidemiol* 2016, doi 10.1002/gepi.21998). F is 2,305 for PCSK9 and 115 for APOC3, so that
bias should be small. Extracts supporting each scale and cohort statement:
`revision_sources\units_and_overlap_extracts.txt`. The GTEx file documentation calls `slope`
only a regression slope; that it refers to the ALT allele is the GTEx convention and was not
confirmed from that documentation. The AGTR1 verdict does not depend on it.

## How to read the pattern

Reported as **2 pass / 1 fail / 2 not harmonizable.** Both controls that passed used a
plasma pQTL for a circulating protein that the drug lowers; the one tissue-eQTL control
failed. With five controls this is illustrative, not a general result about pQTLs and
eQTLs. The two that could not be harmonized are a finding in their own right: a candidate
that is available is not necessarily usable. No rule was changed for ANGPTL3 or REN.

## Consequence for the manuscript

The recovery map counts *candidate* instruments. In this test of five, two reproduced a known
drug direction, one did not, and two could not be analysed. With the main controls (HMGCR
pass, IL6R fail on the pre-specified sign, PCSK9 untestable in eQTLGen), validation of the
MR pipeline remains thin, and no claim in the paper rests on an MR estimate.
