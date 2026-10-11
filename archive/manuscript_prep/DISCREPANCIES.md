# DISCREPANCIES — inconsistent, unverified, missing, stale

> **Status 2026-10-10.** Resolved by decision: D1 ("20" corrected to 17 of 44), D2/D3 (ChEMBL link built), D4 (claim removed), D12 (two named sub-types), D13 (gene records resolved; ACE instrumented), D16 (CV set extended), D17 (Class B retired), D20 (deviation note). Reported, not changed: D18 (label mismatches → Supplementary Table S6), D19 (no phase gradient), D21 (one passing control). Not acted on: D5, D7, D8, D9, D10, D14, D15, D23, D24 are errors in existing documents that were outside the instruction to edit; D11 (the modality heuristic) is dropped from the paper. **New in the rerun:** eQTLGen tested no X-linked gene (`KEY_NUMBERS.md` §6 and §10).

Nothing here was resolved. Each item states what conflicts, where, and what I measured. Items are
ordered by how much they affect the paper. Blocks **[Vnn]** refer to
`verify/verify_output.txt`; reproduce with `python manuscript_prep\verify\verify_numbers.py`.

**What did reproduce:** the flow (29,476 → 29,460 → 29,456), the 630-target set, all four state
counts, the 232 / 36.8 % Class A figure, the 21 launched and 18 Phase III Class A counts, the
enrichment odds ratio, the tissue yields, both control results, the PCSK9 absence, and every
Class C statistic. The atlas CSV matches a from-scratch recomputation on all 630 rows.

---

## A. Affects the headline

**D1 — CV-PRIMARY launched Class A count: documents say 20, strict recomputation gives 17.**
`ANALYSES_1_2_3.md:88`, `CARDIOVASCULAR_ATLAS.md:170`, `LIMITATIONS.md:274`,
`cv_robustness.json` (`classA_launched: 20`). `src/a2_robustness.py:84` restricts to CV-PRIMARY
targets but reads `max_phase_num`, which was computed over CV-EXTENDED pairs. ANGPTL3, APOB and
NPC1L1 are in CV-PRIMARY through a lower-phase C14 indication and are launched only in lipid
indications, so they are counted as launched under a definition that excludes the indications
where they are launched. With phase taken over CV-PRIMARY pairs: **44 launched targets, 17
Class A (38.6 %)** [V05]. The share is stable; the stated count is not.

**D2 — No drug-to-indication link exists anywhere in the repository.**
Minikel `pp.tsv` has no drug names. `src/build_table1.py:38-46` takes the example agent from Open
Targets by target alone. The agent in the existing Table 1 is therefore not the drug that earned
the phase label in the listed indication. Examples from `table1_classA_launched_phase3.csv`:
ADRA1A/B/D "Prazosin, antagonist" against Hypotension; PLG "Aprotinin" against Thrombosis; AGT
"GSK-2586881" against a Launched label; ABCC8 "Glipizide" against Brain Ischemia; SLC6A2 "Doxepin"
against Brain Ischemia. Five rows read "not in Open Targets" (ADRA2B, APOA1, CNTFR, HCN1, HCN2)
[V05]. The brief asks for "the drug(s), indication" per target; that pairing is **UNVERIFIED** and
cannot be produced from current data. The Open Targets 26.06 release manifest also reports
`failure` (`DECISIONS.md:84-93`, deferred).

**D3 — Table 1 shows one indication per target, chosen by file order.**
`src/build_table1.py:62` keeps the first pair encountered at the maximum phase. 13 of 39 rows have
more than one CV indication at that phase [V05]. F2 is shown as "Brain Ischemia" and has six;
EDNRB is shown as "Brain Ischemia" and is also launched for pulmonary hypertension.

**D16 — The cardiovascular set excludes myocardial infarction, angina and every arrhythmia.**
Minikel's top-level map files these under "signs/symptoms", not C14. The atlas documents the
lipid exclusion and fixes it with a token rule (`CARDIOVASCULAR_ATLAS.md:24-33`) but says nothing
about these. Measured with `python manuscript_prep\verify\check_excluded_cardiac_indications.py`:
Myocardial Infarction 129 pairs (7 launched), Atrial Fibrillation 33 (5), Angina Pectoris 33 (2),
Arrhythmias, Cardiac 16, Angina, Unstable 13 (7), Ventricular Fibrillation 6 (4), plus dilated
cardiomyopathy (filed under "congenital"). Across 15 such indications: 259 pairs, 195 targets,
**30 targets not in the atlas at all**. KCNH2 carries a Launched label there and is absent from
the atlas; SCN5A and AVPR2 are in the atlas at a lower phase than they hold in these indications.
The set also *includes* Retinal Vein Occlusion, Retinal Artery Occlusion, Retinal Telangiectasis
and Ischemic Optic Neuropathy. Also: Brain Ischemia and Stroke have identical counts (173 targets,
65 Class A) [V06]. The CV set rule is not on the frozen list in the brief, so this is yours to
decide; I changed nothing.

**D13 — Four launched targets have no class because of gene-record matching, and one is ACE.**
`gene_coords_unmapped.tsv`: ACE and ITGB3 are "ambiguous_2_gene_records" in the GRCh37 release 87
GTF; GUCY1A1 and GUCY1B1 are "not_in_GRCh37_87_GTF" (newer symbols). They are excluded from class
statistics, yet every percentage uses 630 or 52 as the denominator. With mapped denominators,
Class A is 37.5 % of 618 and 43.8 % of 48 launched [V03, V05]. ADRA2B is classed ABSENT but has
no GTEx or HPA values, so its expression is unknown. Six "targets" in the absent group are
non-coding or lack expression data (MIR122, MIR143, MIR145, MIR21, GDF2, ADRA2B).

**D12 — "Untestable" means two different things.**
Class A is ABSENT + TESTED-NULL (232). The enrichment result and Figure 4 compare ABSENT alone
against INSTRUMENTED (`src/cv_atlas.py:205-208`), while Figure 4's axis says "untestable target".
Genome-wide the tested-null group is only weakly tissue-restricted (OR 1.46) and not enriched
for secreted proteins (OR 0.82, p = 0.27) (`classA_tests.json`). Three launched Class A targets
are tested-null, and two of them are highly expressed in blood: ITGA2B 36.3 TPM, IL1B 29.5 TPM.
Figure 1 labels the tested-null box "no usable genetic marker" but the summary text
(`supplementary/README.md:12`) says 232 "cannot be instrumented at all in whole blood".

**D17 — The brief's description of the two classes does not match the repository.**
The brief describes class 2 as "modality mismatch (e.g., tissue-restricted secreted proteins)".
In the repository, tissue-restricted and secreted proteins characterise **Class A**. Class B is
"well instrumented but the exposure cannot represent the drug mechanism" and has **one member,
IL6R, assigned by name** in `src/cv_atlas.py:152-153`. Its size is bounded (62 of 386 on CV
traits), not measured. A paper framed as "two classes" has one class of 232 and one of 1.

**D18 — Phase labels could not be checked, and some look wrong in both directions.**
From my general knowledge, not from the repository, and therefore **UNVERIFIED**: I cannot match
the Launched label to an approved CV drug for AGT, IL1B (atherosclerosis) or RXFP1 (serelaxin);
and HCN4 is labelled Phase III although ivabradine is an approved heart-failure drug. Dependence
on Minikel labels is a frozen design choice; this is recorded for the limitations section.

**D19 — Class A is not concentrated in launched targets.**
Not a conflict between files, but the existing text reads as if approved drugs were singled out.
Class A share by phase: 37.7 / 32.2 / 37.1 / 35.3 / 40.4 % [V04]; all CV 36.8 %; all 2,517
targets 35.0 %.

## B. Affects the Class C result and the colocalization statements

**D20 — The Class C refutation is built on the method the project declared untrustworthy.**
`PROTOCOL.md:135` and `DECISIONS.md:259-262`: "No colocalization output from this project may be
used, cited, tabulated or plotted." The Class C test is a comparison of coloc.abf PP4 values, the
existing Table 2 tabulates colocalization counts, and Phase 2 item 7(d) asks for a Class C
figure. Either the standing rule has an unstated exception for these or it conflicts with the
plan. Related, from [V09]: **no paired locus has PP4 above 0.5 in either arm** (maxima 0.091 and
0.142), so the refutation compares 50 values that all say "no colocalization". The documents
argue the pairing makes this harmless (`PART_3_DEGRADATION.md:103-105`); a reviewer may not agree.
Only 31 of the 50 loci are CV-atlas genes.

**D4 — "Colocalization is unavailable for 96.9 % of cardiovascular target loci" is mostly not
about cardiovascular outcomes.** `CARDIOVASCULAR_ATLAS.md:154-158`, `table2_method_availability.csv`,
`cv_report.json`. The 65 loci are CV-atlas *genes* (`src/cv_atlas_report.py:96`), each paired with
whatever outcome Part C selected: Asthma 17, Heart Failure 12, Idiopathic Pulmonary Fibrosis 8,
Alzheimer Disease 6, Sepsis 4, and others. **Only 12 of the 65 have a CV-EXTENDED indication as
the outcome** [V10]. The three with a genome-wide significant outcome signal are ADORA1 / Atrial
Fibrillation, PPARG / Asthma and TNF / Asthma, none of which is a CV-EXTENDED indication.

**D7 — A withdrawn figure is still asserted in two documents.**
"Credible-set colocalization returned an answer at only 2 of 28 loci (7 %)" appears as a live
finding in `INSTRUMENT_AVAILABILITY.md:213-215` and `PART_3_DEGRADATION.md:75-77, 130`. It was
withdrawn in `FIXES_1_2_3.md:85-87` and `LIMITATIONS.md:297-299`, and Table 2 prints it as
"WITHDRAWN". The DENTIST figures quoted beside it in `INSTRUMENT_AVAILABILITY.md:207-211` (median
4.46, 93 % above 2, 28 loci) come from the same flawed rule; the corrected values are 3.15 and
100 % on 6 loci.

**D8 — Two documents still say the coloc failure "worsens as the exposure resource grows".**
`PROTOCOL.md:33` and `LIMITATIONS.md:108`. That is the Class C power-dependence claim, which
`PROTOCOL.md:140-150` and `LIMITATIONS.md:185` record as refuted. `LIMITATIONS.md` contains both
statements.

**D15 — Counts of withdrawn artefacts disagree.** `PROTOCOL.md:120` says three artefactual PP4
values; `LIMITATIONS.md:333` and `FIXES_1_2_3.md:20` say four. The HMGCR / 1000G value is 1.000
in `PROTOCOL.md:127` and 0.998 in `LIMITATIONS.md:341`.

## C. Errors and gaps in the control and methods record

**D9 — IL6R "lead p = 1.16 × 10⁻⁵⁹" is the MR p-value, not an instrument p-value.**
`CARDIOVASCULAR_ATLAS.md:87-88`, `DECISIONS.md:236`, `INSTRUMENT_AVAILABILITY.md:117`. From
`d1c_summary.json` [V08]: the clumped lead rs10908839 has exposure p = 4.09e-34 (F = 148.3); it
is absent from the CRP GWAS, so MR used rs2229238 with exposure p = 4.27e-31 (**F = 134.5**);
1.16e-59 is the p of the MR estimate. The F is right; the p beside it is from a different
quantity.

**D21 — The paper has one passing positive control and a null validation.**
HMGCR passes; IL6R fails on direction and is reclassified, not rescued (`DECISIONS.md:206-224`);
PCSK9 is absent; named controls are retired (`DECISIONS.md:245-253`). The replacement validation,
Part C, returned concordance of 287 / 571 = 50.3 %, binomial p = 0.933
(`PART_C_CONCORDANCE.md:11-12`; not recomputed here). A reviewer will ask what shows the MR
pipeline works. This does not touch the atlas classification, which needs only instrument
lookup, but any sentence about MR results rests on a single control.

**D14 — PCSK9's best liver cis-eQTL has two values.** `DELIVERABLE_1.md:54` and
`DECISIONS.md:111`: p = 4.57e-6, F = 22.26, rs1165287 (`chr1_55054539_G_A`).
`d1/egene_status.json`: top nominal p = 1.03e-5 at `chr1_55055322_C_CG`. Different variants,
possibly different files or windows. Both are above 5e-8, so the control verdict is unaffected.
Unresolved.

**D10 — The launched table in the atlas write-up lists 19 of 21.** `CARDIOVASCULAR_ATLAS.md:94-112`
omits AGT and SLC5A1; the CSV has both.

**D11 — The "required modality" rule gives clinically odd answers.** The rule
(`src/cv_atlas.py:124-125`) assigns "cardiac eQTL" whenever heart TPM ≥ 1, before checking liver.
APOB (heart 1.8, liver 350), AGT (24.1 / 312.6), AGTR1 (1.4 / 28.2) and APOA1 (17.6 / 8,481) are
all labelled cardiac [V07]. This feeds "43.5 % need cardiac eQTL" and Figure 3. The documents call
the rule a heuristic; they do not show these cases.

**D5 — `INSTRUMENT_AVAILABILITY.md` is partly stale.** Lines 137-147 say Class B "remains
unquantified" and "nothing measured so far bounds how large Class B is"; Analysis 1 bounded it the
next day (`ANALYSES_1_2_3.md:29-41`). Line 79 heads the tissue section "not yet run" while its
yields are reported as results elsewhere. See also D7.

**D22 — The infrastructure-attrition denominator includes pairs whose GWAS file was never
probed.** `src/a3_cascade.py:132` counts a pair as lost when its accession is in the skip list.
Of the 448 CV pairs, 265 are in the skip list, 181 map to an accession that produced results, and
**2 map to an accession Part C never touched** [V10]. Small, but 59.2 % is a floor, not an exact
rate. `ANALYSES_1_2_3.md:129` also describes a skip label (`missing_effect_or_se_columns`) that
the file no longer contains; it now reads `plain_gzip_no_tbi_not_range_queryable` for all 70.

**D23 — `cv_report.json` still carries 57 / 103 with no "sub-rate" label**; the corrected 265 /
448 lives only in `cv_robustness.json`. Anyone reading the JSON alone gets the superseded figure.

**D24 — `supplementary/README.md:54` says "Exactly two commands" and lists three.**

## D. Missing

| item | status |
|---|---|
| CV indication-area grouping (lipids, HF, arrhythmia, hypertension, thrombosis) | **Does not exist.** Only the 44 MeSH terms. Needs a mapping rule from you before Table 3 can be built |
| Drug names per target–indication pair | **Does not exist** (D2) |
| Figures for flow, class by phase, launched/Phase III highlight, Class C | Do not exist; existing Figures 1–4 cover other content |
| Retrieval URL and timestamp for GTEx v10 TPM, Human Protein Atlas, GWAS Catalog associations, Open Targets 26.06 | Not in `logs/retrieval_log.csv` (self-reported in `LIMITATIONS.md:429`) |
| Licence terms of incorporated sources | Self-marked UNVERIFIED (`supplementary/README.md:46-50`) |
| Genome-wide liftover failure rate | Measured on 503 variants in two loci only |
| Version control | `D:\causal2drug` is not a git repository; there is no commit history to cite for "frozen before results" claims, only file timestamps and the documents' own dates |
| Software and package versions | No `requirements.txt` or environment file. This machine has Python 3.12.10, pandas 3.0.1, scipy 1.17.1, numpy 2.4.2; whether those produced the original outputs is **UNVERIFIED** |
| Citation and version of the Minikel dataset | `pp.tsv` git blob SHA-1 matches GitHub per `DELIVERABLE_0.md:179-196`; I confirmed the blob hash locally (`git hash-object data\raw\pp.tsv` → `16d64118…`), not the GitHub side |

## E. Not examined

* `C:\Users\Deven\causal2drug_feasibility\` (predecessor stage) was listed, not read.
* `PREREGISTRATION.md` was read by section heading only; it is the main source for Phase 2
  methods and I will read it in full then.
* Numbers outside the six requested groups (DENTIST, SuSiE, UKB LD, Part B genome-wide tests)
  were not recomputed.
