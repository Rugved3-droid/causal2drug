# NUMBERS — verified key numbers

> **Superseded 2026-10-10 by `KEY_NUMBERS.md`.** This file is the Phase 1 record for the original 630-target set, with the old "Class A / Class B" names and the unresolved gene records. Sections 5 (controls) and 6 (Class C) are unchanged and still current.

Every number below was recomputed in this pass from files on disk. None is copied from a
document without a recomputation behind it.

**One command reproduces all of it** (read-only; takes a few seconds):

```powershell
cd D:\causal2drug
python manuscript_prep\verify\verify_numbers.py > manuscript_prep\verify\verify_output.txt
```

Tags like **[V04]** name the block of `verify/verify_output.txt` where the number is printed.
The PCSK9 raw-file scan is a separate command (§5). No script in `src/` was executed, because
several overwrite their own outputs.

Naming: the repository calls the two classes **Class A** (absent or uninstrumented) and
**Class B** (modality mismatch). This document uses those names.

---

## 1. Flow: programs → deduplicated pairs → CV targets  [V01, V02]

| step | N | source |
|---|---|---|
| Raw rows in Minikel `pp.tsv` | **29,476** | `data/raw/pp.tsv` |
| Distinct `ti_uid` (target–indication pairs) | **29,460** | dedup by `ti_uid` |
| After dropping the 4 conflicting programs — frozen set | **29,456** | `data/derived/pp_dedup_frozen.tsv` |
| Rows removed in total | 20 | |
| Distinct targets / MeSH indications in the frozen set | **2,517 / 1,014** | |
| CV-PRIMARY (MeSH C14): pairs / targets / indications | 1,364 / 594 / 37 | |
| **CV-EXTENDED (atlas definition): pairs / targets / indications** | **1,543 / 630 / 44** | C14 + 8-token lipid/atherosclerosis rule |
| Targets added by the extension | 36 | 7 indications, all lipid disorders |

Checks that passed: the frozen file's `ti_uid` set equals the recomputed set; its sha256
`5cb7b801…c153c` equals the value in `DECISIONS.md:47`; each of the four dropped programs is
present in the raw file with exactly two rows carrying conflicting phase labels; all 1,014
indications map to a top-level area with zero misses; PCSK9 is in the set under both definitions.

The "~630" in the brief is exactly **630**.

## 2. Class counts  [V03]

| state | N | % of 630 |
|---|---|---|
| ABSENT from eQTLGen | 177 | 28.1 % |
| TESTED-BUT-NULL (no cis-eQTL at p < 5e-8) | 55 | 8.7 % |
| **Class A total (absent + tested-null)** | **232** | **36.8 %** |
| Instrumented in whole blood | 386 | 61.3 % |
| …of which **Class B (assigned)** | **1** (IL6R) | 0.2 % |
| UNMAPPED (no unique GRCh37 gene record) | 12 | 1.9 % |
| Class A + Class B | 233 | 37.0 % |

Class A as a share of the 618 mapped targets: 37.5 %.

The shipped atlas (`deliverables/cv/cardiovascular_atlas.csv`, 630 rows × 25 columns) agrees with
this recomputation on state, class and maximum phase for **630 of 630 rows**.

Class B is assigned by name in `src/cv_atlas.py:152-153`, not by a rule. Its bound, from
`cv/classB_upper_bound.csv` [V11]: **62 of 386 (16.1 %)** instrumented targets carry a
protein-altering GWAS Catalog association for a cardiovascular trait; 196 (50.8 %) for any trait.

Genome-wide context, all 2,517 targets [V03]: ABSENT 652, TESTED-NULL 229, INSTRUMENTED 1,576,
UNMAPPED 60.

## 3. Class by highest phase reached  [V04]

Phase is the maximum Minikel `ccatnum` over a target's CV-EXTENDED pairs (`src/cv_atlas.py:110-113`).
`ccat` combines historical and active status. Table also saved as `verify/class_by_phase.csv`.

| highest CV phase | all | ABSENT | TESTED-NULL | **Class A** | **Class A %** | instrumented (not B) | Class B | unmapped |
|---|---|---|---|---|---|---|---|---|
| Preclinical | 297 | 84 | 28 | 112 | 37.7 % | 182 | 0 | 3 |
| Phase I | 90 | 22 | 7 | 29 | 32.2 % | 57 | 0 | 4 |
| Phase II | 140 | 44 | 8 | 52 | 37.1 % | 87 | 0 | 1 |
| Phase III | 51 | 9 | 9 | 18 | 35.3 % | 32 | 1 | 0 |
| **Launched** | **52** | 18 | 3 | **21** | **40.4 %** | 27 | 0 | 4 |
| Total | 630 | 177 | 55 | 232 | 36.8 % | 385 | 1 | 12 |

The Class A share is flat across phases (32–40 %). A descriptive chi-square of Class A versus
instrumented across the five phases gives χ² = 1.47, df = 4, p = 0.83. That test is mine, was not
pre-registered, and is reported only to say the data show no phase gradient.

Under CV-PRIMARY: 44 launched targets, 17 Class A (38.6 %) [V05]. The existing documents say 20;
see DISCREPANCIES D1.

## 4. Class by CV indication area  [V06]

**The repository has no indication-area grouping** (lipids, heart failure, arrhythmia,
hypertension, thrombosis). The atlas has no such column, and no script defines one. Building one
means choosing which MeSH terms belong to which area, which is a new analytic decision, so it was
not done here. What exists and is verified is the breakdown by MeSH indication, 44 terms, full
table in `verify/class_by_indication_mesh.csv`.

A target with several indications is counted once under each. Rows therefore do not sum to 630.

| MeSH indication | targets | Class A | Class A % | launched pairs | launched, Class A |
|---|---|---|---|---|---|
| Brain Ischemia | 173 | 65 | 37.6 % | 4 | 2 |
| Stroke | 173 | 65 | 37.6 % | 4 | 2 |
| Cardiovascular Diseases | 170 | 57 | 33.5 % | 1 | 0 |
| Heart Failure | 151 | 62 | 41.1 % | 10 | 3 |
| Atherosclerosis | 132 | 47 | 35.6 % | 2 | 2 |
| Hypertension, Pulmonary | 105 | 37 | 35.2 % | 6 | 2 |
| Hyperlipidemias | 76 | 31 | 40.8 % | 4 | 1 |
| Hypertension | 66 | 31 | 47.0 % | 10 | 5 |
| Thrombosis | 58 | 21 | 36.2 % | 6 | 2 |
| Peripheral Vascular Diseases | 56 | 19 | 33.9 % | 2 | 0 |
| Hypercholesterolemia | 47 | 21 | 44.7 % | 3 | 1 |
| Coronary Restenosis | 45 | 17 | 37.8 % | 6 | 1 |
| Hypertriglyceridemia | 31 | 14 | 45.2 % | 3 | 1 |
| Acute Coronary Syndrome | 28 | 10 | 35.7 % | 1 | 0 |
| Hypotension | 16 | 12 | 75.0 % | 4 | 4 |
| Hyperlipoproteinemia Type II | 12 | 6 | 50.0 % | 3 | 3 |

At pair level, 577 of 1,543 CV pairs (37.4 %) have a Class A target, and 932 have an instrumented
target.

Two things this table shows that matter for the paper. Brain Ischemia and Stroke have identical
counts (173 / 65), so Minikel appears to map them as one set. And **no arrhythmia indication
appears among the 44**: atrial fibrillation is present in the project only as a Part C outcome
trait. Both are in DISCREPANCIES (D16).

## 5. Positive controls  [V08]

**GTEx v10 liver (original resource) — all three fail.** No control reaches p < 5e-8 and none is
an eGene at FDR 0.05; the resource has 2,603 eGenes of 20,540 genes tested (12.7 %).
Top nominal cis p from `d1/egene_status.json`: HMGCR 1.09e-3, PCSK9 1.03e-5, IL6R 5.97e-4.

**eQTLGen phase I whole blood (frozen resource):**

| | HMGCR → LDL-C | IL6R → CRP |
|---|---|---|
| eQTLGen hits p < 5e-8 in ±100 kb | 277 | 226 |
| independent instruments (r² < 0.01) | 1 (rs6453133) | 2 (rs10908839, rs58229059) |
| SNP used for MR | rs6453133 | rs2229238 (clump lead absent from outcome GWAS) |
| exposure p at that SNP / F | 1.21e-50 / 224.0 | 4.27e-31 / 134.5 |
| variants after harmonization | 276 | 167 (91 same, 76 flipped) |
| **MR β (SE)** | **+0.3891 (0.01204)** | **+0.3762 (0.02309)** |
| MR p | 4.24e-229 | 1.16e-59 |
| direction / pre-registered expectation | POSITIVE / POSITIVE | POSITIVE / **NEGATIVE** |
| **verdict** | **PASS** | **FAIL (direction)** |

Source `deliverables/d1c/d1c_summary.json`. Both Wald ratios were recomputed independently from
the per-variant merge tables and match to the printed precision. At IL6R, 94 of the 102 variants
with exposure p < 1e-10 give a positive Wald ratio.

The IL6R failure stands as a failure. Decision 11 reclassifies its interpretation (the
pre-registered expectation rested on a shedding variant that total mRNA cannot represent), not
its verdict. The paper has **one passing positive control**.

**PCSK9 — absent.** Direct scan of the full eQTLGen cis-eQTL file:

```bash
gzip -dc "data/raw/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
 | awk -F'\t' 'NR>1{n++; if(index($0,"ENSG00000169174")>0)e++; if($9=="PCSK9")s++;
     if($9~/^PCSK[4567]$/)c[$9]++} END{print n, e+0, s+0; for(k in c) print k, c[k]}'
```

Result, saved in `verify/pcsk9_scan.txt`: **127,341,798** data rows; **0** contain
`ENSG00000169174`; **0** have gene symbol PCSK9. PCSK4, PCSK5, PCSK6 and PCSK7 have 9,510, 8,583,
6,525 and 8,377 rows. GTEx median TPM for PCSK9: whole blood 0.29, liver 28.3 (atlas row).

## 6. Class C test  [V09]

**Criterion** (`PROTOCOL.md:62-72`, written 2026-08-30 before the computation): if PP4 increases
with exposure power in the paired head-to-head, the Class C power-dependence claim is refuted,
reported as a refutation, and Class C is dropped.

**Design:** coloc.abf PP4 for the same locus, outcome GWAS, window and harmonization, with the
exposure taken from GTEx v10 whole blood (n = 853) and from eQTLGen (n = 31,684). 60 loci
attempted; **50 paired**; 10 lost because GTEx had too few overlapping variants.

| quantity | value |
|---|---|
| PP4 higher in eQTLGen / lower / ties | **32 / 18 / 0** |
| median PP4, GTEx blood → eQTLGen | 0.00821 → **0.01390** |
| median log10 difference | +0.031 |
| **Wilcoxon signed-rank on log10 PP4** | **W = 408.0, p = 0.0262** |
| exact binomial, 32 of 50 | p = 0.0649 |
| median lead \|z\|, GTEx → eQTLGen | 3.6 → 24.4 |
| Spearman of the difference vs eQTLGen \|z\| | ρ = −0.455, p = 0.0009 |

**Verdict: REFUTED. Class C is dropped.** PP4 rises with exposure power, the opposite of the
prediction.

Three further facts from the same file that the existing write-ups do not state:

* **No locus reaches PP4 > 0.5 in either arm.** The largest PP4 is 0.091 (eQTLGen) and 0.142
  (GTEx). The test compares values that are all in the "no colocalization" range.
* 31 of the 50 paired loci are genes in the 630-target CV atlas; the test is not CV-specific.
* The paired loci span 19 distinct outcome GWAS.

## 7. Supporting numbers that appear in the existing drafts

| number | value | block |
|---|---|---|
| Instrument yield, CV targets: blood / heart LV / atrial / aorta / liver / 4-tissue union | 61.3 / 7.1 / 6.8 / 8.7 / 4.9 / 17.5 % | [V07] |
| Class A enrichment for tissue-restricted expression (ABSENT vs INSTRUMENTED) | OR 10.02, 95 % CI 4.77–21.03, Fisher p 4.05e-15; 2×2 = 156, 8, 253, 130 | [V07] |
| Required modality for the 232 Class A targets | cardiac 101, tissue 50, liver 35, vascular 19, pQTL 15, unknown 12 | [V07] |
| Part C pairs lost to unindexed GWAS files, all areas | 578 of 1,310 (44.1 %), 70 accessions | [V10] |
| Same, CV pairs under Part C rules | 57 of 103 (55.3 %), 4 accessions | [V10] |
| Same, all CV pairs mappable to a EUR harmonised GWAS | 265 of 448 (59.2 %); CV-PRIMARY 218 of 401 (54.4 %) | [V10] |
| R15 loci screened / outcome genome-wide significant | 106 / 6 | [V10] |
| …restricted to CV-atlas genes | 65 screened, 3 significant, 2 with a credible-set result | [V10] — see D4 before using |

## 8. Not verified

| item | why |
|---|---|
| Every upstream step that needs the network (clumping against 1000G, GWAS region queries, `.tbi` probes) | Recomputed from saved outputs only; not re-run |
| Class B bound from the raw GWAS Catalog zip | Verified from `classB_upper_bound.csv`; the zip was not re-parsed |
| Tissue-yield counts from the raw eQTL Catalogue permuted files | Verified from the atlas columns; the permuted files were not re-parsed |
| Whether each "Launched" label corresponds to a currently approved drug for that indication | The repository has no drug-to-indication link. See HEADLINE and D2 |
| DENTIST RMSE figures, SuSiE results, UKB LD result | Outside the six requested number groups; not recomputed |
