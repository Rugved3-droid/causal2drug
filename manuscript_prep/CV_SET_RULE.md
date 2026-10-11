# CV_SET_RULE — the cardiovascular set

**Status: APPROVED 2026-10-10 (option 2 is primary; no post-hoc exclusions; Reperfusion Injury
stays in). Implemented in `analysis/classify.py`.** The earlier proposal text is superseded by
this file.

Reproduce the term lists and counts:

```powershell
cd D:\causal2drug
python manuscript_prep\verify\list_cv_candidates.py     # candidate terms, pair counts
python manuscript_prep\analysis\classify.py             # block [K01]; writes analysis\out\cv_terms_primary.csv
```

## The rule

An indication is cardiovascular if **either** holds:

1. Minikel's own `areas` field (`data/raw/minikel_indic.tsv`) contains `cardiovascular`; or
2. its MeSH term contains one of eight tokens: `hyperlipid`, `hypercholesterol`, `dyslipid`,
   `hyperlipoprotein`, `hypertriglycerid`, `lipoprotein`, `atheroscler`, `arterioscler`.

No term is added or removed by hand. The second clause exists because Minikel tags the lipid
disorders `metabolic` only.

The original atlas read a different Minikel file (`minikel_indic_topl_match.tsv`), which assigns
each indication a single branch. That is why myocardial infarction, angina and the arrhythmias
were missing: Minikel tags them cardiovascular but the single-branch file filed them under
signs/symptoms. The script asserts that every C14-branch term also carries the cardiovascular
tag, so the new rule contains the old set exactly.

## The sets

| set | definition | terms | pairs | targets |
|---|---|---|---|---|
| **Primary** | the rule above | **93** | **2,128** | **726** |
| Sensitivity 1 | original atlas: single branch = C14, or lipid token | 44 | 1,543 | 630 |
| Sensitivity 2 | sensitivity 1 + tagged cardiovascular and filed under C23 (signs/symptoms) | 62 | 1,879 | 663 |
| Sensitivity 3 | single branch = C14 only | 37 | 1,364 | 594 |

How each of the 93 terms gets in: 37 by the C14 branch, 7 by a lipid token, 18 tagged
cardiovascular and filed under signs/symptoms, 31 tagged cardiovascular and filed elsewhere.
The full list with rule, area and counts is `tables/tableS2_cv_indication_terms.csv`.

## The 18 terms from signs/symptoms

Myocardial Infarction (129 pairs), Reperfusion Injury (49), Angina Pectoris (33), Atrial
Fibrillation (33), Arrhythmias, Cardiac (16), Angina, Unstable (13), Cerebral Hemorrhage (13),
Intermittent Claudication (11), Subarachnoid Hemorrhage (10), Ventricular Fibrillation (6),
Cerebral Infarction (6), Edema, Cardiac (5), Tachycardia, Ventricular (4), Tachycardia,
Supraventricular (3), Long QT Syndrome (2), Tachycardia (1), Bradycardia (1), CADASIL (1).

## The 31 additional terms filed outside C14 and signs/symptoms

These are in the primary set and not in sensitivity 2. Together: 249 pairs. **Flag** marks terms
that are plainly not primarily cardiac. **None is removed.**

| MeSH ID | term | Minikel areas | branch | pairs | launched pairs | Table 3 area | flag |
|---|---|---|---|---|---|---|---|
| D003930 | Diabetic Retinopathy | ophthalmology, cardiovascular, endocrine | C19 | 99 | 2 | Other | **not primarily cardiac** |
| D017719 | Diabetic Foot | cardiovascular, dermatology, endocrine | C19 | 37 | 1 | Other | **not primarily cardiac** |
| D001528 | Behcet Syndrome | ophthalmology, cardiovascular, congenital, dermatology | C17 | 16 | 5 | Other | |
| D013683 | Telangiectasia, Hereditary Hemorrhagic | cardiovascular, hematology, congenital | C16 | 11 | 0 | Other | |
| D054179 | Angioedemas, Hereditary | cardiovascular, congenital, dermatology, immune | C20 | 10 | 2 | Other | **not primarily cardiac** |
| D002311 | Cardiomyopathy, Dilated | cardiovascular, congenital | C16 | 9 | 0 | Heart failure / cardiomyopathy | |
| D015140 | Dementia, Vascular | neurology, cardiovascular, psychiatry | F03 | 6 | 0 | Cerebrovascular | |
| D014890 | Granulomatosis with Polyangiitis | respiratory, cardiovascular, immune | C20 | 5 | 2 | Other | |
| D004535 | Ehlers-Danlos Syndrome | cardiovascular, hematology, congenital, dermatology | C17 | 5 | 0 | Other | **not primarily cardiac** |
| D014647 | Varicose Ulcer | cardiovascular, dermatology | C17 | 5 | 0 | Other | |
| D013625 | Takayasu Arteritis | cardiovascular, dermatology | C17 | 4 | 0 | Other | |
| D015267 | Churg-Strauss Syndrome | cardiovascular, hematology, immune | C20 | 4 | 0 | Other | |
| D058065 | Diabetic Cardiomyopathies | cardiovascular, endocrine | C19 | 3 | 0 | Heart failure / cardiomyopathy | |
| D016657 | Cerebral Amyloid Angiopathy | neurology, cardiovascular, metabolic | C18 | 3 | 0 | Cerebrovascular | |
| D017241 | MELAS Syndrome | musculoskeletal, neurology, cardiovascular, congenital, metabolic | C18 | 3 | 0 | Other | **not primarily cardiac** |
| D013771 | Tetralogy of Fallot | cardiovascular, congenital | C16 | 3 | 0 | Other | |
| D000795 | Fabry Disease | neurology, cardiovascular, congenital, metabolic | C18 | 3 | 1 | Other | **not primarily cardiac** |
| D014424 | Turner Syndrome | cardiovascular, congenital, endocrine | C19 | 3 | 1 | Other | **not primarily cardiac** |
| D009080 | Mucocutaneous Lymph Node Syndrome | cardiovascular, hematology, dermatology | C17 | 3 | 1 | Other | |
| D009634 | Noonan Syndrome | musculoskeletal, cardiovascular, congenital, dermatology | C17 | 3 | 0 | Other | **not primarily cardiac** |
| D055953 | Microscopic Polyangiitis | neurology, cardiovascular, immune | C20 | 3 | 2 | Other | |
| D004374 | Ductus Arteriosus, Patent | cardiovascular, congenital | C16 | 2 | 0 | Other | |
| D001260 | Ataxia Telangiectasia | neurology, cardiovascular, congenital, metabolic, immune | C20 | 1 | 0 | Other | **not primarily cardiac** |
| D058088 | Sarcoglycanopathies | musculoskeletal, respiratory, neurology, cardiovascular, congenital | C16 | 1 | 0 | Other | **not primarily cardiac** |
| D006623 | von Hippel-Lindau Disease | neurology, cardiovascular, congenital | C16 | 1 | 1 | Other | **not primarily cardiac** |
| D000799 | Angioedema | cardiovascular, dermatology, immune | C20 | 1 | 0 | Other | **not primarily cardiac** |
| D008382 | Marfan Syndrome | musculoskeletal, cardiovascular, congenital, dermatology | C17 | 1 | 0 | Other | **not primarily cardiac** |
| D016738 | Alagille Syndrome | digestive, cardiovascular, congenital | C16 | 1 | 0 | Other | **not primarily cardiac** |
| D001165 | Arteriovenous Malformations | cardiovascular, congenital | C16 | 1 | 0 | Other | |
| D052120 | Glycogen Storage Disease Type IIb | neurology, cardiovascular, congenital, metabolic | C18 | 1 | 0 | Other | **not primarily cardiac** |
| D006330 | Heart Defects, Congenital | cardiovascular, congenital | C16 | 1 | 1 | Other | |

**15 terms are flagged: 170 pairs.** Diabetic Retinopathy and Diabetic Foot are 136 of those.
46 targets are in the primary set only through a flagged term.

The flag is my judgement: a multisystem syndrome or a diabetic complication where the
cardiovascular tag reflects one vascular feature. The systemic vasculitides (Behçet, the
polyangiitides, Takayasu, Kawasaki) are not flagged, because "Vasculitis" was already in the
original set. Marfan syndrome is flagged although its main risk is aortic; that one could
reasonably go either way.

## What the flagged terms do to the results

They do not move the proportions. Sensitivity 2 contains none of the 31 terms:

| | primary (93 terms) | sensitivity 2 (62 terms) |
|---|---|---|
| No eQTL data, lower bound | 217 of 726 (29.9 %) | 191 of 663 (28.8 %) |
| Upper bound | 280 (38.6 %) | 251 (37.9 %) |
| Launched, ChEMBL-verified, no eQTL data | 20 of 58 (34.5 %) | 17 of 49 (34.7 %) |

They do change which names appear. Three launched targets in the primary headline list are
launched only through one of the 31 terms: **GHR** (Turner syndrome), **GLA** (Fabry disease)
and **PDE4C** (Behçet syndrome). A cardiology reader will not recognise somatropin or migalastat
as cardiovascular drugs. I recommend the paper quote the primary proportion and name examples
from the 17 that survive in sensitivity 2.

## Not changed

Dedup, the four dropped programmes, the lipid tokens, and the fact that Brain Ischemia and
Stroke are two indications with identical target lists.
