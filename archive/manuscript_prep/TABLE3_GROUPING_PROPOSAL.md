# TABLE 3 — proposed indication-area grouping

> **Decided 2026-10-10.** Nine areas plus an "Unspecified CV" row; Coronary Thrombosis in ischaemic heart disease; other calls by the disease the term names. The final placement of all 93 terms is logged in `DECISIONS_LOG.md` and `tables/tableS2_cv_indication_terms.csv`. This file is the proposal as reviewed.

**Status: PROPOSED 2026-10-10, awaiting approval. No class counts have been computed.**
Pair counts and phase labels only. The term list assumes the extended set in `CV_SET_RULE.md`
(option 1: 62 terms, 1,879 pairs); it will be regenerated if you choose another option.

Reproduce:

```powershell
cd D:\causal2drug
python manuscript_prep\verify\table3_grouping_proposal.py   # writes table3_grouping_proposal.csv
```

The script holds the mapping and stops with an error if any term is left unmapped. The CSV has
one row per MeSH term with its ID, Minikel areas, pair count, Phase III and launched pair
counts, both proposed area columns and a note. It also carries the 31 "group B" terms, marked
as not proposed.

---

## The seven areas you listed do not fit two large blocks of the data

| block | pairs | where the seven areas would put it |
|---|---|---|
| Stroke, Brain Ischemia and other cerebrovascular terms | 380 | "Other" |
| Myocardial infarction, angina, acute coronary syndrome, restenosis | 257 | nowhere obvious; "Lipids / atherosclerosis" at a stretch |

With exactly seven areas, "Other" would hold 781 of 1,879 pairs (42 %) and
"Lipids / atherosclerosis" would be 46 % coronary events. I therefore propose **nine areas**:
your seven plus **Ischaemic heart disease** and **Cerebrovascular**. The CSV carries both
schemes (`proposed_area_9`, `area_if_only_7`) so either can be used.

| area | terms | pairs | launched pairs | | with only 7 areas: terms | pairs |
|---|---|---|---|---|---|---|
| Lipids / atherosclerosis | 8 | 311 | 20 | | 13 | 559 |
| Ischaemic heart disease *(added)* | 6 | 257 | 27 | | — | — |
| Heart failure / cardiomyopathy | 6 | 198 | 11 | | 6 | 198 |
| Arrhythmia | 8 | 66 | 11 | | 8 | 66 |
| Hypertension | 1 | 66 | 10 | | 1 | 66 |
| Pulmonary hypertension | 1 | 105 | 6 | | 1 | 105 |
| Thrombosis | 4 | 92 | 16 | | 6 | 104 |
| Cerebrovascular *(added)* | 8 | 380 | 8 | | — | — |
| Other | 20 | 404 | 16 | | 27 | 781 |
| **Total** | **62** | **1,879** | **125** | | 62 | 1,879 |

## Every term, by proposed area

Pairs in parentheses. **Bold** marks the 18 terms the extended rule adds.

**Lipids / atherosclerosis (8):** Atherosclerosis (132), Hyperlipidemias (76),
Hypercholesterolemia (47), Hypertriglyceridemia (31), Hyperlipoproteinemia Type II (12),
Dyslipidemias (7), Hyperlipoproteinemia Type I (5), Hypoalphalipoproteinemias (1).

**Ischaemic heart disease (6):** **Myocardial Infarction (129)**, Coronary Restenosis (45),
**Angina Pectoris (33)**, Acute Coronary Syndrome (28), **Angina, Unstable (13)**, Coronary
Thrombosis (9).

**Heart failure / cardiomyopathy (6):** Heart Failure (151), Cardiomyopathies (22),
Endomyocardial Fibrosis (9), Cardiomyopathy, Hypertrophic (9), **Edema, Cardiac (5)**,
Myocarditis (2).

**Arrhythmia (8):** **Atrial Fibrillation (33)**, **Arrhythmias, Cardiac (16)**,
**Ventricular Fibrillation (6)**, **Tachycardia, Ventricular (4)**,
**Tachycardia, Supraventricular (3)**, **Long QT Syndrome (2)**, **Tachycardia (1)**,
**Bradycardia (1)**.

**Hypertension (1):** Hypertension (66).

**Pulmonary hypertension (1):** Hypertension, Pulmonary (105).

**Thrombosis (4):** Thrombosis (58), Venous Thrombosis (20), Pulmonary Embolism (11),
Thrombophlebitis (3).

**Cerebrovascular (8):** Brain Ischemia (173), Stroke (173), **Cerebral Hemorrhage (13)**,
**Subarachnoid Hemorrhage (10)**, **Cerebral Infarction (6)**, Intracranial Thrombosis (3),
Vasospasm, Intracranial (1), **CADASIL (1)**.

**Other (20):** Cardiovascular Diseases (170), Peripheral Vascular Diseases (56),
**Reperfusion Injury (49)**, Vasculitis (26), Retinal Vein Occlusion (20), Hypotension (16),
**Intermittent Claudication (11)**, Raynaud Disease (6), Vascular Diseases (6), Retinal Artery
Occlusion (6), Optic Neuropathy, Ischemic (6), Venous Insufficiency (6), Aneurysm (6), Aortic
Valve Stenosis (4), Capillary Leak Syndrome (4), Thromboangiitis Obliterans (4),
Erythromelalgia (4), Retinal Telangiectasis (2), Endocarditis (1), Pericarditis (1).

## Placements that could reasonably go the other way

| term | proposed | alternative | why it is a judgement |
|---|---|---|---|
| Coronary Thrombosis (9) | Ischaemic heart disease | Thrombosis | it is both |
| Coronary Restenosis (45) | Ischaemic heart disease | Other | a complication of intervention, not a disease |
| Intracranial Thrombosis (3) | Cerebrovascular | Thrombosis | it is both |
| Hypotension (16) | Other | Hypertension | blood-pressure disorder, but the pharmacology is the reverse; 12 of its 16 targets were Class A in the original set, so where it goes moves the Hypertension row |
| Myocarditis (2) | Heart failure / cardiomyopathy | Other | inflammatory |
| Reperfusion Injury (49) | Other | Ischaemic heart disease | not organ-specific in MeSH |
| Peripheral Vascular Diseases (56) + Intermittent Claudication (11) | Other | a tenth area, "Peripheral arterial" | 67 pairs; your call |
| Cardiovascular Diseases (170) | Other | drop from Table 3, report in a footnote | a generic parent term, the third-largest in the set |

## Counting rules I propose for the table itself

1. **Unit is the target, counted once per area.** A target with several indications in one area
   counts once there. This also removes the Stroke / Brain Ischemia duplication, since the two
   have identical target lists.
2. **A target with indications in several areas appears in each.** Rows therefore do not sum to
   the set total; the table will say so.
3. **Phase within an area** follows the launched definition in `DEFINITIONS.md`, applied to the
   area: the highest phase over the target's indications in that area.

## What I need from you

1. Nine areas or seven.
2. Any of the placements above you want changed.
3. Whether "Cardiovascular Diseases" stays in "Other" or moves to a footnote.
4. Confirmation of the three counting rules.
