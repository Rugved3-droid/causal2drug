# HEADLINE — launched and Phase III cardiovascular targets in an untestable class

> **Superseded 2026-10-10 by `KEY_NUMBERS.md` §4.** On the approved primary set: 58 launched targets with a ChEMBL-verified approved cardiovascular drug; 20 (34.5 %) have no eQTL data and 24 (41.4 %) are at the upper bound; 19 to 23 beyond PCSK9. This file is the Phase 1 record for the original 630-target set.

> **Update 2026-10-10, after the Phase 1 review. This file describes the original 630-target
> set and has not been recomputed.** Four things below are now superseded:
> (1) Class B is no longer a class, so IL6R leaves the untestable counts (`DECISIONS_LOG.md`).
> (2) The four unmapped launched targets are resolved: ACE, ITGB3 and GUCY1A1 are instrumented,
> GUCY1B1 is tested-null, making it 22 of 52 (`UNMAPPED_RESOLUTION.md`).
> (3) My "questionable" note on IL1B was half wrong: ChEMBL lists rilonacept as approved for
> pericarditis, so the target is supported though the atherosclerosis pair is not. RXFP1 and AGT
> stay unsupported, and HCN4's Phase III label is confirmed too low (`drug_links/README.md`).
> (4) The CV set is being extended (`CV_SET_RULE.md`); every count here will be redone on it.

Reproduce every count and every row:

```powershell
cd D:\causal2drug
python manuscript_prep\verify\verify_numbers.py     # block [V05]; rows -> verify\headline_targets.csv
```

## The three answers

**How many launched CV drug targets are untestable by this pipeline, beyond PCSK9?**
**20.** Of the 52 CV targets carrying a Minikel "Launched" label, 21 are Class A; PCSK9 is one of
them. 17 of the 20 are absent from eQTLGen altogether and 3 are tested but have no cis-eQTL at
p < 5e-8.

**What share of launched CV targets is that?**

| numerator | of 52 launched | of the 48 with a gene mapping |
|---|---|---|
| Class A including PCSK9: 21 | **40.4 %** | 43.8 % |
| Class A beyond PCSK9: 20 | **38.5 %** | — |
| Absent from eQTLGen only: 18 | 34.6 % | 37.5 % |

Four further launched targets (ACE, ITGB3, GUCY1A1, GUCY1B1) have no class at all because the
gene symbol could not be matched to one GRCh37 record. Class B adds nothing to the launched
count: its single member, IL6R, reaches Phase III in a CV indication (vasculitis).

**The sharpest sentence the data support:**

> Of 52 cardiovascular drug targets with a launched drug, 21 (40 %) have no genome-wide
> significant cis-eQTL in the largest whole-blood eQTL resource, and 18 of those are not
> measured in it at all, so blood cis-eQTL Mendelian randomization cannot test PCSK9 or the
> targets of ezetimibe, angiotensin receptor blockers, aliskiren, endothelin antagonists and
> direct thrombin inhibitors.

## Three things that limit how hard that sentence can be pushed

1. **Launched targets are not worse off than other targets.** Class A is 40.4 % among launched
   targets, 36.8 % across all 630 CV targets, 32–38 % in every other phase, and 35.0 % across
   all 2,517 Minikel targets (881 of 2,517). The finding is "four in ten approved CV targets are
   invisible", not "approved CV targets are especially invisible".
2. **The count is of genes, not drugs.** Five of the 21 are α-adrenoceptor subtypes and two are
   endothelin receptors. The repository has no drug-to-indication link, so "how many approved
   drugs" cannot be answered from it (DISCREPANCIES D2).
3. **"Launched" is the Minikel label, taken as given.** It is the highest phase over a target's
   CV indications, historical or active, any market. For three of the 21 I could not match the
   label to an approved CV drug I know of (AGT, IL1B, RXFP1). That check comes from my general
   knowledge, not from the repository, and is marked as such below.

---

## Launched CV targets in Class A (21)

Columns: **state** is the eQTLGen state. **CV indications** are every CV-EXTENDED indication in
which the target carries the Launched label (the existing Table 1 shows only one of them).
**Open Targets agent** is the example drug the existing `table1_classA_launched_phase3.csv`
lists for the target; it is linked to the target, **not** to the indication. **Blood / liver /
heart** are GTEx v10 median TPM. **Note** is my own plausibility check and is **UNVERIFIED
against the repository**.

| gene | state | CV indications at Launched | Open Targets agent | blood / liver / heart LV | note (UNVERIFIED, mine) |
|---|---|---|---|---|---|
| **PCSK9** | ABSENT | Atherosclerosis; Hypercholesterolemia; Hyperlipidemias; Hyperlipoproteinemia Type II | Alirocumab | 0.29 / 28.3 / 0.0 | consistent |
| ANGPTL3 | ABSENT | Hyperlipoproteinemia Type II | Evinacumab | 0.02 / 127.7 / 0.0 | consistent |
| APOB | ABSENT | Hyperlipoproteinemia Type II | Mipomersen | 0.02 / 350.4 / 1.8 | consistent; mipomersen later withdrawn |
| APOC3 | ABSENT | Hyperlipoproteinemia Type I; Hypertriglyceridemia | Olezarsen | 0.79 / 8,122 / 0.1 | consistent |
| NPC1L1 | TESTED-NULL | Dyslipidemias | Ezetimibe | 0.03 / 31.1 / 0.1 | consistent |
| AGTR1 | ABSENT | Heart Failure; Hypertension | Losartan | 0.03 / 28.2 / 1.4 | consistent |
| REN | ABSENT | Hypertension | Aliskiren | 0.00 / 0.4 / 0.0 | consistent |
| AGT | ABSENT | Hypertension; Hypotension | GSK-2586881 | 0.16 / 312.6 / 24.1 | **questionable** — the listed agent is not a launched drug |
| EDNRA | ABSENT | Hypertension, Pulmonary | Bosentan | 0.02 / 1.1 / 6.2 | consistent |
| EDNRB | ABSENT | Brain Ischemia; Hypertension, Pulmonary; Stroke | Bosentan | 0.53 / 7.8 / 10.4 | consistent for pulmonary hypertension |
| F2 | ABSENT | Brain Ischemia; Pulmonary Embolism; Stroke; Thrombosis; Venous Insufficiency; Venous Thrombosis | Desirudin | 0.09 / 598.4 / 0.1 | consistent |
| PLG | ABSENT | Thrombosis | Aprotinin | 0.02 / 288.6 / 0.0 | target consistent (thrombolytics); listed agent is an antifibrinolytic |
| ITGA2B | TESTED-NULL | Coronary Restenosis; Coronary Thrombosis | Abciximab | 36.3 / 0.4 / 0.3 | consistent; well expressed in blood |
| ADRA1A | ABSENT | Hypotension | Prazosin | 0.01 / 13.3 / 2.2 | target consistent; listed agent treats hypertension, not hypotension |
| ADRA1B | ABSENT | Hypotension | Prazosin | 0.00 / 6.4 / 0.79 | same |
| ADRA1D | ABSENT | Hypotension | Prazosin | 0.01 / 0.04 / 0.13 | same |
| ADRA2B | ABSENT | Hypertension | not in Open Targets | no GTEx value | no agent retrievable; expression unknown (D13) |
| ADRA2C | ABSENT | Hypertension | Clonidine | 0.05 / 0.33 / 3.4 | consistent |
| SLC5A1 | ABSENT | Heart Failure | Sotagliflozin | 0.02 / 1.1 / 46.7 | consistent |
| IL1B | TESTED-NULL | Atherosclerosis | Rilonacept | 29.5 / 0.5 / 0.1 | **questionable** — I know of no IL-1β agent approved for atherosclerosis; well expressed in blood |
| RXFP1 | ABSENT | Heart Failure | Serelaxin | 0.00 / 0.6 / 0.0 | **questionable** — serelaxin is not approved in the US or EU as far as I know |

TPM values are from `deliverables/cv/cardiovascular_atlas.csv` (`gtex_blood_tpm`,
`gtex_liver_tpm`, `gtex_heart_LV_tpm`). None of the 177 absent CV targets appears in the eQTLGen
tested-gene list under its symbol either, so absence is not an identifier-join artefact.

## Phase III CV targets in an untestable class (18 Class A + 1 Class B)

`hcat|acat` shows whether the Phase III label is historical, active, or both.

| gene | class | state | CV indications at Phase III | Open Targets agent | historical \| active |
|---|---|---|---|---|---|
| LPA | A | ABSENT | Cardiovascular Diseases; Hyperlipidemias | Olpasiran | — \| Phase III |
| APOA1 | A | TESTED-NULL | Acute Coronary Syndrome | not in Open Targets | — \| Phase III |
| HCN4 | A | ABSENT | Heart Failure | Ivabradine | Phase III \| — |
| HCN1 | A | ABSENT | Heart Failure | not in Open Targets | Phase III \| — |
| HCN2 | A | TESTED-NULL | Heart Failure | not in Open Targets | Phase III \| — |
| PDE3A | A | TESTED-NULL | Brain Ischemia; Peripheral Vascular Diseases; Stroke | Enoximone | Phase III \| — |
| AVPR2 | A | ABSENT | Heart Failure | Lypressin | — \| Phase III |
| AVPR1A | A | ABSENT | Hypotension | Lypressin | Phase III \| — |
| ABCC8 | A | ABSENT | Brain Ischemia; Stroke | Glipizide | — \| Phase III |
| ENPEP | A | TESTED-NULL | Hypertension, Pulmonary | Firibastat | — \| Phase III |
| BRD4 | A | TESTED-NULL | Atherosclerosis | Pelabresib | — \| Phase III |
| IL17A | A | ABSENT | Vasculitis | Ixekizumab | — \| Phase III |
| KIT | A | TESTED-NULL | Hypertension, Pulmonary | Imatinib | — \| Phase III |
| PDGFRA | A | TESTED-NULL | Hypertension, Pulmonary | Cediranib | — \| Phase III |
| GRIN1 | A | TESTED-NULL | Brain Ischemia; Stroke | Ketamine | — \| Phase III |
| SLC6A2 | A | ABSENT | Brain Ischemia; Hypotension; Stroke | Doxepin | Phase III \| — |
| SLC6A4 | A | TESTED-NULL | Brain Ischemia; Hypotension; Stroke | Mazindol | Phase III \| — |
| CNTFR | A | ABSENT | Retinal Telangiectasis | not in Open Targets | — \| Phase III |
| **IL6R** | **B** | INSTRUMENTED | Vasculitis | (not in table1) | — \| Phase III |

Phase III totals: 51 targets, 18 Class A (35.3 %), 9 absent and 9 tested-null, plus IL6R.
Phase III or launched combined: 103 targets, 39 Class A (37.9 %); 40 when IL6R is added.

Two notes on this list. The strongest Phase III example for a cardiology reader is **LPA**
(blood 0.00 TPM, liver 14.6): the lipoprotein(a) programmes cannot be tested in blood. And
**HCN4 is labelled Phase III although ivabradine is, to my knowledge, an approved heart-failure
drug**, which suggests the Minikel phase labels err in both directions. That last point is
UNVERIFIED against the repository.

## What the hook does not yet survive

Each of these is detailed in `DISCREPANCIES.md`.

* The CV set excludes myocardial infarction, angina, atrial fibrillation and every arrhythmia
  indication, because Minikel files them under "signs/symptoms" (D16). A cardiology reviewer will
  see that at once.
* ACE is among the 4 unmapped launched targets, for a gene-record ambiguity (D13).
* Three of the 21 are tested-null, not absent, and two of those (ITGA2B, IL1B) are highly
  expressed in blood. The "not expressed in blood" explanation does not cover them (D12).
* Under the strict CV definition the launched Class A count is 17 of 44 (38.6 %), not the 20
  the existing documents report (D1). The share holds; the count does not.
