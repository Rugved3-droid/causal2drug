# COLOC_CLAIM_CORRECTION — the "96.9 % of cardiovascular target loci" figure

> **Decided 2026-10-10: remove entirely.** The claim and the colocalization block of the old Table 2 are not carried into any new table, figure or outline. No replacement sentence is used.

**Dated 2026-10-10. Recommendation: remove the figure. The existing files are not edited.**

Reproduce:

```powershell
cd D:\causal2drug
python manuscript_prep\verify\coloc_claim_check.py    # output saved in verify\coloc_claim_output.txt
```

## The claim as it stands

`deliverables/cv/CARDIOVASCULAR_ATLAS.md:154-158`, `cv/table2_method_availability.csv` rows 1-4,
`cv/cv_report.json`, and the rendered Table 2:

> Of 65 CV target loci screened, only 3 (4.6 %) have a genome-wide significant outcome signal in
> the cis-window, and only 2 (3.1 %) yield a usable credible-set colocalization.
> **Colocalization is unavailable for 96.9 % of cardiovascular target loci.**

## What the 65 loci are

`src/cv_atlas_report.py:96` keeps every R15 locus whose **gene** is in the atlas. It does not
look at the outcome. Each R15 locus pairs one target with the outcome GWAS that Part C selected
for one of that target's indications, which is often not a cardiovascular one.

| outcome trait of the 65 loci | loci |
|---|---|
| Asthma | 17 |
| **Heart Failure** | **12** |
| Idiopathic Pulmonary Fibrosis | 8 |
| Alzheimer Disease | 6 |
| Sepsis | 4 |
| Diabetes Mellitus | 3 |
| Sarcoidosis | 2 |
| 13 other traits, one locus each (among them Atrial Fibrillation and Myocardial Infarction) | 13 |

**53 of the 65 loci test a cardiovascular drug target against a non-cardiovascular disease.**
The three loci counted as having an outcome signal are ADORA1 / Atrial Fibrillation,
PPARG / Asthma and TNF / Asthma.

## The same screen, restricted to cardiovascular outcomes

| definition of "cardiovascular locus" | loci | outcome genome-wide significant | credible-set result |
|---|---|---|---|
| As published: gene is in the atlas, any outcome | 65 | 3 | 2 |
| Outcome is an original-set CV indication | **12** (all Heart Failure) | **0** | **0** |
| Outcome is a CV indication under the proposed extended rule | **16** (12 Heart Failure, 3 Atrial Fibrillation, 1 Myocardial Infarction) | **2** (ADORA1, KCNH2; both Atrial Fibrillation) | 2 |
| All 106 screened loci, any area | 106 | 6 | 4 |

## Why the figure should go, not be re-stated

1. **The corrected denominator is 12 or 16 loci from one, or three, outcome GWAS.** "0 of 12"
   describes the heart-failure GWAS used here. It is not a statement about cardiovascular
   target loci in general, and a percentage on 12 would overstate what is known.
2. **The screen was never designed as a cardiovascular survey.** R15 takes one locus per target,
   chosen by largest variant overlap, from pairs that already passed two Part C eligibility
   rules (an approved-drug direction and an exact trait match).
3. **It sits awkwardly with the standing rule** that colocalization output is not tabulated.
   The exception recorded in `DEVIATION_CLASS_C_FIGURE.md` covers the Class C test only.

## Proposed replacement

Remove the sentence, the "96.9 %" headline and the colocalization block of Table 2 from the
paper. If one sentence on the point is wanted, this is what the data support:

> In a descriptive screen of 106 instrumented target loci, 6 had a genome-wide significant
> outcome association in the cis-window. Of the 16 loci paired with a cardiovascular outcome,
> 2 did, both for atrial fibrillation. The screen covered too few cardiovascular outcomes to
> estimate how often colocalization is available for cardiovascular targets.

The second sentence depends on the extended CV set being approved. Under the original set it
reads "Of the 12 loci paired with a cardiovascular outcome (all heart failure), none did."

## What is not affected

* The infrastructure-attrition figures (265 of 448 CV pairs lost to unindexed GWAS files) are a
  separate result about file packaging and are not colocalization output. They were verified in
  Phase 1 and will be recomputed on the extended set.
* The Class C test and its numbers.

## Decision needed

Remove entirely (my recommendation), or keep the replacement sentence in the supplement.
