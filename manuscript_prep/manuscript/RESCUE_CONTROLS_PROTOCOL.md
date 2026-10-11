# RESCUE_CONTROLS_PROTOCOL — positive controls using candidate instruments from the recovery map

**Written 2026-10-11, before any outcome association for any of these variants was retrieved
and before any Mendelian randomization estimate was computed.** What was known at the time of
writing: the recovery map (`recovery_map.csv`), that is, which variant is the lead qualifying
cis variant for each target in each resource and its exposure p-value; the column headers of the
exposure tables; and the outcome-selection and harmonization rules of the main pipeline
(`deliverables/d1/PREREGISTRATION.md`, R2, R2a, R6, R10, R12.2). This file is dated, not
externally timestamped; the project is not under version control. It is not edited after
results are seen; anything that has to change is added under "Amendments", dated, with the
reason.

## Question

Does a candidate instrument that the recovery map found outside eQTLGen whole blood reproduce
the known direction of effect of an approved drug on its target? Five controls, set by the
project owner.

## Controls and expected direction

The exposure is genetically predicted **higher** level of the target (plasma protein, or tissue
expression). Each drug lowers the level or blocks the action of its target, so the expected
sign of the MR estimate is the opposite of the drug's effect on the outcome.

| # | target → outcome | approved drug and its effect | expected sign of MR β (per unit higher target) |
|---|---|---|---|
| 1 | PCSK9 → LDL cholesterol | PCSK9 inhibitors lower LDL cholesterol | **positive** |
| 2 | ANGPTL3 → triglycerides | evinacumab (ANGPTL3 inhibitor) lowers triglycerides | **positive** |
| 3 | APOC3 → triglycerides | volanesorsen (lowers APOC3) lowers triglycerides | **positive** |
| 4 | AGTR1 → systolic blood pressure | angiotensin receptor blockers lower blood pressure | **positive** |
| 5 | REN → systolic blood pressure | aliskiren (renin inhibitor) lowers blood pressure | **positive** |

Stated in advance: for REN the exposure is circulating renin *concentration*, which rises by
feedback when the system is blocked or pressure falls, and for AGTR1 it is expression in one
tissue. Either could plausibly give the "wrong" sign for reasons other than a broken pipeline.
That does not change the criterion: a wrong sign is a fail.

## Instrument: taken from the recovery map, no new selection

Rule, same as the main analysis and the recovery map: a variant with p < 5 × 10⁻⁸ within the
gene body ± 100 kb in the resource's own build.

- **Which resource.** If the target's recovery-map rescue is a plasma cis-pQTL, the pQTL is
  used; if it is a cardiovascular-relevant-tissue eQTL only, the eQTL is used. Where a target
  has both, the plasma pQTL is used (it is the same molecule the drug acts on or lowers).
  By the recovery map this gives: PCSK9, ANGPTL3, APOC3, REN → pQTL; AGTR1 → tissue eQTL
  (it is assayed by none of the three proteomic studies).
- **Which pQTL study.** The largest study whose published table states the effect allele, the
  effect size and its standard error explicitly. UKB-PPP (alleles in the variant identifier,
  `BETA` with respect to A1, `SE`) and Fenland (`EA`, `NEA`, `Effect`, `SE`) do; the deCODE
  table gives minor/major alleles, a beta and a p-value but no standard error and no stated
  effect allele in its header, so deCODE is used only if neither other study has an
  instrument, and then only if its effect allele can be established from the table's own
  legend. Order: UKB-PPP, then Fenland, then deCODE.
- **Which variant.** Exactly the lead qualifying variant recorded for that study in
  `recovery_map.csv` (`<study>_lead_variant`; for UKB-PPP the sheet is in `ukbppp_note`).
- **Which tissue (AGTR1).** The cardiovascular-relevant GTEx v8 tissue with the lowest lead
  p-value, as recorded in `recovery_map.csv` (`gtex_best_cv_tissue`,
  `gtex_best_cv_variant_b38`). Exposure effect: GTEx `slope` and `slope_se`, which are for the
  ALT allele of the variant identifier.
- One variant per control. No variant is substituted if the recorded one cannot be used.

## Outcome GWAS: the main pipeline's rule, unchanged

Rule R2 with amendment R2a, applied by the pipeline's own function
(`src/select_gwas.py`, imported, with its output written under `manuscript/rescue_controls/`):
GWAS Catalog studies release v1.0.3.1 mapped to the trait's EFO term; harmonised summary
statistics listed; a tabix index present; European ancestry only; largest parsed sample size;
ties by lowest accession.

| outcome | EFO term | study |
|---|---|---|
| LDL cholesterol | EFO_0004611 | GCST90239658 (already selected by this rule for the main controls) |
| triglycerides | EFO_0004530 | whatever the rule returns |
| systolic blood pressure | EFO_0006335 | whatever the rule returns |

The outcome region is read by the pipeline's own tabix client (`src/tabix_http.py`). If the
selected file cannot be range-queried, the control is reported as "outcome not retrievable".
No other study is substituted.

## Harmonization: the main pipeline's code, unchanged

`src/harmonize.py`, function `harmonize_one` (rule R6), imported and called as the main
pipeline calls it: match on GRCh38 chromosome and position; aligned, allele-flipped and
strand-flipped cases handled; palindromic variants resolved by allele frequency only if both
frequencies are below 0.42 or both above 0.58 and agree, otherwise dropped; indels and
multi-allelic records dropped. Positions: UKB-PPP and GTEx are already GRCh38; Fenland is
GRCh37 and is lifted with the pipeline's chain file (`hg19ToHg38.over.chain.gz`, pyliftover),
a failed lift being a failure to harmonize. **If the variant is absent from the outcome file,
is an indel, mismatches on alleles or is an unresolved palindrome, the control is reported as
"not harmonizable" with the reason. It is not replaced by another variant, study or tissue.**

## Estimator

One variant per control, so the **Wald ratio**, as in the main pipeline (R12.2):
β = β_outcome / β_exposure; SE = |SE_outcome / β_exposure|; 95% CI = β ± 1.96 SE; two-sided p
from the normal distribution. F = (β_exposure / SE_exposure)². Inverse-variance weighting is
not used because one variant is taken per target. Units: outcome units of the selected GWAS
per unit of the exposure as scaled in its source table (reported with each estimate).

## Pass criterion (R10 of the main pipeline, with the eGene clause not applicable)

A control **passes** only if all hold: (1) the instrument has p < 5 × 10⁻⁸ within ± 100 kb in
its resource; (2) F ≥ 10; (3) MR p < 0.05; (4) the sign is the expected sign above. Anything
else is a **fail**, or "not harmonizable" / "outcome not retrievable" where no estimate
exists. Nothing is changed to turn a fail into a pass. All five are reported whatever happens.

## What is reported

For each control: resource, variant, alleles, exposure β, SE, p and F; outcome study, N,
β, SE, p; harmonization status; MR β, 95% CI, p; expected and observed sign; pass or fail.
Output: `manuscript/rescue_controls/` and `manuscript/RESCUE_CONTROLS_RESULTS.md`.

## What this does and does not show

A pass shows that this one candidate instrument reproduces a known direction. It does not
validate the other candidate instruments in the recovery map, and a fail does not show that
the target is irrelevant.

## Amendments

None at the time of writing.
