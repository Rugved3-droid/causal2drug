# DECISIONS_LOG — decisions taken after Phase 1

All decisions are the project owner's, taken 2026-10-10 in two rounds. This file is the record
for `manuscript_prep/`. The project's own `DECISIONS.md` and `PROTOCOL.md` are not edited.

## Round 1 (after the Phase 1 review)

| # | decision | outcome | where |
|---|---|---|---|
| 1 | Extend the CV set; keep the original 630 as a sensitivity analysis | rule proposed, then approved in round 2 | `CV_SET_RULE.md` |
| 2 | Resolve gene-record mapping for unmapped targets, ACE first | done; ACE instrumented | `UNMAPPED_RESOLUTION.md` |
| 3 | Define "tested-null" with code; one "launched" definition; reconcile 21/52 and 17/44 | done | `DEFINITIONS.md` |
| 4 | Class B is no longer a class; IL6R is a worked example | done | below |
| 5 | Target → drug → indication link from ChEMBL | built | `drug_links/` |
| 6 | Correct or remove the "96.9 % of cardiovascular loci" figure | removed in round 2 | `COLOC_CLAIM_CORRECTION.md` |
| 7 | Indication-area grouping for Table 3 | proposed, then approved in round 2 | below |
| 8 | Dated deviation note for the Class C figure; supplementary only | done | `DEVIATION_CLASS_C_FIGURE.md` |

## Round 2 (the three gates)

| decision | implemented as |
|---|---|
| **CV set.** Option 2 is primary: every indication Minikel tags cardiovascular (93 terms, 726 targets). Reperfusion Injury stays. No post-hoc exclusions. The 31 extra terms are listed and non-cardiac ones flagged, not removed. Original 630 = sensitivity 1. | `analysis/classify.py` `SETDEF`. Sensitivity 2 (663 targets) and sensitivity 3 (C14 only, 594) were added by me; see the note under "Open points". |
| **Naming.** Class A becomes "uninstrumentable", with sub-types "no eQTL data" and "no qualifying instrument", reported separately everywhere. The sum is an explicitly labelled upper bound; no-eQTL-data alone is the lower bound. | Status values `NO_EQTL_DATA`, `NO_QUALIFYING_INSTRUMENT`; every table has both columns plus an "upper bound" column. |
| Supplement carries the three no-instrument sub-types. | `tables/tableS5_no_qualifying_instrument_subtypes.*` |
| **Cis-window sensitivity** at ± 500 kb and ± 1 Mb; ± 100 kb stays primary. | `classify.py` [K07]; `tables/tableS4_sensitivity_cis_window.*` |
| **Launched** = per-set definition. Correct the "20" in all docs. | Definition in `classify.py`. The 20 was corrected in three existing Markdown files; see "Edits to existing files". |
| **Drop the old Class B bound.** Modality mismatch is one sentence with the IL6R example. | Not computed, not tabulated. `classB_upper_bound.csv` stays in the old bundle, unused. |
| **Table 3.** Nine areas plus an "Unspecified CV" row for the generic "Cardiovascular Diseases" term. Coronary Thrombosis → ischaemic heart disease. Other judgement calls assigned by the disease the term names; every placement logged. | `classify.py` `AREA`; log below. |
| **Gene records.** HGNC/Ensembl fixes approved; MIR29B1 and MIR92A1 named as unresolved. | Applied to all 60 unmapped symbols. **MIR195 is a third unresolved target in the primary set.** |
| **Drugs.** Add APOB and APOC3 antisense drugs. Do not edit Minikel labels. Mismatches go in a supplementary table. Headline claims about approved CV drugs use the ChEMBL-verified count. | Manual mRNA-target map in `verify/build_chembl_links.py`; `tables/tableS6_phase_label_vs_chembl.*`; `KEY_NUMBERS.md` §4. |
| **Colocalization 96.9 % claim:** remove entirely. | It appears in no new table, figure or outline. The old `table2_method_availability.*` is not carried forward. |

## Class B retired

* One category of untestable target remains. The paper does not use "Class A / Class B".
* IL6R is **instrumented** and counted as such in every table.
* Modality mismatch is one sentence in the Discussion. Facts available for it
  (`verify_numbers.py` block [V08]): 226 eQTLGen variants at p < 5e-8 within ± 100 kb; MR used
  rs2229238 (exposure p = 4.27e-31, F = 134.5); estimate +0.376 (SE 0.023), p = 1.16e-59; the
  pre-registered expectation was negative, based on the shedding variant rs2228145, so the
  control failed on direction and stays failed.

## Table 3 placement log

Rule for judgement calls: **the disease the term names, not the mechanism.** The full
term-to-area list is `tables/tableS2_cv_indication_terms.csv`.

| term | area | basis |
|---|---|---|
| Cardiovascular Diseases | Unspecified CV | your instruction |
| Coronary Thrombosis | Ischaemic heart disease | your instruction |
| Coronary Restenosis | Ischaemic heart disease | names coronary disease |
| Intracranial Thrombosis | Cerebrovascular | names an intracranial vascular disease; same logic as Coronary Thrombosis |
| Myocarditis | Heart failure / cardiomyopathy | names a disease of the myocardium |
| Edema, Cardiac | Heart failure / cardiomyopathy | names cardiac oedema |
| Hypotension | Other | names hypotension, which is not hypertension |
| Reperfusion Injury | Other | names no organ |
| Peripheral Vascular Diseases; Intermittent Claudication | Other | peripheral arterial disease; no such area was requested |
| Cardiomyopathy, Dilated; Diabetic Cardiomyopathies | Heart failure / cardiomyopathy | name a cardiomyopathy |
| Dementia, Vascular; Cerebral Amyloid Angiopathy; CADASIL | Cerebrovascular | name a cerebrovascular disease |
| Cerebral Hemorrhage; Subarachnoid Hemorrhage; Cerebral Infarction; Vasospasm, Intracranial | Cerebrovascular | as named |
| Tetralogy of Fallot; Ductus Arteriosus, Patent; Heart Defects, Congenital | Other | congenital structural heart disease; none of the nine areas names it |
| Aortic Valve Stenosis; Endocarditis; Pericarditis | Other | valvular, infective and pericardial disease; no area names them |
| Vasculitis and the named vasculitides (Behçet, Takayasu, Churg-Strauss, Granulomatosis with Polyangiitis, Microscopic Polyangiitis, Mucocutaneous Lymph Node Syndrome, Thromboangiitis Obliterans) | Other | vasculitis |
| Retinal Vein Occlusion; Retinal Artery Occlusion; Retinal Telangiectasis; Optic Neuropathy, Ischemic; Diabetic Retinopathy | Other | ocular vascular disease |
| Diabetic Foot; Varicose Ulcer; Venous Insufficiency; Raynaud Disease; Erythromelalgia; Aneurysm; Arteriovenous Malformations; Vascular Diseases; Capillary Leak Syndrome; Telangiectasia, Hereditary Hemorrhagic | Other | peripheral or unspecified vascular disease |
| Angioedemas, Hereditary; Angioedema; Ehlers-Danlos; Fabry; Turner; MELAS; Noonan; Marfan; von Hippel-Lindau; Ataxia Telangiectasia; Alagille; Sarcoglycanopathies; Glycogen Storage Disease IIb | Other | multisystem syndromes |
| Thrombosis; Venous Thrombosis; Pulmonary Embolism; Thrombophlebitis | Thrombosis | as named |

Result: "Other" holds 46 of the 93 terms. That follows from having no area for peripheral,
ocular, congenital, valvular or vasculitic disease.

Counting rules, as proposed and not objected to: a target is counted once per area; a target
with indications in several areas appears in each; phase within an area is the highest phase
over the target's indications in that area.

## Edits to existing files

The original brief said to modify no existing file. Round 2 said to "correct the '20' in all
docs". I treated the later, specific instruction as authorising that one correction and nothing
else.

| file | change |
|---|---|
| `deliverables/cv/ANALYSES_1_2_3.md` | table cell 20 → "17 (corrected 2026-10-10; was 20)"; one sentence corrected with a dated note |
| `deliverables/cv/CARDIOVASCULAR_ATLAS.md` | "20 vs 21" → "17 of 44 vs 21 of 52" with a dated note |
| `deliverables/cv/LIMITATIONS.md` | table cell 20 → "17 (corrected 2026-10-10; was 20)" |

Byte-for-byte originals are in `manuscript_prep/originals_before_correction/`.

**Not edited:** `deliverables/cv/cv_robustness.json` still holds `"classA_launched": 20` for
CV-PRIMARY. It is the output of `src/a2_robustness.py`, and editing it by hand would make it
disagree with the script that writes it. Fixing the script is a code change you have not asked
for.

## Open points

1. **"Both sensitivity sets."** Only sensitivity 1 was named. I used the 663-target set as
   sensitivity 2 and added C14-only as sensitivity 3.
2. **MIR195** is unresolved in addition to the two you named.
3. **X-linked genes.** eQTLGen tested none. 29 of the 217 no-eQTL-data targets in the primary
   set are X-linked (`KEY_NUMBERS.md` §6). They are correctly "no eQTL data", but the reason is
   the resource's design and not blood expression. Whether to report them as a third named
   sub-type is your call; for now they are a column in Supplementary Table S1 and a table of
   their own (S10).
4. **Launched examples.** GHR, GLA and PDE4C are in the primary launched list through Turner
   syndrome, Fabry disease and Behçet syndrome. I recommend naming examples only from the 17
   that survive in sensitivity 2.

## Round 3 (2026-10-10)

| decision | implemented as |
|---|---|
| Keep S1 (630), S2 (663) and S3 (C14 only) as sensitivity sets, all in the supplement. | Table S3; `classify.py` `SETDEF`. |
| X-linked is a third named sub-type, inside "no eQTL data". Report the autosomal-only lower bound beside the overall one everywhere. | Column `no_eqtl_data_subtype`; every bound table carries "autosomal gene not tested", "X-linked, not analysed" and "Autosomal-only lower bound"; Figure 1 and Figure 2 show the split. Definition: autosomal no-eQTL-data targets ÷ targets not on X, Y or MT. |
| Leave `cv_robustness.json`; mark it superseded in `README.md`. | Done. |
| Add a strict drug-verification rule; report strict and permissive side by side; output disagreements as `drug_links/MANUAL_REVIEW.csv`; do not finalise headline numbers before review. | `verify/build_chembl_links.py` "strict rule"; `KEY_NUMBERS.md` §4; drug-verified rows in every table and figure are labelled provisional. |
| Examples in tables and figures only from strict-verified targets; never GHR, PDE4C or other non-cardiac additions. | Column `example_eligible` = strict-verified and same phase in S2. Figure 3 draws only eligible targets. Complete listings (Tables S1, S6, S7) still contain every target, with the flag. |

**How I read "strict".** The instruction was: "a target counts as verified only if a linked
approved drug carries an approved cardiovascular indication". The permissive rule already
required an approved drug with an approved cardiovascular indication, so a literal reading
changes nothing. I tightened the two loose words. *Linked*: the drug must be annotated to this
gene directly (single protein, specific complex, mRNA), not through a family or group target.
*Cardiovascular*: the indication's MeSH descriptor must sit only in cardiovascular branches.
This is my interpretation and is flagged as open conflict N1 in `KEY_NUMBERS.md`.

**Effect.** Launched verified targets in the primary set fall from 58 to 35. Sixteen drop for
the link (adrenoceptor and endothelin subtypes, sGC subunits, CACNA1C, PDE4 subtypes, TUBB) and
seven for the indication (GHR, GLA, KLKB1, SERPING1, C5AR1, EPAS1, VEGFA).

**Example eligibility is a rule of mine too.** "The other non-cardiac additions" was made
mechanical as "not at the same phase in S2", which removes every target that owes its phase to
one of the 31 terms filed outside C14 and signs/symptoms.

## Round 4 — 2026-10-11 (recovery map and relabel)

| # | Decision (owner) | What was done |
|---|---|---|
| 4.1 | Fenland workbook supplied by the owner | All four recovery-map resources read; `manuscript/DATA_ACCESS.md` updated |
| 4.2 | Symbol match in GTEx only after a documented Ensembl ID change | `manuscript/RECOVERY_RULE.md` amendment 8 (post-hoc). ADRA2B locked → fixable-tissue; GDF2 and MIR145 stay locked. Evidence: `manuscript/recovery_out/gtex_id_change_resolution.csv` |
| 4.3 | Relabel the sub-type "X-linked (27) and mitochondrial (2)" | `analysis/classify.py`, `tables/make_tables.py`, `figures/make_figures.py`, `KEY_NUMBERS.md`, `LIMITATIONS.md`, `METHODS_FACTS.md`, `OUTLINE.md`, `figures/README.md` edited; classifier, tables and Figures 1–2 regenerated. Previous versions in `originals_before_relabel/`. No count changes; "14 of 29 expressed in blood" becomes "12 of 27 X-linked and 2 of 2 mitochondrial" |
| 4.4 | Reporting tiers: plausible rescue / off-target-tissue rescue / locked (never assayed; assayed, no instrument) | `manuscript/RECOVERY_RULE.md` amendment 9 (post-hoc, reporting only) |
| 4.5 | Figure 3 becomes the recovery map; Table 2 the 13 named targets | The former Figure 3 (expression of launched and Phase III targets) is now Supplementary Figure S2; the former Tables 2 and 3 move to the supplement |

## Round 5 — 2026-10-11 (revision before submission)

| # | Decision (owner) | What was done |
|---|---|---|
| 5.1 | Positive controls with candidate instruments, pre-specified | `manuscript/RESCUE_CONTROLS_PROTOCOL.md` written first; 2 pass (PCSK9, APOC3), 1 fail (AGTR1), 2 not harmonizable (ANGPTL3, REN); nothing substituted |
| 5.2 | Audit the IL6R control; do not reclassify | `manuscript/IL6R_AUDIT.md`: harmonization correct; instrument r² 0.16 with rs2228145; control remains a fail |
| 5.3 | "Candidate instrument availability" language | Title, text, tier display labels, Tables 2–3, S13–S17, Figure 3. Keys in `recovery_map.csv` unchanged |
| 5.4 | "No evidence of a phase gradient", with CIs | `analysis/phase_ci.py`; Table S21 |
| 5.5 | Heart failure as a measurement-coverage gap | Results and Discussion |
| 5.6 | IL1B not a named example; Table 2 site-of-action and approval-status columns; flag APOB | Figure S2 redrawn without IL1B; Table 2 and S18 from FDA and EMA records; APOB and PLG flagged |
| 5.7 | Zenodo preprint, official STROBE-MR wording, GTEx v8 sample sizes | `manuscript/PRIOR_ART.md`; supplement checklist; Tables S16, S20 |
