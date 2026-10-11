# drug_links — target → drug → indication from ChEMBL, for manual review

> **MANUAL_REVIEW.csv — the file to review (round 3, 2026-10-10).**
>
> A strict rule was added beside the permissive one. **Strict:** the same approved drug must
> (1) be linked directly to this gene — ChEMBL target type SINGLE PROTEIN or PROTEIN COMPLEX,
> or the mRNA map — and (2) have an approved indication whose MeSH descriptor is under C14 or a
> lipid term **and** sits in no disease branch other than C14, C23, C10, C08 (lipid terms may
> also be in C18, C16). **Permissive:** any link, any approved indication under C14 or a lipid
> term.
>
> | Minikel launched | permissive | strict |
> |---|---|---|
> | Primary set (66) | 58 | 35 |
> | S1 (52) | 46 | 32 |
> | S2 (55) | 49 | 35 |
>
> `MANUAL_REVIEW.csv` has **327 rows for the 32 targets** verified by the permissive rule and
> not by the strict one (23 launched, 9 Phase III). One row per drug × approved cardiovascular
> indication, sorted so the strongest evidence for each target comes first.
>
> | column | use |
> |---|---|
> | `chembl_link_types`, `link_is_direct` | how ChEMBL attaches the drug to the gene |
> | `approved_indication`, `mesh_tree_numbers`, `indication_strictly_cv` | the indication and where MeSH files it |
> | `why_row_fails_strict` | the reason, per row |
> | `same_indication_as_minikel_launched` | True when this is the indication Minikel labels Launched |
> | `REVIEWER_accept_as_approved_cv_drug_for_this_target`, `REVIEWER_note` | **empty — yours to fill** |
>
> Two groups, which need different judgements:
> - **Link (20 targets):** ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, CACNA1C, EDNRA,
>   EDNRB, GUCY1A1, GUCY1B1, PDE4A, PDE4B, PDE4C, PDE4D, TUBB; Phase III: ADORA1, CACNA2D1,
>   GRIN1, SCN9A. The drug is real and cardiovascular; ChEMBL records it against the receptor
>   family or complex group. The question is whether a non-selective drug verifies each
>   subtype gene.
> - **Indication (12 targets):** C5AR1, EPAS1, GHR, GLA, KLKB1, SERPING1, VEGFA; Phase III:
>   BDKRB2, F11, F12, IL6R, SLC10A2. The link is direct; the indication is filed by MeSH under
>   a non-cardiovascular branch as well (for example hereditary angioedema, Turner syndrome,
>   vasculitis).
>
> Per-target result and reason: `chembl_target_summary.csv`, columns `verified_strict`,
> `verified_permissive`, `strict_drugs`, `strict_fail_reason`. Per-row flags in
> `chembl_target_drug_indication.csv`: `link_direct`, `ind_strict_cv`, `row_strict`.
> **No headline number is final until this file is reviewed.**

> **Rebuilt 2026-10-10 with antisense and RNAi drugs.** A hand-made map from ChEMBL's 62 antisense/RNAi mechanisms to genes (`verify/antisense_mechanisms.tsv`; map in `verify/build_chembl_links.py`) now links mipomersen to APOB, volanesorsen and olezarsen to APOC3, pelacarsen and olpasiran to LPA and inclisiran to PCSK9. **Current counts:** 8,151 link rows; 3,199 approved-indication rows, 2,359 flagged not cardiovascular; 1,688 mechanism rows, 1,397 drugs; 10 targets with no mechanism. Launched targets supported by an approved CV drug: 46 of 52 (S0), 49 of 55 (S0A), **58 of 66 (S0AB, the primary set)**. APOB and APOC3 are now supported; LPA has Phase 3 drugs only. Still unsupported in the primary set: AGT, EGFR, ITGB2, LPL, PPARG, RXFP1, SLC12A2, TNF. Counts in the body below predate this rebuild. Mismatches between Minikel labels and ChEMBL are in `tables/tableS6_phase_label_vs_chembl.csv`.

Built 2026-10-10 from **ChEMBL 37** (released 2026-05-01) through the ChEMBL REST API, with
MeSH tree numbers from the NLM MeSH RDF service. Every API response is cached under
`verify/cache/`, so the build re-runs offline and gives the same files.

```powershell
cd D:\causal2drug
python manuscript_prep\verify\build_chembl_links.py   # log saved in verify\chembl_build_output.txt
```

## Scope

Every target whose highest Minikel phase is Phase III or Launched under **any** of the three
candidate CV sets, so nothing needs rebuilding whichever set you approve. No instrument state
was read.

| candidate set | targets at Phase III or launched |
|---|---|
| S0 — original 630 set | 103 |
| S0A — original + rule R3 (`CV_SET_RULE.md`, option 1) | 115 |
| S0AB — original + every indication Minikel tags cardiovascular (option 2) | 141 |

Each file carries the target's Minikel phase under all three sets (`minikel_phase_S0`,
`_S0A`, `_S0AB`); a blank means the target is below Phase III in that set.

## Files

| file | rows | what it is |
|---|---|---|
| `chembl_target_summary.csv` | 141 | **Start here.** One row per target: how many drugs act on it, how many are approved, how many have an approved cardiovascular indication, their names, and a `review_priority`. |
| `chembl_target_drug_indication.csv` | 8,129 | One row per target × drug × indication, for indications at phase 3 or approved. This is the file with the flag you asked for. |
| `chembl_targets_without_mechanism.csv` | 14 | Targets for which ChEMBL has no mechanism row. |
| `chembl_build_info.json` | — | ChEMBL version and build time. |

## How the link is made

gene symbol → ChEMBL human targets that contain that gene (single protein, complex or family;
`target_type` says which) → `mechanism` rows → parent molecule → `drug_indication` rows.

An indication is treated as cardiovascular (`ind_is_cv`) when its MeSH descriptor has any tree
number under **C14**, or its heading matches one of the atlas's eight lipid tokens. This is the
same idea as the CV set rule, applied to ChEMBL's MeSH terms. Three further columns say whether
the descriptor is one of the Minikel indications in S0, S0A or S0AB.

**`FLAG_approved_not_cv`** is true when ChEMBL records the indication as approved
(`max_phase_for_indication` = 4) and it is not cardiovascular.

## What came out

| | count |
|---|---|
| Targets in scope | 141 |
| …with at least one ChEMBL human target | 138 |
| …with no ChEMBL mechanism row | 14 |
| Mechanism rows / distinct drugs | 1,674 / 1,384 |
| Rows in the link file | 8,129 across 121 targets |
| …approved-indication rows | 3,191 |
| …of which **flagged: approved, not cardiovascular** | **2,359** |

Does ChEMBL support Minikel's "Launched" label with an approved drug that has an approved
cardiovascular indication?

| set | launched targets | yes | phase 3 CV indication only | no CV indication | no mechanism in ChEMBL |
|---|---|---|---|---|---|
| S0 | 52 | 44 | 2 | 2 | 4 |
| S0A | 55 | 47 | 2 | 2 | 4 |
| S0AB | 66 | 56 | 3 | 3 | 4 |

Not supported in S0 and S0A: **AGT, APOB, APOC3, ITGB2, LPL, PPARG, RXFP1, SLC12A2**.
S0AB adds EGFR and TNF.

## The launched targets that Phase 1 listed as untestable

For the 21 launched targets with no blood instrument in the original set:

| result | targets |
|---|---|
| Approved drug with an approved CV indication (17) | PCSK9, ANGPTL3, NPC1L1, AGTR1, REN, EDNRA, EDNRB, F2, PLG, ITGA2B, ADRA1A, ADRA1B, ADRA1D, ADRA2B, ADRA2C, SLC5A1, IL1B |
| No mechanism row in ChEMBL (2) | APOB, APOC3 |
| Phase 3 CV indication only (1) | RXFP1 (serelaxin) |
| No approved or phase 3 CV indication (1) | AGT |

Three things in this table correct or sharpen what I wrote in `HEADLINE.md`:

* **IL1B.** I marked it questionable. ChEMBL lists rilonacept as approved for **pericarditis**,
  which is cardiovascular, so the target is supported. The Minikel indication carrying the
  Launched label is Atherosclerosis, and the only approved cardiovascular row ChEMBL has for
  an IL-1β drug is pericarditis (canakinumab and gevokizumab reach phase 3 for a
  cardiovascular indication). So the target passes and the specific target–indication pair
  does not.
* **RXFP1 and AGT** stay unsupported, as I suspected.
* **APOB and APOC3 are a limit of the link, not evidence against the label.** Their drugs are
  antisense oligonucleotides, and ChEMBL attaches those mechanisms to an RNA target, not to the
  protein the gene symbol resolves to. They need to be checked by hand.

For the four launched targets that were unmapped: ACE, ITGB3, GUCY1A1 and GUCY1B1 are all
supported (ACE inhibitors; abciximab, eptifibatide, tirofiban; nitrates and sGC stimulators).

## Phase III labels that look too low

ChEMBL has an approved drug with an approved CV indication for these targets, which Minikel
labels Phase III in the CV set:

* S0: ADORA2A, AVPR1A, AVPR2, DPP4, GRIN1, **HCN4**, IL6R, PDE3A, SLC6A2
* S0A adds: ADORA1, CACNA2D1, F11, SCN9A

HCN4 is the clear case: ivabradine is approved for heart failure and angina, dronedarone for
atrial fibrillation. Some of the others rest on a weak ChEMBL row (see below) and need your eye.

## Known weaknesses — why this needs a manual pass

1. **The match is at the level of "any cardiovascular indication", not the same indication.**
   A target counts as supported if any approved drug acting on it has any approved C14
   indication. Use `ind_in_S0` / `ind_in_S0A` to tighten that to the Minikel indications.
2. **ChEMBL indication rows are noisy.** Examples in this build: amantadine approved for
   "Intracranial Arteriosclerosis" (which is what makes GRIN1 look supported), and
   dextroamphetamine for "Arteriosclerosis". I have not removed any row.
3. **Complex and family targets widen the net.** A drug annotated against a receptor family is
   linked to every member gene. `target_type` and `n_chembl_targets_for_drug_gene` show this;
   the five adrenoceptor genes share most of their drug lists for this reason.
4. **Non-selective drugs inflate counts.** ADRA1A has 35 approved drugs, including carvedilol
   and epinephrine. Counts of drugs per target are not counts of drugs developed for that
   target.
5. **Oligonucleotide and gene-directed drugs are missed** (APOB, APOC3, and LPA at Phase III).
6. **`max_phase_for_indication` is ChEMBL's curation**, with no country or date. "Approved"
   means approved somewhere, at some time; `drug_withdrawn` marks withdrawals.
7. **Three targets have no ChEMBL human target at all** by gene symbol (among them CNTFR and
   LPA).

## Suggested review order

1. `chembl_target_summary.csv`, `review_priority = HIGH` (10 rows): launched in Minikel, no
   approved CV drug found in ChEMBL.
2. `review_priority = MEDIUM` (27 rows): Phase III in Minikel, no approved or phase 3 CV
   indication found.
3. In the link file, filter `FLAG_approved_not_cv = True` for the targets you care about.
