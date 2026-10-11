# manuscript_prep — index

> **Status, 2026-10-11: final.** Drug verification is complete (`drug_links/REVIEW_DECISION.md`, `manuscript/APPROVED_35_CHECK.csv`); nothing is provisional. Older notes below and in `METHODS_FACTS.md`, `LIMITATIONS.md` and `HEADLINE.md` that say otherwise predate this. A strict verification rule now sits beside the permissive one; "no eQTL data" has two named sub-types (autosomal, not tested; X-linked, not analysed) and the autosomal-only lower bound is reported with the overall one. Current figures: `KEY_NUMBERS.md`.

> **Round 4, 2026-10-11.** The manuscript draft, supplement, cover letter and recovery map are in `manuscript/` (start with `manuscript/OPEN_ITEMS.md`). The paper is now framed around the recovery map; the overall percentage is cited as known. The sub-type formerly called "X-linked (29)" is "X-linked (27) and mitochondrial (2)". The former Figure 3 is Supplementary Figure S2 and the former Tables 2 and 3 are Supplementary Tables S11 and S12. `OUTLINE.md`, `JOURNAL_FIT.md` and `PRIOR_ART_CHECK.md` are superseded.

Publication preparation for the cardiovascular drug-target atlas (project `D:\causal2drug`).
Built 2026-10-10. Start with `KEY_NUMBERS.md`.

## Read first

| file | what it is |
|---|---|
| `KEY_NUMBERS.md` | **The corrected numbers**, each with its reproducing command. Lists the conflicts found in the rerun. |
| `JOURNAL_FIT.md` | Which tier the evidence supports, and why. |
| `OUTLINE.md` | Article outline, bullets with numbers placed. |
| `LIMITATIONS.md` | Limitations, each with its measured size. |
| `METHODS_FACTS.md` | Sources, versions, parameters, software; gaps marked [MISSING]. |
| `PRIOR_ART_CHECK.md` | Overlapping work found by web search, and the limits of that search. |
| `DECISIONS_LOG.md` | Every decision taken after Phase 1, the Table 3 placement log, edits made to existing files, open points. |

## Results

| folder | content |
|---|---|
| `tables/` | Tables 1–3 and S1–S10 as CSV, Markdown and one formatted workbook; `README.md` |
| `figures/` | Figures 1–3 and Supplementary Figure S1 as PDF, SVG, PNG; `README.md` has the legends |
| `drug_links/` | ChEMBL 37 target → drug → indication tables; **`MANUAL_REVIEW.csv` is the file waiting for you** (327 rows, 32 targets where the strict and permissive rules disagree); `README.md` |
| `analysis/` | `classify.py`, its captured output, and `out/` (per-target tables, `key_numbers.json`) |

## Rules and definitions

| file | status |
|---|---|
| `CV_SET_RULE.md` | approved; the 93-term primary set and the 31 flagged-or-not additions |
| `DEFINITIONS.md` | approved; what "no qualifying instrument" means, with code; the launched definition |
| `UNMAPPED_RESOLUTION.md` | approved; gene-record fixes |
| `DEVIATION_CLASS_C_FIGURE.md` | dated deviation note for Supplementary Figure S1 |
| `COLOC_CLAIM_CORRECTION.md` | why the 96.9 % claim was removed |
| `TABLE3_GROUPING_PROPOSAL.md` | the proposal as reviewed; final placements are in `DECISIONS_LOG.md` |

## Phase 1 records (original 630-target set; superseded where marked)

`INVENTORY.md`, `NUMBERS.md`, `HEADLINE.md`, `DISCREPANCIES.md`.

## Reproduce everything

```powershell
cd D:\causal2drug
python manuscript_prep\verify\verify_numbers.py       # Phase 1 checks, controls, Class C
python manuscript_prep\verify\build_chembl_links.py   # ChEMBL link (uses verify\cache)
python manuscript_prep\analysis\classify.py           # classification and key numbers
python manuscript_prep\tables\make_tables.py
python manuscript_prep\figures\make_figures.py
```

The scripts read the project's raw data by absolute path and are not runnable from the zip
alone. `verify/cache/` (API responses, 36 MB) and `originals_before_correction/` are in the
project folder; the cache is left out of the zip.

## Superseded project file, left as it is

`deliverables/cv/cv_robustness.json` holds `"classA_launched": 20` for CV-PRIMARY. **That value
is superseded.** The correct figure under the per-set launched definition is 17 of 44 (as
shipped) or 18 of 44 after the gene-record fixes. The file is the output of
`src/a2_robustness.py` and was left unedited on instruction; do not cite the 20 from it.

## Existing project files changed

Three Markdown files under `deliverables/cv/` had the erroneous "20" corrected, with a dated
note in place. Originals are in `originals_before_correction/`. Nothing else outside this
folder was modified. See `DECISIONS_LOG.md`.
