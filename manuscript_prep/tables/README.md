# Tables

```powershell
cd D:\causal2drug
python manuscript_prep\analysis\classify.py      # inputs
python manuscript_prep\tables\make_tables.py     # every table below
```

Each table exists as `.csv` (data), `.md` (formatted, with title and footnotes) and as a sheet
in `manuscript_tables.xlsx` (formatted, 13 sheets). Markdown versions of tables over 80 rows
show the first 40; the CSV and workbook hold all rows.

| file stem | content | rows |
|---|---|---|
| `table1_flow_and_cohort` | **Table 1.** Derivation of the set; composition under the primary and three sensitivity definitions | 24 |
| `table2_status_by_phase` | **Table 2.** Instrument status by highest phase, primary set; strict and permissive drug-verified rows | 10 |
| `table3_status_by_indication_area` | **Table 3.** Instrument status by indication area (nine areas plus Unspecified CV) | 11 |
| `tableS1_per_target_classification` | **S1.** Every target in the primary set: status, sub-type, chromosome, window sensitivity, phases per set, expression, tissue instruments, areas, ChEMBL result | 726 |
| `tableS2_cv_indication_terms` | **S2.** The 93 MeSH indications: rule, area, pair counts, flags | 93 |
| `tableS3_sensitivity_cv_set` | **S3.** Bounds overall, by phase and for verified launched targets (strict, permissive) under S1, S2, S3 | 33 |
| `tableS4_sensitivity_cis_window` | **S4.** Bounds at ± 100 kb, ± 500 kb, ± 1 Mb | 18 |
| `tableS5_no_qualifying_instrument_subtypes` | **S5.** The 63 targets with no qualifying instrument, by sub-type | 63 |
| `tableS6_phase_label_vs_chembl` | **S6.** Launched targets not strict-verified, and Phase III targets with an approved CV drug in ChEMBL | 46 |
| `tableS7_launched_targets` | **S7.** All 66 launched targets with status and ChEMBL drugs | 66 |
| `tableS8_gene_record_resolution` | **S8.** The 60 symbols without a unique gene record and their resolution | 60 |
| `tableS9_enrichment` | **S9.** Tissue restriction and secretion, by group and set | 16 |
| `tableS10_no_eqtl_data_by_chromosome` | **S10.** No-eQTL-data targets: autosomal versus X/Y/MT | 5 |

Not carried forward from the old bundle: `table2_method_availability` (colocalization, removed)
and `classB_upper_bound.csv` (Class B retired).

Checks built into the pipeline: the classifier stops if it cannot reproduce the shipped state
of every already-mapped target, if any primary-set term lacks an indication area, or if the
C14 branch is not contained in the cardiovascular tag.

## Conventions in every bound table (round 3)

- "No eQTL data" is followed by its two sub-types, "autosomal gene not tested" and "X-linked,
  not analysed", and by the **autosomal-only lower bound** (autosomal no-eQTL-data targets ÷
  targets not on X, Y or MT).
- Rows about approved drugs come in a **strict** and a **permissive** version and are labelled
  **provisional** until `drug_links/MANUAL_REVIEW.csv` is reviewed.
- `example_eligible` (Tables S1, S7) marks the only targets that may be named as examples.
