# Data dictionary

Supplementary bundle for the cardiovascular target untestability atlas. Every column of every shipped file is documented here. The counts, permitted values and ranges below are **computed from the shipped files at build time**, so they cannot drift from the data.

## Provenance classes

Every column carries one of three labels:

| class | meaning |
|---|---|
| **MEASURED** | An experimental measurement from an external resource, carried in unchanged. We did not compute it. |
| **DERIVED** | Computed by this project's code from other columns or from a raw source, under a rule stated in the column entry. |
| **ANNOTATION** | A categorical label curated by an external resource and carried in unchanged. |

A **DERIVED** column is only as good as its rule; the rule is given in full for each one, and the code that applies it is named in README.md.

---

## File 1 — the atlas

`cardiovascular_atlas.csv` — **630 rows x 25 columns**

Provenance classes in this file: **MEASURED** 5, **DERIVED** 16, **ANNOTATION** 4

### `gene`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | Minikel genetic_support v1.0 pp.tsv |
| in this file | 630 / 630 populated; 630 distinct; e.g. A2M, ABCA1, ABCC8 |

HGNC gene symbol of the drug target, as given by the source. Primary key; unique.

### `ensembl_id`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | Ensembl GRCh37.87 GTF (2026-08-26) |
| in this file | 618 / 630 populated; 618 distinct; e.g. ENSG00000004478, ENSG00000005381, ENSG00000005421 |

Symbol matched to a unique Ensembl gene record. Empty when no unique GRCh37 record exists (these rows have state = UNMAPPED).

### `cv_primary`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | Minikel areas.tsv x indic_topl_match.tsv |
| in this file | 630 / 630 populated; values: `False`, `True` |

TRUE if >=1 of the target's indications maps to MeSH top-level C14 (cardiovascular). The strict cardiovascular definition.

### `cv_extended`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | as cv_primary, plus a fixed token list |
| in this file | 630 / 630 populated; values: `True` |

TRUE if cv_primary, OR an indication term matches one of {hyperlipid, hypercholesterol, dyslipid, hyperlipoprotein, hypertriglycerid, lipoprotein, atheroscler, arterioscler}. This is the atlas definition; TRUE for every row.

### `n_cv_pairs`

| | |
|---|---|
| type | integer (count) |
| class | **DERIVED** |
| source | Minikel pp.tsv (deduplicated by ti_uid) |
| in this file | 630 / 630 populated; range: 1 to 13 |

Number of distinct target-indication pairs for this target that satisfy cv_extended.

### `max_clinical_phase`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | Minikel pp.tsv ccatnum |
| in this file | 630 / 630 populated; values: `Launched`, `Phase I`, `Phase II`, `Phase III`, `Preclinical` |

Highest clinical phase reached across this target's cardiovascular indications.

### `max_phase_num`

| | |
|---|---|
| type | integer (ordinal 1-5) |
| class | **DERIVED** |
| source | as max_clinical_phase |
| in this file | 630 / 630 populated; values: `1`, `2`, `3`, `4`, `5` |

Numeric encoding of max_clinical_phase: 5 = Launched, 4 = Phase III, 3 = Phase II, 2 = Phase I, 1 = Preclinical. Sorts correctly; use this for filtering.

### `approved_drug_direction`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | Open Targets 26.06 drug_mechanism_of_action |
| in this file | 171 / 630 populated; values: `DOWN`, `UP` |

UP or DOWN: whether the approved drug agonises or antagonises the target. Empty where the mechanism is absent or ambiguous (459 of 630 rows).

### `state`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | eQTLGen phase I cis-eQTL (2026-08-26) |
| in this file | 630 / 630 populated; values: `ABSENT`, `INSTRUMENTED`, `TESTED_NULL`, `UNMAPPED` |

Instrument availability in whole blood. ABSENT = gene not present in the eQTLGen file at all. TESTED_NULL = gene tested but no cis-eQTL reaches p < 5e-8. INSTRUMENTED = >=1 genome-wide significant cis-eQTL survives LD clumping. UNMAPPED = no unique GRCh37 gene record, so the question was never asked.

### `untestability_class`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | this study |
| in this file | 233 / 630 populated; values: `CLASS_A`, `CLASS_B` |

CLASS_A = ABSENT or TESTED_NULL (no instrument exists). CLASS_B = instrumented, but the drug mechanism operates at a level total mRNA cannot represent. Empty for instrumented targets with no Class B assignment and for UNMAPPED rows. CLASS_B is assigned at N = 1 (IL6R) and is under-counted by construction; see LIMITATIONS.md and classB_upper_bound.csv.

### `n_independent_instruments`

| | |
|---|---|
| type | integer (count) |
| class | **DERIVED** |
| source | eQTLGen + 1000 Genomes phase 3 EUR LD |
| in this file | 386 / 630 populated; range: 1 to 30 |

Independent genome-wide significant cis-eQTLs after LD clumping at r2 < 0.01 within +/-100 kb of the gene. Empty unless state = INSTRUMENTED.

### `best_eqtl_p`

| | |
|---|---|
| type | float (p-value) |
| class | **MEASURED** |
| source | eQTLGen phase I cis-eQTL |
| in this file | 386 / 630 populated; range: 3.2717e-310 to 4.0974e-08 |

Smallest cis-eQTL p-value for this gene in whole blood. Carried in unchanged. Empty unless state = INSTRUMENTED.

### `instrumented_blood`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | eQTLGen phase I (n = 31,684) |
| in this file | 630 / 630 populated; values: `False`, `True` |

TRUE if state = INSTRUMENTED. The convenient-tissue reference for the yield comparison.

### `instrumented_heart_LV`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | eQTL Catalogue R8 QTD000266 (GTEx heart left ventricle, n = 481) |
| in this file | 630 / 630 populated; values: `False`, `True` |

TRUE if the permuted-file top cis p-value < 5e-8. UPPER BOUND: uses a +/-1 Mb window (wider than the pre-registered +/-100 kb) and is unclumped.

### `instrumented_heart_atrial`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | eQTL Catalogue R8 (GTEx heart atrial appendage, n = 473) |
| in this file | 630 / 630 populated; values: `False`, `True` |

As instrumented_heart_LV. Upper bound.

### `instrumented_artery_aorta`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | eQTL Catalogue R8 (GTEx artery aorta, n = 492) |
| in this file | 630 / 630 populated; values: `False`, `True` |

As instrumented_heart_LV. Upper bound.

### `instrumented_liver`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | eQTL Catalogue R8 (GTEx liver, n = 269) |
| in this file | 630 / 630 populated; values: `False`, `True` |

As instrumented_heart_LV. Upper bound.

### `instrumented_cardiac_union`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | the four columns above |
| in this file | 630 / 630 populated; values: `False`, `True` |

TRUE if instrumented in ANY of heart LV, heart atrial, artery aorta or liver. Upper bound, inheriting the window caveat.

### `tissue_restricted`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | Human Protein Atlas RNA tissue specificity |
| in this file | 630 / 630 populated; values: `NA`, `NO`, `YES` |

YES = tissue-enriched, group-enriched or tissue-enhanced. NO = low tissue specificity. NA = gene absent from the HPA table. Used as the exposure in the Class A enrichment test.

### `secreted`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | Human Protein Atlas secretome annotation |
| in this file | 630 / 630 populated; values: `NA`, `NO`, `YES` |

YES = annotated as secreted. NO = not annotated as secreted. NA = gene absent from the HPA table.

### `gtex_blood_tpm`

| | |
|---|---|
| type | float (TPM) |
| class | **MEASURED** |
| source | GTEx v10 median TPM by tissue |
| in this file | 610 / 630 populated; range: 0 to 2997.85 |

Median transcripts per million in whole blood. Carried in unchanged. Empty where the gene is absent from the GTEx matrix (20 of 630 rows).

### `gtex_heart_LV_tpm`

| | |
|---|---|
| type | float (TPM) |
| class | **MEASURED** |
| source | GTEx v10 median TPM by tissue |
| in this file | 610 / 630 populated; range: 0 to 39841.1 |

Median TPM in heart left ventricle. Carried in unchanged.

### `gtex_liver_tpm`

| | |
|---|---|
| type | float (TPM) |
| class | **MEASURED** |
| source | GTEx v10 median TPM by tissue |
| in this file | 610 / 630 populated; range: 0 to 23227.2 |

Median TPM in liver. Carried in unchanged.

### `gtex_artery_aorta_tpm`

| | |
|---|---|
| type | float (TPM) |
| class | **MEASURED** |
| source | GTEx v10 median TPM by tissue |
| in this file | 610 / 630 populated; range: 0 to 10811.6 |

Median TPM in artery aorta. Carried in unchanged.

### `required_modality`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | this study (declared heuristic) |
| in this file | 630 / 630 populated; 32 distinct; e.g. UNKNOWN, cardiac eQTL, liver eQTL |

For CLASS_A targets, the measurement that could in principle instrument them. Rule, applied in order: cardiac eQTL if heart LV or atrial TPM >= 1; else liver eQTL if liver TPM >= 1; else vascular eQTL if aorta TPM >= 1; else pQTL if HPA-secreted; else 'tissue eQTL (<max-TPM tissue>)'; else UNKNOWN. IL6R is assigned from its established mechanism. This says where a signal COULD exist, not that one does. Non-CLASS_A rows carry the same rule's output for completeness.

> `cardiovascular_atlas.parquet` holds the identical content with types applied. It is verified at build time to match the CSV in row count and column order.

---

## File 2 — the Class B upper bound

`classB_upper_bound.csv` — **196 rows x 8 columns**

Provenance classes in this file: **MEASURED** 1, **DERIVED** 3, **ANNOTATION** 4

### `gene`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | as atlas.gene |
| in this file | 196 / 196 populated; 196 distinct; e.g. A2M, ABCA1, ABL1 |

Instrumented cardiovascular target carrying >=1 protein-altering GWAS association.

### `max_clinical_phase`

| | |
|---|---|
| type | string |
| class | **DERIVED** |
| source | as atlas.max_clinical_phase |
| in this file | 196 / 196 populated; values: `Launched`, `Phase I`, `Phase II`, `Phase III`, `Preclinical` |

Highest clinical phase for this target's cardiovascular indications.

### `n_protein_altering_assoc`

| | |
|---|---|
| type | integer (count) |
| class | **DERIVED** |
| source | GWAS Catalog associations, ontology-annotated full release |
| in this file | 196 / 196 populated; range: 1 to 825 |

Number of associations whose MAPPED_GENE contains this gene as an exact token AND whose CONTEXT is protein-altering.

### `top_variant`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | GWAS Catalog SNPS column |
| in this file | 196 / 196 populated; 196 distinct; e.g. rs1042522, rs1042714, rs1045485 |

rsID of the strongest such association (smallest p-value).

### `consequence`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | GWAS Catalog CONTEXT column |
| in this file | 196 / 196 populated; values: `frameshift_variant`, `inframe_insertion`, `missense_variant`, `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`, `stop_gained` |

VEP functional consequence of top_variant, as computed and versioned by the Catalog. Counted as protein-altering: missense_variant, stop_gained, frameshift_variant, start_lost, stop_lost, inframe_deletion, inframe_insertion, protein_altering_variant, splice_donor_variant, splice_acceptor_variant, coding_sequence_variant.

### `trait`

| | |
|---|---|
| type | string |
| class | **ANNOTATION** |
| source | GWAS Catalog MAPPED_TRAIT column |
| in this file | 196 / 196 populated; 153 distinct; e.g. 5-methyluridine (ribothymidine) measurement, ABL1/DNMBP protein level ratio in blood, ACPA-positive rheumatoid arthritis, rheumatoid arthritis, rheumatoid factor seropositivity measurement |

Reported trait for top_variant.

### `p_value`

| | |
|---|---|
| type | float (p-value) |
| class | **MEASURED** |
| source | GWAS Catalog P-VALUE column |
| in this file | 196 / 196 populated; range: 0 to 6e-06 |

p-value of top_variant. Carried in unchanged.

### `cv_relevant_trait`

| | |
|---|---|
| type | boolean |
| class | **DERIVED** |
| source | this study (fixed token list) |
| in this file | 196 / 196 populated; values: `False`, `True` |

TRUE if >=1 of this gene's protein-altering associations is for a cardiovascular trait, by fixed token match on the trait string. TRUE rows are the primary bound (16.1%); all 196 rows are the loose bound (50.8%).

---

## Sources, versions and access dates

| source | version | retrieved (UTC) | licence position |
|---|---|---|---|
| eQTLGen phase I cis-eQTL | 2019-12-11 release | 2026-08-26 | Open summary statistics; no explicit licence on the download page |
| eQTL Catalogue | R8 (QTD000131/136/141/251/256/266) | 2026-08-26 / 2026-08-27 | EMBL-EBI; open access |
| GWAS Catalog associations (full, ontology-annotated) | v1.0.3.1 release | 2026-08-31 | EMBL-EBI terms of use; open access |
| GTEx median TPM by tissue | v10 | 2026-08-30 (file mtime) | Open access for summary-level median expression |
| Human Protein Atlas | proteinatlas.tsv | 2026-08-30 (file mtime) | CC BY-SA 3.0 |
| Open Targets | 26.06 | 2026-08-30 (file mtime) | CC0 1.0 |
| Minikel genetic_support | v1.0 (GitHub tag) | 2026-08-26 | Repository licence |
| Ensembl GRCh37 GTF | GRCh37.87 | 2026-08-26 | Ensembl: no restriction on use |
| 1000 Genomes phase 3 / high coverage | 20130502 / 2504-sample panels | 2026-08-26 | Open access |
| UCSC liftOver chain | hg19ToHg38 | 2026-08-27 | Free for academic use |

**Licence caveat.** The licence column states the position as commonly documented for each resource. **None was re-verified against the source's own terms page during this build — treat the column as UNVERIFIED** and check before redistributing.

**Retrieval caveat.** URL, timestamp, byte count and SHA-256 are logged in `logs/retrieval_log.csv` for 22 of the raw inputs. Four — GTEx v10 median TPM, Human Protein Atlas, the GWAS Catalog associations archive and the Open Targets 26.06 extract — are **not in that log**; their retrieval URLs are **UNVERIFIED** and the dates above are filesystem mtimes, not logged retrieval times.
