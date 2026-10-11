# METHODS_FACTS — what a reviewer will ask for

> **Round 3, 2026-10-10.** Numbers about approved drugs are **provisional** until `drug_links/MANUAL_REVIEW.csv` is reviewed. A strict verification rule now sits beside the permissive one; "no eQTL data" has two named sub-types (autosomal, not tested; X-linked, not analysed) and the autosomal-only lower bound is reported with the overall one. Current figures: `KEY_NUMBERS.md`.

Each fact names the file it comes from. **[MISSING]** marks something not found in the
repository. **[UNVERIFIED]** marks something stated in a project document that I did not
re-check. Retrieval times are UTC, from `logs/retrieval_log.csv` unless noted.

## 1. Data sources

| source | version / file | retrieved | size, sha256 (first 16) | used for |
|---|---|---|---|---|
| Minikel et al. target–indication table | `ericminikel/genetic_support`, tag **v1.0**, `data/pp.tsv` | 2026-08-26 02:52 | 3,162,501 B, `db1bda5367759fb2`; git blob SHA-1 `16d641186c39…` matches GitHub per `DELIVERABLE_0.md:179-196` | targets, indications, phase labels |
| Minikel indication maps | same tag: `areas.tsv`, `indic.tsv`, `indic_topl_match.tsv`, `mesh_2023_topl_maps.tsv` | 2026-08-26 23:15 | in log | cardiovascular set (`indic.tsv` column `areas`) |
| eQTLGen phase I cis-eQTL, full | `2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | 2026-08-26 23:16 | 3,880,187,026 B, `1133a89fe33b795d` | tested-gene universe (19,250 genes; 127,341,798 rows) |
| eQTLGen phase I cis-eQTL, FDR < 0.05 | `2019-12-11-cis-eQTLsFDR0.05-…txt.gz` | 2026-08-26 23:15 | 322,775,879 B, `8d963046d7b74cf3` | instruments (10,507,664 rows; max p 2.019e-05) |
| eQTLGen allele frequencies | `2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz` | 2026-08-27 23:36 | 240,045,342 B, `edffd43094598d08` | β and SE from Z (controls only) |
| Ensembl gene annotation | `Homo_sapiens.GRCh37.87.gtf.gz` | 2026-08-26 23:23 | 42,400,095 B, `56bc520f5e9cbc14`; matches Ensembl `CHECKSUMS` (`DELIVERABLE_1B.md`) | gene coordinates |
| HGNC | REST `rest.genenames.org/fetch/symbol/…` | 2026-10-10 | cached in `verify/cache/hgnc/` (60 symbols) | resolving unmapped symbols. HGNC release identifier **[MISSING]** — the REST service returns none |
| 1000 Genomes phase 3 (GRCh37) | release 20130502 VCFs; panel `integrated_call_samples_v3.20130502.ALL.panel` | 2026-08-26 23:19 | panel 55,156 B, `b4023dc6ee2d62ee` | LD for clumping; 503 EUR samples |
| 1000 Genomes high coverage (GRCh38) | `20130606_g1k_3202_samples_ped_population.txt` | 2026-08-26 22:40 | 97,483 B | LD for control and Class C analyses; 525 unrelated EUR **[UNVERIFIED count]** |
| eQTL Catalogue | release **R8**; GTEx v10 permuted files QTD000256 (heart LV, n = 481), QTD000251 (atrial appendage, 473), QTD000131 (aorta, 492), QTD000266 (liver, 269) | 2026-08-26 / 27 | in log | tissue instrument yield (upper bounds) |
| eQTL Catalogue | QTD000356, GTEx v10 whole blood, n = 853 | **[MISSING]** not in retrieval log; queried by range | — | Class C comparison arm |
| GTEx median expression | `gtex_v10_median_tpm.gct.gz`, v10 | file mtime 2026-08-30; **URL [MISSING]** | 8,846,936 B | blood and tissue TPM |
| Human Protein Atlas | `proteinatlas.tsv.zip` | file mtime 2026-08-30; **URL and HPA version [MISSING]** | 7,460,135 B | tissue specificity, secretome |
| GWAS Catalog studies | `gwas-catalog-download-studies-v1.0.3.1.txt` | 2026-08-26 02:49 | 129,139,393 B, `c66e0f409aca64ef` | outcome GWAS selection |
| GWAS Catalog harmonised list | `harmonised_list.txt` | 2026-08-26 02:47 | 11,959,625 B, `0b85b1d792096280` | which studies have harmonised files |
| UCSC liftOver chain | `hg19ToHg38.over.chain.gz` | 2026-08-27 23:26 | 227,698 B, `5c0598e500ceb5a7` | GRCh37 → GRCh38 for MR |
| ChEMBL | **release 37** (2026-05-01), REST API | 2026-10-10 | cached in `verify/cache/chembl/` | drug–target–indication link |
| MeSH | NLM RDF `id.nlm.nih.gov/mesh/<ID>.json` | 2026-10-10 | cached in `verify/cache/mesh/` (984 descriptors); **MeSH year [MISSING]** | C14 membership of ChEMBL indications |
| Open Targets | 26.06 `drug_molecule`, `drug_mechanism_of_action` | file mtime 2026-08-30; **URL [MISSING]** | — | **not used in the manuscript analysis** (Part C and the old Table 1 only) |

eQTLGen sample size: n = 31,684, as stated throughout the project documents and consistent with
the per-variant `NrSamples` values in the control files (for example 31,684 at rs6453133).
Number of cohorts and the consortium citation: **[MISSING]** from the repository.

Citation for the Minikel dataset: **[MISSING]** from the repository. A web search returns
Minikel et al., "Refining the impact of genetic evidence on clinical success", *Nature* 2024,
doi 10.1038/s41586-024-07316-0; confirm against the repository README before citing.

## 2. Genome build

* Target classification is done entirely on **GRCh37**: eQTLGen positions and the Ensembl
  release-87 GTF.
* The two MR controls lift eQTLGen variants to **GRCh38** to meet GWAS Catalog harmonised
  files: `pyliftover` 0.4.1 with the UCSC chain above. 503 variants lifted, 0 failed
  (`DELIVERABLE_1C.md`). A variant that fails is dropped and counted. Genome-wide failure rate
  **[UNVERIFIED]**.

## 3. Study set

| step | rule | N | source |
|---|---|---|---|
| Raw rows | `pp.tsv` | 29,476 | `classify.py` [K01] |
| Deduplicate | one row per `ti_uid`, first row kept | 29,460 | `DECISIONS.md` decision 1 |
| Exclude | all rows of 4 `ti_uid` with conflicting highest phase: DAO-D022124, HRH1-D000086382, IFNAR2-D012141, MC1R-D010787 | 29,456 | frozen file sha256 `5cb7b8015bdc4756…` |
| Cardiovascular, primary | Minikel `areas` contains `cardiovascular`, or MeSH term contains one of 8 lipid tokens | 93 indications, 2,128 pairs, 726 targets | `CV_SET_RULE.md` |
| Sensitivity sets | 630, 663 and 594 targets | | `CV_SET_RULE.md` |

Three of the retained duplicate keys differ in non-outcome columns (`DECISIONS.md:49-54`); no
such column is used.

Phase: Minikel `ccat` (combined historical and active highest phase), coded Preclinical,
Phase I, Phase II, Phase III, Launched. A target's phase in a set is its highest phase over its
indications in that set.

## 4. Gene records

1. A symbol maps if it matches exactly one `gene` feature on a primary chromosome in Ensembl
   GRCh37 release 87 (`src/d1b_gene_coords.py`). 2,457 of 2,517 map.
2. For the other 60, take the Ensembl gene ID of the approved HGNC record and accept it if it
   is a primary-chromosome gene in the same GTF. 45 resolve; 15 do not.
3. Unresolved in the primary set: MIR29B1, MIR92A1, MIR195.

The join to eQTLGen is by Ensembl gene ID, not symbol (eQTLGen symbols are Ensembl v71
nomenclature, `DELIVERABLE_1B.md` §3).

## 5. Instrument definition

| parameter | value | source |
|---|---|---|
| Exposure resource | eQTLGen phase I, whole blood, cis-eQTL | `DECISIONS.md` decision 7 |
| Significance | p < 5 × 10⁻⁸, frozen | decision 6 |
| Cis-window, primary | gene body ± 100,000 bp (GRCh37) | pre-registration R1, R12.6 |
| Cis-window, sensitivity | ± 500 kb, ± 1 Mb | added 2026-10-10 |
| eQTLGen's own cis range | variants within 1 Mb of the gene centre **[UNVERIFIED against eQTLGen methods]**; largest gene-body distance seen for any target is 812,919 bp | `classify.py` [K07] |
| Tested-gene universe | 19,250 genes; all with a GRCh37.87 record are autosomal (19,123); 127 IDs are absent from GRCh37.87 | [K10] |

**Status.** Instrumented: ≥ 1 variant meeting the criterion. No eQTL data: Ensembl ID not in the
tested universe, with sub-types — autosomal gene not tested; X-linked (27) and mitochondrial (2),
not analysed by eQTLGen (relabelled 2026-10-11; no primary-set target is Y-linked). Autosomal-only lower bound: autosomal no-eQTL-data targets ÷ targets
not on X, Y or MT. No qualifying instrument: tested, no variant meets the criterion. Code:
`analysis/classify.py`, function `status`; original: `src/partb_classA.py:50-58`.

**Clumping** (used for instrument counts in the existing atlas and for MR, not for status):
greedy, lowest p first, drop r² ≥ 0.01; LD from 1000 Genomes phase 3 GRCh37, 503 EUR samples,
by HTTP range query (`src/d1b_clump.py`). Strength criterion F ≥ 10, with F = z².

Sub-types of "no qualifying instrument": T1, a p < 5e-8 variant exists outside the window; T2,
the gene has a row in the FDR < 0.05 file but none at p < 5e-8; T3, no row in that file.

## 6. Annotations

* **Tissue restricted:** Human Protein Atlas "RNA tissue specificity" in {Tissue enriched,
  Group enriched, Tissue enhanced}; not restricted = "Low tissue specificity"; others excluded.
* **Secreted:** HPA "Secretome location" non-empty, or "Protein class" contains "secreted".
* **Blood expression:** GTEx v10 `Whole_Blood` median TPM.
* **Other-tissue instruments:** top nominal cis p < 5e-8 in the eQTL Catalogue permuted file,
  ± 1 Mb, unclumped. These are upper bounds.
* Tests: Fisher exact; odds ratio with Woolf 95 % CI (Haldane correction if a cell is zero).

## 7. Drug verification

ChEMBL 37. Gene symbol → human targets containing the gene → `mechanism` → parent molecule →
`drug_indication`. Antisense and RNAi drugs are linked through a hand-made map from ChEMBL's 62
such mechanisms (`verify/antisense_mechanisms.tsv`). Scope: targets at Phase III or launched.
Minikel labels were not changed.

| rule | drug | link to the gene | approved indication (`max_phase_for_ind` 4) |
|---|---|---|---|
| Permissive | `max_phase` 4 | any ChEMBL target containing the gene | MeSH tree number under C14, or a lipid token |
| **Strict** | `max_phase` 4 | target type SINGLE PROTEIN or PROTEIN COMPLEX, or the mRNA map; not PROTEIN FAMILY, PROTEIN COMPLEX GROUP, SELECTIVITY GROUP, CHIMERIC PROTEIN, PROTEIN-PROTEIN INTERACTION | as permissive, and every disease branch of the descriptor is one of C14, C23, C10, C08 (lipid terms may also sit in C18, C16) |

The same drug must satisfy all three conditions. Which rule the paper uses is **undecided**
pending review of `drug_links/MANUAL_REVIEW.csv`. A target is named as an example only if it is
strict-verified and holds the same phase in sensitivity set 2.

## 8. Positive controls (MR)

| | value | source |
|---|---|---|
| Estimator | lead-SNP Wald ratio; SE = \|SE_out / β_exp\| | pre-registration R12.2 |
| Exposure β, SE | from Z by the SMR approximation (Zhu et al. 2016) using eQTLGen allele frequencies | `DELIVERABLE_1C.md` |
| LDL-C outcome | GCST90239658, Graham et al. 2021, European, N = 1,320,016 | `DELIVERABLE_1C.md` |
| CRP outcome | GCST90179146, Koskeridis et al. 2022, N = 575,531 | `DELIVERABLE_1C.md` |
| Outcome selection | EFO match, harmonised file with `.tbi`, European-only, largest N | R2, R2a |
| Harmonisation | match on chromosome and position; allele flip; strand flip; palindromic SNPs resolved by frequency only if both < 0.42 or both > 0.58, else dropped; indels dropped | R6; `src/harmonize.py`; 10 of 10 branch tests pass |
| Pass criteria | instrument at p < 5e-8, F ≥ 10, MR p < 0.05, sign as pre-registered | R10 |
| Original resource | GTEx v10 liver (QTD000266, n = 269): no control had an instrument | `DELIVERABLE_1.md` |

Outcome GWAS "per trait" applies only to these two controls. The atlas itself uses no outcome
GWAS: it classifies targets by instrument availability.

## 9. Class C test

Pre-registered in `PREREGISTRATION.md` R14.4 and `PROTOCOL.md:62-72`. 571 Part C pairs → 402 with
≥ 100 overlapping variants → 106 one per target → 60 sampled, seed 20260830 → 50 paired.
coloc.abf implemented directly (Giambartolomei et al. 2014); priors p1 = p2 = 1e-4, p12 = 1e-5;
prior effect SD 0.15. Exposures: GTEx v10 whole blood (n = 853) and eQTLGen (n = 31,684), same
outcome, window and harmonisation. Tests: Wilcoxon signed-rank on log10 PP4; exact binomial.

## 10. Software

| | version | note |
|---|---|---|
| Python | 3.12.10 | this machine, 2026-10-10 |
| numpy / scipy / pandas | 2.4.2 / 1.17.1 / 3.0.1 | |
| matplotlib / openpyxl | 3.11.1 / 3.1.5 | figures, workbook |
| pyliftover / requests / pyarrow | 0.4.1 / 2.32.5 / 25.0.1 | used by `src/` |
| tabix, pysam | not used | `src/tabix_http.py` is a pure-Python tabix client |
| coloc, SuSiE | not the R packages | reimplemented in `src/` |

These are the versions installed now. **Whether they are the versions that produced the August
outputs is [UNVERIFIED]:** the repository has no `requirements.txt`, lock file or environment
export, and is not under version control. The manuscript pipeline (`manuscript_prep/`) was run
with exactly these versions.

## 11. Pre-registration and deviations

| item | where | note |
|---|---|---|
| Pre-registration | `deliverables/d1/PREREGISTRATION.md`, R1–R15 | dated in the text; **no external timestamp** (no registry, no commit history) |
| Resource changed from GTEx liver to eQTLGen | amendment R3a, decision 7 | before any re-run |
| Threshold frozen at 5e-8 after a control was seen to pass at 1e-5 | decision 6 | disclosed |
| PCSK9 retired, not replaced | decision 8 | |
| Class C refutation criterion | `PROTOCOL.md:62-72` | written before the test |
| Class C reported despite the no-colocalization rule | `DEVIATION_CLASS_C_FIGURE.md` | 2026-10-10 |
| CV set, sub-type names, window sensitivity, gene-record rule, ChEMBL verification | `DECISIONS_LOG.md` | all post hoc, 2026-10-10; not pre-registered |

## 12. Reproduction

```powershell
cd D:\causal2drug
python manuscript_prep\verify\build_chembl_links.py   # ChEMBL link (offline from cache)
python manuscript_prep\analysis\classify.py           # classification, key numbers
python manuscript_prep\tables\make_tables.py          # all tables
python manuscript_prep\figures\make_figures.py        # all figures
```

Raw inputs are not in `manuscript_prep/`; paths are absolute to `D:\causal2drug`. Licences of
the incorporated sources are **[UNVERIFIED]** (`deliverables/supplementary/README.md:46-50`).
Data and code availability statement: **[MISSING]** — no public repository exists yet.
