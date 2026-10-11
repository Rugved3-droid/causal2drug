# INVENTORY — Causal2Drug / cardiovascular target untestability atlas

Project root: `D:\causal2drug` (not a git repository; 539 files before `manuscript_prep/` was added).
Status key: **FINAL** = feeds the paper as-is · **INTERMEDIATE** = an input to a final item, or a
record of a completed step · **STALE** = superseded, withdrawn, or contradicted by a later file.

Reproduce the listing: `Get-ChildItem D:\causal2drug -Recurse -File | Select FullName,Length,LastWriteTime`

Timestamps are local (UTC−4): `MANIFEST.json` records `built_utc 2026-09-01T00:52:56Z` for files
stamped 08-31 20:52. There are **no notebooks** (`*.ipynb`: 0 files).

Checkpoint directories (one small file per locus or chromosome) are listed as one line per
directory. Everything else is listed individually.

---

## 1. Root documents

| file | what it is | status |
|---|---|---|
| `DECISIONS.md` | Frozen decisions 1–13 plus Gate D notice. Dedup rule, p<5e-8, eQTLGen switch, PCSK9 retirement, class definitions. | FINAL (binding), but its last section predates Gate D closure — see `PROTOCOL.md` |
| `INSTRUMENT_AVAILABILITY.md` | Running record of instrument availability: genome-wide census, PCSK9 anchor, class definitions, `.tbi` attrition, Class C refutation. | **PARTLY STALE** — Class C section is current; the "Class B remains unquantified", "Lane B not yet run" and "2 of 28 loci (7 %)" passages were superseded on 08-31 (DISCREPANCIES D5, D7) |
| `PROTOCOL.md` | Gate status table, Gate D closure, pre-registered Class C refutation criterion and outcome. | FINAL for the Class C criterion and verdict; line 33 ("worsens as the exposure resource grows") is STALE (D8) |

## 2. `src/` — scripts (41)

| script | what it does | status |
|---|---|---|
| `cv_atlas.py` | Builds the 630-row atlas, class counts, tissue yield, enrichment; writes `cv/cardiovascular_atlas.csv`, `cv_summary.json`. | FINAL |
| `cv_atlas_report.py` | Rewrites the parquet; launched/Phase III Class A lists, modality summary, CV method availability; writes `cv_report.json`. | FINAL |
| `a1_classB_bound.py` | Class B upper bound from GWAS Catalog protein-altering associations. | FINAL |
| `a2_robustness.py` | Every headline result under CV-PRIMARY vs CV-EXTENDED. | FINAL |
| `a3_cascade.py` | Audits the infrastructure-attrition denominator (57/103 → 265/448). Print-only. | FINAL |
| `build_table1.py` | Existing Table 1: launched + Phase III Class A targets with Open Targets example agents. | FINAL as a script; its drug and indication columns are unreliable (D2, D3) |
| `make_table2.py` | Existing Table 2: method availability, rendered PDF/SVG/PNG + CSV. | FINAL |
| `make_figures.py` | Existing Figures 1–4. | FINAL |
| `make_supplementary.py` | Assembles `deliverables/supplementary/` with dictionary, README, manifest. | FINAL |
| `partb_classA.py` | Genome-wide Class A enrichment tests (HPA tissue-restricted / secreted, GTEx TPM). | FINAL (feeds atlas) |
| `partc_concordance.py` | Part C: cis-MR direction vs approved-drug mechanism across 1,310 pairs. | INTERMEDIATE (null result; supplies attrition counts) |
| `part3_headtohead.py` | **Class C refutation test** — paired coloc.abf PP4, GTEx blood vs eQTLGen. | FINAL (documented negative result) |
| `part3_analyse.py` | Statistics for R14.3–R14.5 on Part 3 outputs. | INTERMEDIATE |
| `part3_degradation.py` | LD-dependent Part 3 arms under rule R14.1. | STALE — rule R14.1 flawed, outputs withdrawn |
| `fix2_r15.py` | Re-specified locus rule R15 (outcome must be genome-wide significant). | FINAL (method-availability numbers) |
| `part2_ukb_ld.py` | Single-locus test of UK Biobank LD at HMGCR. | INTERMEDIATE (Gate D diagnosis) |
| `coloc_diagnose.py`, `coloc_signals.py` | Gate D diagnosis at HMGCR–LDL-C. | INTERMEDIATE (Gate D diagnosis) |
| `susie_coloc.py`, `run_susie_real.py` | SuSiE-RSS + coloc.susie implementation and real-locus run. | INTERMEDIATE (Gate D; no output may be used as evidence) |
| `test_coloc.py`, `test_susie_coloc.py`, `test_coloc_assertion.py`, `test_harmonize.py` | Synthetic validation suites. | FINAL as tests |
| `control_chain_v2.py` | Positive controls on eQTLGen (HMGCR PASS, IL6R FAIL). | FINAL |
| `control_chain.py` | Positive controls on GTEx v10 liver (all three fail). | STALE as a pipeline; kept as the record of why the resource changed |
| `harmonize.py` | Allele harmonization module (R6). | FINAL |
| `tabix_http.py` | Pure-Python tabix over HTTP range requests. | FINAL (post-fix) |
| `audit_prefix.py` | Recomputes numbers produced before the tabix upper-bound fix. | INTERMEDIATE (audit record) |
| `d1b_scan_eqtlgen.py`, `d1b_tested_universe.py`, `d1b_gene_coords.py`, `d1b_hits_by_ensg.py` | eQTLGen scan, tested-gene universe, GRCh37 coordinates, ENSG-join re-extraction. | FINAL (instrument pipeline) |
| `d1b_clump.py`, `d1b_clump_delta.py` | LD clumping per target; reconciliation onto the ENSG join. | FINAL (instrument pipeline) |
| `d1b_census.py` | Genome-wide state census (ABSENT / TESTED-NULL / INSTRUMENTED). | FINAL |
| `tissue_power.py` | Instrument yield per GTEx tissue (Lane B upper bounds). | FINAL |
| `select_gwas.py` | Mechanical outcome-GWAS selection (R2). | FINAL |
| `pp_duplicates.py` | Duplicate analysis of `pp.tsv`. | FINAL (supports Decision 1) |
| `check_tbi.py`, `dirsize.py` | `.tbi` HEAD probes; FTP directory sizing. | INTERMEDIATE (Deliverable 0 utilities) |

## 3. `deliverables/cv/` — the atlas (12)

| file | what it is | status |
|---|---|---|
| `cardiovascular_atlas.csv` / `.parquet` | **The resource**: 630 rows × 25 columns, one per CV target. | FINAL — row-by-row recomputation matches (verify [V03]) |
| `CARDIOVASCULAR_ATLAS.md` | Narrative of the atlas results. | FINAL with errors: IL6R "lead p" mislabelled (D9), launched table lists 19 of 21 (D10), coloc 96.9 % framing (D4) |
| `ANALYSES_1_2_3.md` | Class B bound, CV-definition robustness, attrition audit. | FINAL |
| `LIMITATIONS.md` | Existing 12-section limitations document. | FINAL, with one stale sentence (D8) |
| `cv_summary.json` | Set sizes, state counts, enrichment 2×2. | FINAL |
| `cv_report.json` | Class counts, modality, CV Part C and coloc availability. | FINAL, but carries the superseded 57/103 sub-rate without labelling it |
| `cv_robustness.json` | Both CV definitions, all headline metrics. | FINAL |
| `classB_bound.json`, `classB_upper_bound.csv` | Class B bound: 196 flagged rows. | FINAL |
| `table1_classA_launched_phase3.csv` | 21 launched + 18 Phase III Class A targets. | FINAL counts; drug/indication columns need rework (D2, D3) |
| `table2_method_availability.csv` | Method-availability table data. | FINAL |

## 4. `deliverables/supplementary/` (6)

| file | what it is | status |
|---|---|---|
| `cardiovascular_atlas.csv` / `.parquet`, `classB_upper_bound.csv` | Byte-identical copies of the `cv/` files (sha256 match, verify [V03]). | FINAL |
| `DATA_DICTIONARY.md` | Every column: type, values, source, derivation, licences. | FINAL (licences self-marked UNVERIFIED) |
| `README.md` | Bundle description. Says "Exactly two commands" then lists three. | FINAL (minor, D24) |
| `MANIFEST.json` | sha256 + byte counts for the bundle. | FINAL |

## 5. `deliverables/figures/` (15)

| file | what it is | status |
|---|---|---|
| `figure1_untestability_cascade.{pdf,svg,png}` | 630 → ABSENT 177 / TESTED-NULL 55 / MEASURABLE 386; "21 with an approved medicine". | FINAL; not the participant-flow diagram Phase 2 asks for |
| `figure2_tissue_yield.{pdf,svg,png}` | Instrument yield: blood 61.3 % vs heart LV 7.1 % etc. | FINAL |
| `figure3_closed_loop.{pdf,svg,png}` | Required modality for the 232 vs what that modality delivers. | FINAL; inherits the modality-heuristic problem (D11) |
| `figure4_classA_enrichment_forest.{pdf,svg,png}` | OR forest: CV atlas 10.0, strict 10.6, all targets 6.9. | FINAL; "untestable" here means ABSENT only, not Class A (D12) |
| `table2_method_availability.{pdf,svg,png}` | Rendered Table 2. | FINAL |

None of the four requested Phase 2 figures (flow, class by phase, launched/Phase III highlight,
Class C) exists yet.

## 6. Deliverables d0 – d1j (analysis history)

| path | what it is | status |
|---|---|---|
| `d0/DELIVERABLE_0.md` | Access verification: GWAS Catalog, eQTL Catalogue, 1000G, `pp.tsv` checksum. | INTERMEDIATE (methods provenance) |
| `d0/pp_duplicates.txt`, `pp_first20.txt` | Duplicate report and head of `pp.tsv`. | INTERMEDIATE |
| `d0/tbi_probe_ldl.csv`, `tbi_probe_random200.csv` | `.tbi` availability probes (74.49 % LDL census; 86.50 % random 200). | INTERMEDIATE |
| `d0/extracts/` (4) | HMGCR window extracts + provenance (GTEx liver, LDL GWAS). | INTERMEDIATE |
| `d1/PREREGISTRATION.md` | Pre-registration R1–R15 with amendments. | FINAL (methods source of record) |
| `d1/DELIVERABLE_1.md`, `d1_summary.json`, `egene_status.json` | GTEx v10 liver controls — all three fail. | INTERMEDIATE (justifies Decision 7) |
| `d1/gwas_selection_EFO_0004458.csv`, `_0004611.csv` | Outcome GWAS ranking for CRP and LDL-C. | INTERMEDIATE |
| `d1/extracts/` (6), `d1/harmonized/` (3), `d1/provenance/` (3) | GTEx-liver control extracts, merges, provenance. | STALE as results (resource replaced); kept as record |
| `d1b/DELIVERABLE_1B.md`, `census.json`, `census_stdout.txt` | Genome-wide instrument census: 652 / 229 / 1,576 / 60. | FINAL |
| `d1b/clumped_instruments_final.tsv` | 1,576 instrumented targets with independent SNPs. | FINAL (atlas input) |
| `d1b/clumped_instruments.tsv` | Clump output before the ENSG-join reconciliation. | INTERMEDIATE |
| `d1b/clumped_instruments_BUGGY.tsv`, `clumped_instruments_final_BUGGY.tsv`, `clump_ckpt_BUGGY_truncated/` (22) | Pre-fix outputs from the tabix truncation bug. | **STALE** — named BUGGY, do not use |
| `d1b/clump_ckpt/` (22) | Per-chromosome clump checkpoints, post-fix. | INTERMEDIATE |
| `d1b/eqtlgen_hits_by_ensg.tsv`, `eqtlgen_hits_targets.tsv` | eQTLGen significant hits per target (ENSG join / symbol join). | INTERMEDIATE (ENSG join is primary) |
| `d1b/eqtlgen_tested_universe.tsv`, `eqtlgen_sig_gene_universe.tsv` | Genes tested / genes with a significant hit. | FINAL (defines ABSENT) |
| `d1b/gene_coords_grch37.tsv`, `gene_coords_unmapped.tsv` | GRCh37 coordinates; 60 unmapped with reasons. | FINAL (see D13 on ACE, ITGB3) |
| `d1b/mapping_methodB_absent.txt`, `mapping_methodB_tested_no_instrument.txt` | Per-target ABSENT and TESTED-NULL lists, primary join. | FINAL |
| `d1b/mapping_absent_from_eqtlgen.txt`, `mapping_tested_no_instrument.txt` | Same lists under the symbol join (method A). | STALE — superseded by method B |
| `d1b/clump_bytes.json` | Bytes over the wire for clumping. | INTERMEDIATE |
| `d1c/DELIVERABLE_1C.md`, `d1c_summary.json` | **Positive controls on eQTLGen**: HMGCR PASS, IL6R FAIL. | FINAL |
| `d1c/HMGCR_LDL_harmonized_merge.tsv`, `IL6R_CRP_harmonized_merge.tsv`, `*_provenance.json`, `prefix_audit.json` | Per-variant merges and provenance for the controls. | FINAL |
| `d1d/PART_A_COLOC_DIAGNOSIS.md`, `coloc_diagnosis.json`, `hmgcr_gtexblood_merge.tsv` | Gate D diagnosis: ABF peakiness at HMGCR. | INTERMEDIATE (Gate D record) |
| `d1e/PART_B_CLASS_A.md`, `classA_tests.json`, `classA_annotation.tsv` | Genome-wide Class A enrichment; per-target state + HPA + TPM. | FINAL (atlas input) |
| `d1f/PART_C_CONCORDANCE.md`, `partc_pairs.tsv`, `partc_results.tsv`, `partc_discordant.tsv`, `target_direction.json` | Part C concordance — null result. | INTERMEDIATE |
| `d1f/partc_skipped_accessions.tsv` | 70 unindexed accessions, 578 pairs lost. | FINAL (attrition source) |
| `d1f/ckpt/` (85) | Part C per-accession checkpoints. | INTERMEDIATE |
| `d1g/PART_A_SUSIE_BLOCKER.md`, `susie_controls.json` | coloc.susie fails at HMGCR (PP4 0.148). | INTERMEDIATE (Gate D record) |
| `d1h/PART_3_DEGRADATION.md` | **Class C refutation write-up.** | FINAL for R14.4; its R14.5 "2 of 28 (7 %)" and RMSE rows were later WITHDRAWN (D7) |
| `d1h/part3_headtohead.tsv`, `ckpt_h2h/` (60) | Class C paired data: 60 attempted, 50 paired. | FINAL |
| `d1h/part3_results.tsv`, `ckpt/` (43), `loci.tsv` | R14.1 LD-dependent arm: 28 of 60 loci completed. | **STALE** — flawed selection rule, superseded by R15 |
| `d1h_ftp_qts.txt` | eQTL Catalogue FTP study listing (QTS000001–042). | INTERMEDIATE |
| `d1i/PART_2_UKB_LD.md`, `part2_ukb_ld.json` | UKB LD at HMGCR; RMSE 7.80 → 7.12. | INTERMEDIATE (Gate D record) |
| `d1j/FIXES_1_2_3.md`, `r15_results.tsv`, `ckpt/` (106) | Artefact assertion, R15 re-run (6 of 106 loci qualify). | FINAL (method availability) |
| `lit/PART_1_LITERATURE.md` | Literature check on coloc failure mechanisms, dated 2026-08-30. | INTERMEDIATE — will need refreshing in Phase 2 step 11 |

## 7. `data/` and `logs/`

| path | what it is | status |
|---|---|---|
| `data/raw/pp.tsv` | Minikel target–indication table, 29,476 rows. | FINAL input |
| `data/derived/pp_dedup_frozen.tsv` | **Frozen analysis set**, 29,456 pairs; sha256 matches `DECISIONS.md`. | FINAL |
| `data/raw/minikel_areas.tsv`, `minikel_indic.tsv`, `minikel_indic_topl_match.tsv`, `minikel_mesh_2023_topl_maps.tsv` | Minikel indication → therapeutic-area maps. | FINAL input |
| `data/raw/2019-12-11-cis-eQTLsFDR-ProbeLevel-…txt.gz` (3.88 GB), `…FDR0.05…txt.gz`, `eqtlgen_SNP_AF.txt.gz` | eQTLGen phase I full / significant cis-eQTL files and allele frequencies. | FINAL input |
| `data/raw/Homo_sapiens.GRCh37.87.gtf.gz`, `hg19ToHg38.over.chain.gz` | Gene coordinates; liftover chain. | FINAL input |
| `data/raw/QTD000131/136/141/251/256/266.permuted.tsv.gz` | eQTL Catalogue R8 GTEx v10 permuted files (aorta, coronary, tibial, atrial, LV, liver). | FINAL input (tissue yield) |
| `data/raw/gtex_v10_median_tpm.gct.gz`, `proteinatlas.tsv.zip` | GTEx median TPM; Human Protein Atlas. | FINAL input (retrieval URL not logged) |
| `data/raw/ot2606/` (3 parquet) | Open Targets 26.06 drug molecule + mechanism of action. | FINAL input (release manifest reports `failure`, Decisions open item 4) |
| `data/raw/gwas-catalog-associations-full.zip`, `gwas-catalog-download-studies-v1.0.3.1.txt`, `harmonised_list.txt` | GWAS Catalog associations, studies, harmonised file list. | FINAL input |
| `data/raw/dataset_metadata_r8.tsv`, `dataset_metadata_upcoming.tsv` | eQTL Catalogue dataset metadata. | FINAL input |
| `data/raw/20130606_g1k_3202_samples_ped_population.txt`, `integrated_call_samples_v3.20130502.ALL.panel` | 1000 Genomes sample panels. | FINAL input |
| `data/raw/ukbld/` (7, ≈3.7 GB) | UKB LD block chr5:74–77 Mb and extracted arrays. | INTERMEDIATE (Gate D only) |
| `data/derived/ldl_*`, `harmonised_accessions.txt` (6) | Deliverable 0 LDL study census working files. | INTERMEDIATE |
| `logs/retrieval_log.csv` | URL, UTC time, bytes, sha256 per raw file. | FINAL, incomplete (GTEx, HPA, GWAS Catalog associations, Open Targets not logged) |
| `logs/clump_run.log`, `partc.log`, `part3.log`, `pysam_install_attempt.log`, `minikel_tree.json` | Run logs; Minikel GitHub tree listing. | INTERMEDIATE |

## 8. Outside the project root

| path | what it is | status |
|---|---|---|
| `C:\Users\Deven\causal2drug_feasibility\` | Predecessor feasibility stage (FEAS001–004, `AUDIT_02.md`, `FEASIBILITY_VERDICT.md`, ~39 `.tbi` indexes). | **Not read in this pass.** Dated 08-24/25, before every frozen decision; treat as STALE unless you say otherwise |

## 9. Created in this pass (`manuscript_prep/`)

| file | what it is |
|---|---|
| `verify/verify_numbers.py` | Read-only recomputation of every number in `NUMBERS.md`. Imports nothing from `src/`. |
| `verify/verify_output.txt` | Its captured output, blocks tagged [V01]–[V11]. |
| `verify/pcsk9_scan.txt` | Direct scan of the full 127,341,798-row eQTLGen file. |
| `verify/class_by_phase.csv`, `class_by_indication_mesh.csv`, `headline_targets.csv` | Tables behind `NUMBERS.md` §3–§4 and `HEADLINE.md`. |

Added 2026-10-10 after the Phase 1 review (see `DECISIONS_LOG.md` for what each answers):

| file | what it is |
|---|---|
| `DECISIONS_LOG.md` | The eight post-review decisions and their status; the Class B retirement. |
| `CV_SET_RULE.md` | Proposed rule for the extended CV set. **Awaiting approval.** |
| `DEFINITIONS.md`, `tested_null_detail.csv` | What "tested-null" means, with code; the single "launched" definition. **Awaiting review.** |
| `TABLE3_GROUPING_PROPOSAL.md`, `table3_grouping_proposal.csv` | Every MeSH term with pair counts and proposed indication areas. **Awaiting approval.** |
| `UNMAPPED_RESOLUTION.md`, `unmapped_resolution.csv` | Gene-record fixes for the 12 unmapped targets. |
| `COLOC_CLAIM_CORRECTION.md` | Correction of the "96.9 % of cardiovascular loci" figure. |
| `DEVIATION_CLASS_C_FIGURE.md` | Dated deviation note for the Class C supplementary figure. |
| `drug_links/` | ChEMBL 37 target → drug → indication tables and README. |
| `verify/list_cv_candidates.py`, `table3_grouping_proposal.py`, `describe_tested_null.py`, `launched_reconciliation.py`, `resolve_unmapped.py`, `gtf_lookup_unmapped.py`, `coloc_claim_check.py`, `build_chembl_links.py`, `check_excluded_cardiac_indications.py` | The read-only scripts behind those files; each has a saved `*_output.txt`. |
| `verify/cache/` | Cached HGNC, ChEMBL and MeSH API responses, so the scripts re-run offline. |

Added 2026-10-10 in Phase 2 (see `README.md` for the index):

| file | what it is |
|---|---|
| `README.md`, `KEY_NUMBERS.md` | Index; corrected numbers on the approved definitions. |
| `analysis/classify.py`, `analysis/out/` | Corrected classification on the primary and sensitivity sets; per-target tables; `key_numbers.json`. |
| `tables/` | Tables 1-3 and S1-S10 (CSV, Markdown, workbook) and `make_tables.py`. |
| `figures/` | Figures 1-3 and S1 (PDF, SVG, PNG) and `make_figures.py`. |
| `METHODS_FACTS.md`, `OUTLINE.md`, `LIMITATIONS.md`, `JOURNAL_FIT.md`, `PRIOR_ART_CHECK.md` | Phase 2 documents. |
| `originals_before_correction/` | Copies of the three existing Markdown files as they were before the "20" correction. |

Existing project files changed in Phase 2: `deliverables/cv/ANALYSES_1_2_3.md`,
`CARDIOVASCULAR_ATLAS.md` and `LIMITATIONS.md` (the "20" correction only).
