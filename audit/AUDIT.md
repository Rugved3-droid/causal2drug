# AUDIT — second-implementation reproduction of the CV set, eQTLGen status, recovery map and two controls

**Who ran it:** a second implementation, written in a separate AI-assisted session from the rule
documents alone, without access to the original code or per-target outputs. Tier and category
names below were updated on 2026-10-11 to the current ones; no number was changed.

Audit date 2026-10-11. All code is in `audit/`; nothing in the original project was edited.

## Summary

1. **CV set, eQTLGen status: exact match, 726 of 726 targets.** 726 targets; no eQTL data 217 (188 autosomal / 27 X-linked / 2 mitochondrial); no qualifying instrument 63; unresolved 3. Rebuilt from the raw eQTLGen file, GTF and HGNC, not from the project's intermediate files.
2. **Recovery map tiers: exact match, 280 of 280 targets**, after I fixed two errors in my own code. Candidate: CV-relevant tissue or plasma protein 151, candidate: other tissue only 72, locked never assayed 35, locked assayed 22. All GTEx cells match on every target.
3. **HMGCR → LDL-C: match.** β +0.3891 (SE 0.0120), p 4.2 × 10⁻²²⁹, F 224.0.
4. **PCSK9 → LDL-C: MISMATCH.** Original β +0.376 (SE 0.0047, F 1,723); mine +0.387 (SE 0.0048, F 2,305). Both pass. **Cause: an error in the original.**
5. The error: `build_recovery_map.py:370` and `:417` turn every p below 10⁻³⁰⁰ into 0.0, so "lowest p" becomes a tie broken by sheet or row order. For 30 UKB-PPP targets the recorded sheet is the discovery sheet (ST9) although the combined sheet (ST10) has the lower p; for 9 of them and 2 deCODE targets the recorded lead variant is not the lowest-p one. No category or tier changes.
6. **A second, smaller original error, no effect on any tier:** Fenland records CCN1 and GRK2 as "not assayed" (they are in the table under CYR61 and ADRBK1).
7. **One rule ambiguity, no effect on any tier:** whether a deCODE aptamer marked "Included in analysis = No" counts as assayed (ALPP, CDK5, FGA, ITGA2B).
8. **My two errors (fixed, first-run output kept):** chromosome X coded "23" vs "X" in UKB-PPP (6 targets); Fenland gene body looked up by symbol only (VEGFD).
9. **Limit of this audit:** the rule documents quote the totals (726, 217, 63, 280, 57 locked), so I was blind to per-target results and to the control estimates, not to the headline counts.

## Step 1 — rule documents read, and inputs

Read before any code was written: `PROTOCOL.md`, `manuscript_prep/CV_SET_RULE.md`, `DEFINITIONS.md`,
`DECISIONS_LOG.md`, `manuscript/RECOVERY_RULE.md` (with clarifications 1–7 and amendments 8–9),
`manuscript/RESCUE_CONTROLS_PROTOCOL.md`, `DEVIATION_CLASS_C_FIGURE.md`. Three more were needed
because the listed documents defer to them for rules: `DECISIONS.md` (dedup, Decision 1),
`deliverables/d1/PREREGISTRATION.md` (R2, R2a, R6, R12.2) and `UNMAPPED_RESOLUTION.md` (HGNC rule);
`manuscript/DATA_ACCESS.md` for file locations.

Not opened before Step 3: `src/`, `manuscript_prep/analysis/`, `recovery_scripts/`, any results CSV,
`KEY_NUMBERS.md`, the manuscript.

| raw input | path | exists |
|---|---|---|
| Minikel pairs | `data/raw/pp.tsv` (29,476 rows) | yes |
| Minikel indications with `areas` | `data/raw/minikel_indic.tsv` (1,014) | yes |
| Minikel single-branch file (sensitivity sets only; not used) | `data/raw/minikel_indic_topl_match.tsv` | yes |
| Ensembl GRCh37 release 87 GTF | `data/raw/Homo_sapiens.GRCh37.87.gtf.gz` | yes |
| eQTLGen full cis-eQTL file | `data/raw/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | yes |
| eQTLGen allele frequencies | `data/raw/eqtlgen_SNP_AF.txt.gz` | yes |
| GWAS Catalog studies v1.0.3.1; harmonised list | `data/raw/gwas-catalog-download-studies-v1.0.3.1.txt`; `harmonised_list.txt` | yes |
| hg19→hg38 chain | `data/raw/hg19ToHg38.over.chain.gz` | yes |
| GTEx v8 eQTL archive (49 tissues) | `manuscript/recovery_data/GTEx_Analysis_v8_eQTL.tar` | yes |
| UKB-PPP, deCODE, Fenland workbooks | `recovery_data/41586_2023_6592_MOESM3_ESM.xlsx`, `ferkingstad2021_MOESM4.xlsx`, `pietzner2021_supplementary_tables.xlsx` | yes |
| LDL-C outcome GWAS GCST90239658 (remote, tabix) | EBI FTP, harmonised file and `.tbi` | reachable 2026-10-11 |
| HGNC | project cache `verify/cache/` exists; I downloaded the complete set afresh (`audit/out/hgnc_complete_set.txt`) | yes |

The rule documents also name four **derived** files (`d1b/gene_coords_grch37.tsv`,
`eqtlgen_tested_universe.tsv`, `eqtlgen_hits_by_ensg.tsv`, `analysis/out/per_target_primary.csv`).
They exist. I did not use them; I rebuilt their content from the raw files.

## Step 2 — my reproduction, with commands

Run from `D:\causal2drug`. Tags are printed by the scripts; logs are in `audit/out/aNN_output.txt`.

| # | number | value | command | tag |
|---|---|---|---|---|
| a | pairs after dedup and 4 dropped programmes | 29,456 (2,517 genes, 1,014 indications) | `python audit/a01_cv_set_and_coords.py` | A01 |
| a | CV terms / pairs / **targets**, primary set | 93 / 2,128 / **726** | same | A02 |
| a | symbols without a unique GTF record / fixed by HGNC / unresolved | 15 / 12 / 3 (MIR195, MIR29B1, MIR92A1) | same | A06 |
| b | eQTLGen rows / genes tested | 127,341,798 / 19,250 | `bash audit/a00_eqtlgen_scan.sh` then `python audit/a04_classify.py` | B01 |
| b | instrumented | 443 | `python audit/a04_classify.py` | B03 |
| b | **no eQTL data** | **217 = 188 autosomal + 27 X-linked + 2 mitochondrial** | same | B04 |
| b | **no qualifying instrument**; upper bound | **63**; 280 | same | B05 |
| b | autosomal-only lower bound | 188 of 697 (27.0 %) | same | B07 |
| b | launched (66) | 39 instrumented, 23 no eQTL data, 4 no qualifying instrument | same | B08 |
| c | pQTL sheets dumped; GTEx rows extracted | — | `python audit/a02_dump_pqtl_sheets.py`; `python audit/a03_gtex_extract.py` | — |
| c | in scope | 280 | `python audit/a05_recovery_map.py` | C01 |
| c | categories | candidate-tissue 128, candidate-protein 24, candidate-both 71, locked 57 | same | C05 |
| c | **tiers** | **candidate: CV-relevant tissue or plasma protein 151; candidate: other tissue only 72; locked never assayed 35; locked assayed, no instrument 22** | same | C06 |
| c | CV-tissue instrument; plasma cis-pQTL instrument | 91; 95 (UKB-PPP 76, deCODE 68, Fenland 58) | same | C07 |
| d | outcome study by R2 + R2a | GCST90239658, N 1,320,016, European only | `python audit/a06_controls.py` | D01 |
| d | **PCSK9 → LDL-C** | UKB-PPP ST10 rs11591147, β_exp −1.12113 (SE 0.02335), F 2,304.8; **Wald β +0.3868 (SE 0.0048; 95 % CI 0.3774–0.3963)**, z ≈ 80, p underflows; aligned; PASS | same | D02, D03 |
| d | **HMGCR → LDL-C** | eQTLGen rs6453133, z 14.9667, N 31,684, EAF 0.313; **Wald β +0.3891 (SE 0.0120; 95 % CI 0.3655–0.4127)**, p 4.24 × 10⁻²²⁹, F 224.0; aligned; PASS | same | D05, D06 |

The tier and category values in the table are from the corrected run. The first run, before I had
opened any original output, gave candidate-tissue 130 / both 69 and candidate: CV-relevant tissue or plasma protein 150 / candidate: other tissue only 73
(`audit/out/a05_output_FIRST_RUN_before_my_bugfix.txt`); see mismatches M5 and M6.

## Step 3 — comparison with `KEY_NUMBERS.md` and `manuscript/recovery_map.csv`

`python audit/a07_compare.py` — per-target differences are written to
`audit/out/compare_status.csv` and `audit/out/compare_recovery.csv`. Per-target status was compared
against `analysis/out/per_target_primary.csv`, the control against
`manuscript/rescue_controls/rescue_controls_results.csv`.

### Matches

| item | compared | result |
|---|---|---|
| Target set (KEY_NUMBERS §1) | 726 symbols | identical; 93 terms and 2,128 pairs also match |
| Status, sub-type, Ensembl ID, chromosome, phase | 726 targets × 5 fields | 0 differences |
| Number of qualifying eQTLGen variants; lead variant | 726 targets | 0 differences |
| KEY_NUMBERS §2 primary row | 443 / 217 / 188 / 27 + 2 / 188 of 697 / 63 / 280 / 3 | all match |
| KEY_NUMBERS §3 by phase | 5 rows × 8 columns | all match (`a08_extra_checks.py`, B09) |
| KEY_NUMBERS §9 | 3 unresolved in primary | match |
| Recovery scope | 280 symbols | identical |
| Category and reporting tier | 280 targets | 0 differences (corrected run) |
| GTEx: ID used (incl. amendment 8 for ADRA2B, GDF2, MIR145), tissues tested, non-blood tissues with an instrument, CV tissues, whole-blood flag, best CV tissue and variant | 280 targets × 7 fields | 0 differences |
| UKB-PPP assayed; UKB-PPP, deCODE, Fenland instrument yes/no | 280 targets | 0 differences (corrected run) |
| Fenland lead variant | 58 targets | 0 differences |
| HMGCR → LDL-C (KEY_NUMBERS §8) | β, SE, p, F | match to the printed precision |
| PCSK9 outcome side | study, harmonization, β_out −0.433696, SE 0.005426 | match |

### Mismatches

| # | what | targets | cause |
|---|---|---|---|
| **M1** | UKB-PPP sheet recorded for the lead: original ST9, mine ST10 | 30: ACHE, AGT, ALPP, ANGPTL3, ANGPTL4, AOC3, APOE, CDH17, CPB2, F11, F9, GHR, GLA, IL12B, IL5RA, ITGB6, KDR, LAMP2, LPA, MBL2, MMP10, NELL1, PCSK9, PDGFRB, PLG, PON1, SCG3, SERPINA4, TYRO3, VEGFD | **Error in the original.** `build_recovery_map.py:370` computes `10 ** -lp if lp < 300 else 0.0`. Any p below 10⁻³⁰⁰ becomes 0.0, the ST9 and ST10 rows tie, and the sort keeps ST9 because it is read first. Clarification 6 says the lowest p is recorded; in all 30 the ST10 log10(p) is larger (PCSK9: 502.3 vs 375.8). |
| **M2** | UKB-PPP lead variant itself differs | 9 of the 30: ANGPTL3, CDH17, F11, F9, GHR, GLA, ITGB6, SCG3, TYRO3 | Same error. Example GHR: original rs4073476 (ST9, 10⁻¹³⁹⁸); lowest p is rs6873545 (ST10, 10⁻²¹³⁷). |
| **M3** | deCODE lead variant differs | MBL2 (original rs1800451, 10⁻⁴⁷⁷; lowest rs7899547, 10⁻³⁰¹⁸); SERPINA4 (original rs5511, 10⁻⁴⁶⁶; lowest rs5512, 10⁻⁵⁰⁵) | Same error, `build_recovery_map.py:417`. |
| **M4** | **PCSK9 → LDL-C estimate**: original +0.3761 (SE 0.0047, F 1,722.8); mine +0.3868 (SE 0.0048, F 2,304.8) | PCSK9 | **Consequence of M1.** The protocol takes "exactly the lead qualifying variant recorded in `recovery_map.csv`" with the sheet in `ukbppp_note`; that note says ST9 because of the tie. Same variant (rs11591147), same outcome row; dividing my outcome β by the ST9 exposure β gives 0.3761 exactly (tag E06). Verdict PASS either way. |
| M5 | UKB-PPP instrument: mine no, original yes (first run) | ACE2, CD40LG, F9, GLA, LAMP2, VEGFD; category changed for ACE2 and LAMP2, tier for LAMP2 | **Error in my code.** ST9/ST10 code chromosome X as "23", ST3 as "X"; I compared the strings. Fixed in `a05_recovery_map.py`. |
| M6 | Fenland instrument and assayed: mine no / UNVERIFIED, original yes (first run) | VEGFD | **Error in my code.** I looked the ST3 gene body up by symbol; ST3 lists VEGFD as FIGF and the row matches by Ensembl ID. Fixed. |
| M7 | Fenland assayed: original no, mine yes | CCN1, GRK2 | **Error in the original** (small). The original collects Ensembl-ID matches only from cis-labelled ST2 rows (`build_recovery_map.py`, loop after `[R09]`). Both proteins are in ST3 as CYR61 and ADRBK1 and in ST2 with the target's Ensembl ID on trans rows. Clarification 5 allows an Ensembl-ID match. No tier changes: both are assayed by deCODE. |
| M8 | deCODE assayed: original no, mine yes | ALPP, CDK5, FGA, ITGA2B | **Rule ambiguity.** All four are in ST01 with "Included in analysis = No". The rule says "assayed or not assayed (where the table or its companion protein list allows this to be read)" and does not say which column decides. The original's reading (excluded aptamers are not assayed) is the more defensible one. No tier changes under either reading (tag C08). |

### What M1–M4 change downstream

- Categories, tiers and every count in the recovery map: nothing.
- `recovery_map.csv` columns `ukbppp_lead_variant`, `ukbppp_p`, `ukbppp_note` for 30 targets and `decode_lead_variant` for 2: the recorded lead is not the lowest-p qualifying variant.
- PCSK9 control: β +0.376 becomes +0.387 if the combined-cohort row is used. Still a pass.
- ANGPTL3 control: the lead changes from rs763015008 to rs775677524. Both are indels, so "not harmonizable" stands. I did not run it.
- APOC3, REN, AGTR1: REN's lead and AGTR1's tissue and variant are the same in my map. APOC3 is not assayed by UKB-PPP; its deCODE and Fenland leads (rs964184) match.
- If the owner intended the discovery cohort to be the exposure, that is a legitimate choice, but it is not what the rule text says and it was not made on purpose: the sheet was picked by a tie.

## Rule gaps I had to fill

None of these produced a mismatch.

1. **Z to β for eQTLGen.** No rule document states the conversion. I used β = z / √(2f(1−f)(N+z²)), f from eQTLGen's allele-frequency file. It reproduces the original to four digits.
2. **HMGCR cis-window build.** R1 defines the window on GRCh38 from Ensembl REST; the status rule uses the GRCh37 GTF. I used GRCh37 ± 100 kb and lifted the lead. Same lead (rs6453133).
3. **Amendment 8, GDF2 and MIR145.** The amendment notes a strand difference and a gene-model difference but says they do not gate the match. I applied the match to all three; GDF2 and MIR145 stay locked because they have no instrument. Same as the original.
4. **Fenland build check.** I count 4,038 of 4,055 aptamer–gene pairs identical to GRCh37.87; `DATA_ACCESS.md` says 3,917 of 3,933 single-gene aptamers. Different counting unit, same conclusion.
5. **HGNC snapshot.** Mine is today's complete set, not the project's cache. All 726 Ensembl IDs agree.

## Extra check

`python audit/a08_extra_checks.py` [F01]: no target in "locked: never assayed" appears in any of the
three protein lists under an HGNC previous symbol or by UniProt accession. The two hits (SLC6A2 via
"NAT1"/"NET1") are symbol collisions with other genes, not the same protein.

## Not audited

Sensitivity sets S1–S3, the all-2,517 row, sub-types T1–T3, cis-window sensitivity, drug
verification (KEY_NUMBERS §4), enrichment odds ratios, the other three rescue controls, IL6R, and
anything in the manuscript text.
