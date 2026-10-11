# NUMBER_TRACE

Every number in `manuscript.md`, with the output file and line that carries it and the command that regenerates that file. Built by `python manuscript_prep/manuscript/draft/make_number_trace.py`, which searches each file for the stated value; a value it cannot find is written **UNVERIFIED**.

Paths are relative to `D:\causal2drug\manuscript_prep\`. Run commands from `D:\causal2drug`. Numbers marked in the manuscript depend on the pending manual review of `drug_links/MANUAL_REVIEW.csv`.

| # | where in the manuscript | number(s) | file : line | regenerate with |
|---|---|---|---|---|
| 1 | Abstract; Results 1 | 726 targets; 443 instrumented (61.0%); 217 no eQTL data (29.9%); 63 no qualifying instrument (8.7%); upper bound 280 (38.6%); 3 unresolved; autosomal-only 188 of 697 (27.0%) | `analysis/classify_output.txt` : 34 | `python manuscript_prep\analysis\classify.py` |
| 2 | Methods; Results 4 | 27 X-linked and 2 mitochondrial targets | `analysis/classify_output.txt` : 34 | `python manuscript_prep\analysis\classify.py` |
| 3 | Methods | 29,476 rows; 29,460 pairs; 29,456 pairs; 2,517 targets; 1,014 indications | `tables/table1_flow_and_cohort.md` : 8 | `python manuscript_prep\tables\make_tables.py` |
| 4 | Methods | 29,456 pairs after removing 4 programmes | `tables/table1_flow_and_cohort.md` : 7 | `python manuscript_prep\tables\make_tables.py` |
| 5 | Methods | 93 indications, 2,128 pairs, 726 targets; sensitivity sets 630, 663, 594 | `analysis/classify_output.txt` : 6 | `python manuscript_prep\analysis\classify.py` |
| 6 | Methods | sensitivity sets 630, 663 and 594 targets | `tables/table1_flow_and_cohort.md` : 12 | `python manuscript_prep\tables\make_tables.py` |
| 7 | Methods; Results 8 | 127,341,798 rows in the full eQTLGen file | `verify/pcsk9_scan.txt` : 2 | `see header of verify\pcsk9_scan.txt` |
| 8 | Results 1 | lower bound 27.6% to 29.9%, upper bound 36.9% to 38.6% across set definitions | `tables/table1_flow_and_cohort.md` : 14 | `python manuscript_prep\tables\make_tables.py` |
| 9 | Results 1 | upper bounds 38.6%, 37.6%, 37.9%, 36.9% | `tables/table1_flow_and_cohort.md` : 19 | `python manuscript_prep\tables\make_tables.py` |
| 10 | Results 1 | all 2,517 targets: 26.3% and 35.7% | `KEY_NUMBERS.md` : 62 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 11 | Results 1 | ±1 Mb window moves 1 target, IL1B | `analysis/classify_output.txt` : 156 | `python manuscript_prep\analysis\classify.py` |
| 12 | Results 2 | no eQTL data by phase: 30.0%, 28.6%, 33.3%, 18.7%, 34.8% | `tables/table2_status_by_phase.md` : 5 | `python manuscript_prep\tables\make_tables.py` |
| 13 | Results 2 | phase I 28.6% | `tables/table2_status_by_phase.md` : 6 | `python manuscript_prep\tables\make_tables.py` |
| 14 | Results 2 | phase II 33.3% | `tables/table2_status_by_phase.md` : 7 | `python manuscript_prep\tables\make_tables.py` |
| 15 | Results 2 | phase III 18.7% | `tables/table2_status_by_phase.md` : 8 | `python manuscript_prep\tables\make_tables.py` |
| 16 | Results 2 | launched 34.8% | `tables/table2_status_by_phase.md` : 9 | `python manuscript_prep\tables\make_tables.py` |
| 17 | Abstract; Results 2 | chi-square 6.33, df 4, P=0.18 | `KEY_NUMBERS.md` : 81 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 18 | Results 3 | 35 launched targets pass the ChEMBL strict rule | `tables/table2_status_by_phase.md` : 12 | `python manuscript_prep\tables\make_tables.py` |
| 19 | Results 3 | 33 launched targets with a drug ever approved by FDA or EMA; 9 no eQTL data (27.3%); 3 no qualifying instrument; upper bound 12 (36.4%) | `tables/table2_status_by_phase.md` : 11 | `python manuscript_prep\tables\make_tables.py` |
| 20 | Results 3 | 33 of 35 have an FDA or EMA approval record | `manuscript/revision_sources/approved_35_output.txt` : 4 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 21 | Results 3 | MTOR and PLG have no approval record | `manuscript/revision_sources/approved_35_output.txt` : 5 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 22 | Results 3 | same count from the classifier | `analysis/classify_output.txt` : 78 | `python manuscript_prep\analysis\classify.py` |
| 23 | Results 3 | permissive rule: 20 and 4 of 58 | `tables/table2_status_by_phase.md` : 13 | `python manuscript_prep\tables\make_tables.py` |
| 24 | Results 3 | names of the 13 targets | `KEY_NUMBERS.md` : 136 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 25 | Results 3 | IL1B, ITGA2B, NPC1L1 tested without a qualifying instrument | `KEY_NUMBERS.md` : 137 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 26 | Abstract; Results 4 | 188 autosomal no-eQTL-data targets; 182 with a GTEx value; 180 below 1 TPM | `tables/tableS10_no_eqtl_data_by_chromosome.md` : 5 | `python manuscript_prep\tables\make_tables.py` |
| 27 | Results 4 | 12 of 27 X-linked targets expressed in blood at ≥1 TPM | `tables/tableS10_no_eqtl_data_by_chromosome.md` : 5 | `python manuscript_prep\tables\make_tables.py` |
| 28 | Results 4 | tissue restriction OR 7.47 (4.12–13.55) | `analysis/classify_output.txt` : 162 | `python manuscript_prep\analysis\classify.py` |
| 29 | Results 4 | autosomal only OR 12.18 (5.58–26.62) | `KEY_NUMBERS.md` : 206 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 30 | Results 4 | no qualifying instrument OR 1.39 (0.77–2.51) | `analysis/classify_output.txt` : 163 | `python manuscript_prep\analysis\classify.py` |
| 31 | Methods | GTEx v8: 49 tissues | `manuscript/recovery_map_summary.txt` : 6 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 32 | Methods | UKB-PPP threshold P<1.7×10⁻¹¹ (workbook sheet title) | `manuscript/DATA_ACCESS.md` : 48 | `Fenland workbook, sheet ST1, row 'Participants [N]'` |
| 33 | Methods | 3 targets matched after an identifier change; 1 changed category | `manuscript/recovery_map_summary.txt` : 260 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 34 | Abstract; Results 5 | 280 in scope; 223 fixable (79.6%); 57 locked | `manuscript/recovery_map_summary.txt` : 21 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 35 | Abstract; Results 5 | tiers: 151 (53.9%), 72 (25.7%), 35 never assayed, 22 assayed; locked 57 (20.4%) | `manuscript/recovery_map_summary.txt` : 150 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 36 | Results 5 | GTEx non-blood instrument for 199 targets | `manuscript/recovery_map_summary.txt` : 268 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 37 | Results 5 | plasma cis-pQTL instrument for 95 targets | `manuscript/recovery_map_summary.txt` : 283 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 38 | Results 5 | candidate: CV-relevant tissue or plasma protein: 56 tissue only, 60 protein only, 35 both | `manuscript/recovery_map_summary.txt` : 182 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 39 | Results 5 | tibial artery 44 | `manuscript/recovery_map_summary.txt` : 273 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 40 | Results 5 | aorta 34 | `manuscript/recovery_map_summary.txt` : 271 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 41 | Results 5 | atrial appendage 29 | `manuscript/recovery_map_summary.txt` : 270 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 42 | Results 5 | kidney cortex 0 | `manuscript/recovery_map_summary.txt` : 275 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 43 | Results 5 | off-target tissues: testis 21, tibial nerve 17, fibroblasts 16, subcutaneous adipose 14 | `manuscript/recovery_map_summary.txt` : 256 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 44 | Results 5 | 27 of 72 on a single tissue; 7 on testis alone | `manuscript/recovery_map_summary.txt` : 257 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 45 | Results 5 | never assayed 35 (list); 9 microRNA, 2 mitochondrial | `manuscript/recovery_map_summary.txt` : 258 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 46 | Results 5 | assayed, no instrument 22 (list) | `manuscript/recovery_map_summary.txt` : 259 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 47 | Results 5 | 25 of 35 never-assayed targets tested in ≥1 GTEx tissue (35 minus the 10 in no egenes file) | `manuscript/recovery_map_summary.txt` : 284 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 48 | Results 5 | no qualifying instrument: locked 14 of 63 (22.2%) | `manuscript/recovery_map_summary.txt` : 154 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 49 | Results 5 | autosomal no data: locked 37 of 188 (19.7%) | `manuscript/recovery_map_summary.txt` : 151 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 50 | Results 5 | X-linked: 23 of 27 with an instrument, 4 locked | `manuscript/recovery_map_summary.txt` : 152 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 51 | Results 6 | locked by phase: 27.0%, 16.7%, 15.2%, 24.0%, 3.7% | `manuscript/recovery_map_summary.txt` : 158 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 52 | Results 6 | phase I 16.7% | `manuscript/recovery_map_summary.txt` : 159 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 53 | Results 6 | phase II 15.2% | `manuscript/recovery_map_summary.txt` : 160 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 54 | Results 6 | phase III 24.0% | `manuscript/recovery_map_summary.txt` : 161 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 55 | Results 6 | launched 3.7% | `manuscript/recovery_map_summary.txt` : 162 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 56 | Abstract; Results 6; Discussion | heart failure / cardiomyopathy: 70 targets; 25 (35.7%), 24 (34.3%), locked 21 (30.0%), 16 never assayed | `manuscript/recovery_map_summary.txt` : 168 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 57 | Results 6 | thrombosis: 0 of 22 locked; 17 candidate: CV-relevant tissue or plasma protein | `manuscript/recovery_map_summary.txt` : 172 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 58 | Results 7 | 13 targets of Table 2: 8, 4 and 1 by tier (the text reports the 12 named, without IL1B: 8, 3, 1) | `manuscript/recovery_map_summary.txt` : 179 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 59 | Results 7 | tissues and pQTL status of each named target | `manuscript/display/table2_named_targets.md` : 5 | `python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py` |
| 60 | Results 7 | AVPR2: locked, never assayed; tested in 45 GTEx tissues | `manuscript/recovery_map.csv` : 227 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 61 | Results 7 | HCN4: locked, never assayed; tested in 25 GTEx tissues | `manuscript/recovery_map.csv` : 235 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 62 | Results 7 | AVPR1A, SLC6A2, PDE3A: candidate in a CV-relevant tissue or plasma protein (key 'candidate: CV-relevant tissue or plasma protein' in the CSV) | `manuscript/recovery_map.csv` : 121 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 63 | Results 8 | HMGCR β 0.389, P=4.2×10⁻²²⁹ | `KEY_NUMBERS.md` : 226 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 64 | Results 8 | IL6R β 0.376, P=1.2×10⁻⁵⁹, fails on direction | `KEY_NUMBERS.md` : 227 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 65 | Results 8 | PCSK9: 0 of 127,341,798 rows | `KEY_NUMBERS.md` : 228 | `python manuscript_prep\analysis\classify.py (summarised in KEY_NUMBERS.md)` |
| 66 | Limitations | 26 of 63 with a cis-eQTL at FDR<0.05 | `analysis/classify_output.txt` : 39 | `python manuscript_prep\analysis\classify.py` |
| 67 | Limitations | 46 of 93 indications in the 'other' area | `LIMITATIONS.md` : 91 | `python manuscript_prep\analysis\classify.py (summarised in LIMITATIONS.md)` |
| 68 | Introduction | Storm 2021: 2,786 of 4,863 | `manuscript/PRIOR_ART.md` : 19 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 69 | Introduction | Liang 2025: 2,888 of 4,479 | `manuscript/PRIOR_ART.md` : 17 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 70 | Introduction | Liu 2024: 2,525 of 4,479 | `manuscript/PRIOR_ART.md` : 18 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 71 | Introduction | Wang 2024: 2,531 of 4,317 | `manuscript/PRIOR_ART.md` : 20 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 72 | Introduction | Sadler 2023: ~12,550 of ~19,150; HCN4 | `manuscript/PRIOR_ART.md` : 22 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 73 | Introduction | Gaziano 2021: 1,016 of 1,263; 49 tissues | `manuscript/PRIOR_ART.md` : 21 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 74 | Introduction | eQTLGen: 19,250 autosomal genes expressed in blood and tested | `manuscript/PRIOR_ART.md` : 91 | `quotation from the cited paper's full text; see manuscript\prior_art_search\fulltext` |
| 75 | Supplement | Class C: 50 paired, 32 higher, 18 lower, medians 0.00821 and 0.01390, Wilcoxon P=0.0262, binomial P=0.0649 | `DEVIATION_CLASS_C_FIGURE.md` : 58 | `python manuscript_prep\verify\verify_numbers.py` |
| 76 | Supplement | Fenland build check: 3,917 of 3,933 identical to GRCh37.87; 11 to GENCODE v26 | `manuscript/recovery_map_summary.txt` : 13 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 77 | Supplement | X coverage: GTEx 1,293 X-linked genes | `manuscript/recovery_map_summary.txt` : 6 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 78 | Supplement | amendment 8 before: 222 fixable, 58 locked | `manuscript/recovery_out/recovery_map_summary_before_amendment8.txt` : 17 | `first run of build_recovery_map.py, kept as a record` |
| 79 | Abstract; Results 2 | no eQTL data, preclinical 30.0% (95% CI 25.3-35.2) | `analysis/phase_ci_output.txt` : 2 | `python manuscript_prep\analysis\phase_ci.py` |
| 80 | Results 2 | phase I 28.6% (20.3-38.6) | `analysis/phase_ci_output.txt` : 3 | `python manuscript_prep\analysis\phase_ci.py` |
| 81 | Results 2 | phase II 33.3% (26.7-40.7) | `analysis/phase_ci_output.txt` : 4 | `python manuscript_prep\analysis\phase_ci.py` |
| 82 | Results 2 | phase III 18.7% (11.5-28.9) | `analysis/phase_ci_output.txt` : 5 | `python manuscript_prep\analysis\phase_ci.py` |
| 83 | Results 2 | launched 34.8% (24.5-46.9) | `analysis/phase_ci_output.txt` : 6 | `python manuscript_prep\analysis\phase_ci.py` |
| 84 | Abstract; Results 2 | launched minus preclinical +4.8 points (95% CI -6.8 to +17.8) | `analysis/phase_ci_output.txt` : 7 | `python manuscript_prep\analysis\phase_ci.py` |
| 85 | Results 2 | upper bound difference +1.9 points (-10.4 to +15.0) | `analysis/phase_ci_output.txt` : 15 | `python manuscript_prep\analysis\phase_ci.py` |
| 86 | Abstract; Results 8 | PCSK9 candidate control: rs11591147, MR beta +0.387 (0.377 to 0.396); PASS | `manuscript/rescue_controls/run_log.txt` : 4 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 87 | Results 8 | PCSK9 control verdict | `manuscript/rescue_controls/run_log.txt` : 2 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 88 | Abstract; Results 8 | APOC3 candidate control: rs964184, MR beta +1.171 (1.153 to 1.190); PASS | `manuscript/rescue_controls/run_log.txt` : 9 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 89 | Results 8 | APOC3 control verdict | `manuscript/rescue_controls/run_log.txt` : 7 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 90 | Abstract; Results 8 | AGTR1 candidate control: MR beta -0.037 (-0.268 to 0.194), P=0.75; FAIL | `manuscript/rescue_controls/run_log.txt` : 12 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 91 | Results 8 | AGTR1 control verdict | `manuscript/rescue_controls/run_log.txt` : 10 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 92 | Abstract; Results 8 | ANGPTL3 control not harmonizable (indel) | `manuscript/rescue_controls/run_log.txt` : 5 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 93 | Results 8 | ANGPTL3 reason: not_snp | `manuscript/rescue_controls/run_log.txt` : 6 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 94 | Abstract; Results 8 | REN control not harmonizable (absent from outcome file) | `manuscript/rescue_controls/run_log.txt` : 13 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 95 | Results 8 | REN reason | `manuscript/rescue_controls/run_log.txt` : 14 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 96 | Results 8 | IL6R audit: MR variant rs2229238; eQTLGen lead rs10908839 | `manuscript/il6r_audit/il6r_audit_output.txt` : 1 | `python manuscript_prep\manuscript\il6r_audit\il6r_audit.py` |
| 97 | Results 8 | IL6R audit: r2 = 0.16 between rs2229238 and rs2228145 | `manuscript/il6r_audit/il6r_audit_output.txt` : 30 | `python manuscript_prep\manuscript\il6r_audit\il6r_audit.py` |
| 98 | Results 8 | IL6R audit: 358Ala allele, expression z -3.17 | `manuscript/il6r_audit/il6r_audit_output.txt` : 14 | `python manuscript_prep\manuscript\il6r_audit\il6r_audit.py` |
| 99 | Results 8 | IL6R audit: harmonization verified by hand | `manuscript/il6r_audit/il6r_audit_output.txt` : 53 | `python manuscript_prep\manuscript\il6r_audit\il6r_audit.py` |
| 100 | Limitations | GTEx v8 eQTL samples: kidney cortex 73 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 35 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 101 | Limitations | liver 208 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 37 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 102 | Limitations | coronary artery 213 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 6 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 103 | Limitations | adrenal gland 233 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 4 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 104 | Limitations | heart atrial appendage 372 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 33 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 105 | Limitations | heart left ventricle 386 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 34 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 106 | Limitations | aorta 387 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 5 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 107 | Limitations | tibial artery 584 | `manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv` : 7 | `GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)` |
| 108 | Results 3 | APOB: mipomersen FDA discontinued, EMA refused | `manuscript/display/tableS18_regulatory_status.md` : 15 | `python manuscript_prep\manuscript\revision_sources\fetch_drug_status.py` |
| 109 | Results 3 | PLG: no record in either source (anistreplase) | `manuscript/display/tableS18_regulatory_status.md` : 28 | `python manuscript_prep\manuscript\revision_sources\fetch_drug_status.py` |
| 110 | Results 3 | PLG: no record in either source (streptokinase) | `manuscript/display/tableS18_regulatory_status.md` : 29 | `python manuscript_prep\manuscript\revision_sources\fetch_drug_status.py` |
| 111 | Results 3; Table 2 | IL1B: rilonacept label says IL-1 alpha and IL-1 beta trap | `manuscript/display/table2_named_targets.md` : 15 | `python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py` |
| 112 | Table S22 | UKB-PPP discovery cohort 34,557, European ancestry (no longer the PCSK9 exposure sample) | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 62 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 113 | Results 8; Table S22 | UKB-PPP full cohort 52,363, the PCSK9 exposure sample (body text of the paper; the file is one line) | `manuscript/prior_art_search/fulltext/PMC10567551.txt` : 1 | `python manuscript_prep\manuscript\prior_art_search\papers.py ctx PMC10567551 \"52,363\"` |
| 114 | Results 8; Table S22 | UKB-PPP protein levels inverse-rank normalized | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 59 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 115 | Results 8; Table S22 | Graham 2021: triglycerides log-transformed, residuals inverse-normalized | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 5 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 116 | Table S22 | Graham 2021: UK Biobank participants in the datasets behind the summary statistics (body text; the file is one line) | `manuscript/prior_art_search/fulltext/PMC8730582.txt` : 1 | `python manuscript_prep\manuscript\prior_art_search\papers.py ctx PMC8730582 \"UK Biobank\"` |
| 117 | Results 8; Table S22 | Keaton 2024: effects in mmHg | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 26 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 118 | Table S22 | Keaton 2024 cohorts: UKB, ICBP, MVP, BioVU | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 39 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 119 | Table S22 | Koskeridis 2022: CRP statistics from UK Biobank and CHARGE, MTAG | `manuscript/revision_sources/units_and_overlap_extracts.txt` : 51 | `python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py` |
| 120 | Results 8 | F statistics 2,305 (PCSK9) and 115 (APOC3) | `manuscript/rescue_controls/run_log.txt` : 3 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 121 | Results 8 | F 115 (APOC3) | `manuscript/rescue_controls/run_log.txt` : 8 | `python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py` |
| 122 | Results 3 | PLG: no streptokinase or anistreplase row in four Purple Book releases | `manuscript/revision_sources/purple_book/purple_book_check_output.txt` : 9 | `python manuscript_prep\manuscript\revision_sources\purple_book\purple_book_check.py` |
| 123 | Results 3 | PLG not restored after the Purple Book and Federal Register search | `manuscript/revision_sources/purple_book/purple_book_check_output.txt` : 15 | `python manuscript_prep\manuscript\revision_sources\purple_book\purple_book_check.py` |
| 124 | Results 3; Table 2 | APOB: mipomersen FDA original approval 2013-01-29, NDA 203568; label indication HoFH | `manuscript/APPROVED_35_CHECK.csv` : 102 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 125 | Tables 1, S3, S11 footnote | NR3C1 counted through myocarditis | `manuscript/APPROVED_35_CHECK.csv` : 180 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 126 | Tables 1, S3, S11 footnote | MS4A1 counted through systemic vasculitis | `manuscript/APPROVED_35_CHECK.csv` : 174 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 127 | Results 7; Table S23 | 11 named targets; 9 with a matching GTEx tissue; 1 with a candidate in it (AGTR1) | `manuscript/recovery_out/tableS23_output.txt` : 1 | `python manuscript_prep\manuscript\recovery_scripts\make_table_s23.py` |
| 128 | Results 7; Table S23 | AGTR1: candidate in tibial artery, lead p 2.97e-14 | `manuscript/recovery_out/tableS23_output.txt` : 2 | `python manuscript_prep\manuscript\recovery_scripts\make_table_s23.py` |
| 129 | Methods; Supplement 11 | audit: category identical for all 280 targets | `manuscript/recovery_out/compare_with_audit_output.txt` : 2 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 130 | Methods; Supplement 11 | audit: reporting tier identical for all 280 targets | `manuscript/recovery_out/compare_with_audit_output.txt` : 3 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 131 | Supplement 11 | audit: UKB-PPP lead variant, sheet and -log10 p identical after the fix | `manuscript/recovery_out/compare_with_audit_output.txt` : 10 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 132 | Supplement 11 | no category or tier changed by the fixes | `manuscript/recovery_out/compare_with_audit_output.txt` : 21 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 133 | Supplement 11 | UKB-PPP lead variant changed for 9 targets | `manuscript/recovery_out/compare_with_audit_output.txt` : 22 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 134 | Supplement 11 | deCODE lead variant changed for 2 targets | `manuscript/recovery_out/compare_with_audit_output.txt` : 32 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 135 | Supplement 11 | UKB-PPP sheet changed for 30 targets | `manuscript/recovery_out/compare_with_audit_output.txt` : 36 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 136 | Supplement 11 | Fenland assayed changed for CCN1 and GRK2 | `manuscript/recovery_out/compare_with_audit_output.txt` : 37 | `python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py` |
| 137 | Supplement 11 | deCODE 'not included' reading: 4 targets, 0 tier changes | `manuscript/recovery_map_summary.txt` : 265 | `python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py` |
| 138 | Methods; Supplement 11 | audit: eQTLGen status exact match, 726 of 726 | `../audit/AUDIT.md` : 11 | `second-implementation audit, D:\causal2drug\audit (not part of manuscript_prep)` |
| 139 | Supplement 11 | audit: PCSK9 +0.387, F 2,305; HMGCR match | `../audit/AUDIT.md` : 14 | `second-implementation audit, D:\causal2drug\audit (not part of manuscript_prep)` |
| 140 | Supplement 7; Table S22 | PCSK9 sensitivity, European discovery cohort: MR beta +0.376 (0.367 to 0.385) | `manuscript/rescue_controls/sensitivity_pcsk9_discovery_output.txt` : 3 | `python manuscript_prep\manuscript\rescue_controls\sensitivity_pcsk9_discovery.py` |
| 141 | Supplement 7; Table S22 | PCSK9 sensitivity: F 1,723 | `manuscript/rescue_controls/sensitivity_pcsk9_discovery_output.txt` : 2 | `python manuscript_prep\manuscript\rescue_controls\sensitivity_pcsk9_discovery.py` |
| 142 | Methods | 327 evidence rows, 32 targets reviewed | `drug_links/manual_review_output.txt` : 1 | `python manuscript_prep\drug_links\complete_manual_review.py` |
| 143 | Supplement 4 | 308 rows: drug recorded against a family or group | `drug_links/manual_review_output.txt` : 2 | `python manuscript_prep\drug_links\complete_manual_review.py` |
| 144 | Supplement 4 | 19 rows: direct link, indication outside cardiovascular branches | `drug_links/manual_review_output.txt` : 3 | `python manuscript_prep\drug_links\complete_manual_review.py` |
| 145 | Methods; Supplement 4 | no row accepted; strict rule retained | `drug_links/manual_review_output.txt` : 7 | `python manuscript_prep\drug_links\complete_manual_review.py` |
| 146 | Supplement 4; Results 7 | AVPR2: tolvaptan FDA marketed, label text names heart failure | `manuscript/APPROVED_35_CHECK.csv` : 108 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 147 | Supplement 4 | abciximab: FDA label indication (percutaneous coronary intervention) | `manuscript/APPROVED_35_CHECK.csv` : 156 | `python manuscript_prep\manuscript\revision_sources\check_approved_35.py` |
| 148 | Table 2 | ITGA2B: 2 of 3 linked drugs marketed | `manuscript/display/table2_named_targets.md` : 16 | `python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py` |
| 149 | Abstract; Results 1; Discussion | all 6 lipid targets lack an eQTLGen instrument: 5 with no eQTL data, NPC1L1 tested without one | `manuscript/recovery_out/lead_findings_output.txt` : 7 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 150 | Results 1 | PCSK9 no eQTL data | `manuscript/recovery_out/lead_findings_output.txt` : 1 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 151 | Results 1 | LPA no eQTL data, phase III | `manuscript/recovery_out/lead_findings_output.txt` : 6 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 152 | Results 1 | NPC1L1 tested without a qualifying instrument | `manuscript/recovery_out/lead_findings_output.txt` : 5 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 153 | Results 1 | HMGCR, CETP, LDLR, LPL instrumented | `manuscript/recovery_out/lead_findings_output.txt` : 8 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 154 | Abstract; Results 5; Discussion | heart failure: 21 locked = 4 microRNA + 1 mitochondrial + 16 nuclear protein-coding; 11 of 16 never assayed | `manuscript/recovery_out/lead_findings_output.txt` : 12 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 155 | Results 4 | all areas: 57 locked = 9 microRNA + 2 mitochondrial + 46 protein-coding; 24 never assayed, 22 assayed | `manuscript/recovery_out/lead_findings_output.txt` : 13 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 156 | Results 5; Discussion | DMD locked (assayed, no candidate) | `manuscript/recovery_out/lead_findings_output.txt` : 14 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 157 | Results 5; Discussion | PLN locked (never assayed) | `manuscript/recovery_out/lead_findings_output.txt` : 15 | `python manuscript_prep\manuscript\recovery_scripts\lead_findings.py` |
| 158 | Results 6 | liver: no candidate for PCSK9 | `manuscript/recovery_out/tableS23_output.txt` : 8 | `python manuscript_prep\manuscript\recovery_scripts\make_table_s23.py` |
| 159 | Results 6 | liver: no candidate for ANGPTL3, APOC3, APOB | `manuscript/recovery_out/tableS23_output.txt` : 4 | `python manuscript_prep\manuscript\recovery_scripts\make_table_s23.py` |

## Numbers not carried by a single output line

| where | number | source | how | note |
|---|---|---|---|---|
| Methods | Fenland 10,708 participants | manuscript/recovery_data/pietzner2021_supplementary_tables.xlsx, sheet ST1 | python -I manuscript_prep/manuscript/recovery_scripts/inspect_xlsx.py <workbook> 5 ST1 | row 'Participants [N]' reads 10708 |
| Introduction; Methods | eQTLGen tested 19,250 genes | analysis/classify_output.txt : 185 | python manuscript_prep/analysis/classify.py | the line gives 19,123 autosomal + 127 not in GRCh37.87 = 19,250; the figure is also quoted from the eQTLGen paper in PRIOR_ART.md |
| Introduction; Methods | eQTLGen n = 31,684; UKB-PPP 54,219; deCODE 35,559 | abstracts of the cited papers in PubMed | python manuscript_prep\manuscript\prior_art_search\papers.py lookup \"<title>[ti]\" (abstract text) | verified against the PubMed abstract on 2026-10-11; not stored in a project file |
| Results 5 | 23 of 72 off-target rescues had been measured in a proteomic study | manuscript/recovery_map.csv | count rows with reporting_tier = 'candidate: other tissue only' and non-empty pqtl_studies_assayed | computed on 2026-10-11 from the CSV; the count is not printed in recovery_map_summary.txt |
| Results 5 | 9 microRNA genes among the never-assayed locked targets | manuscript/recovery_map_summary.txt [R30] | count of names beginning MIR in the [R30] never-assayed list | counted from the printed list |
| Results 7 | 5 protein-only rescues among the 13 are secreted proteins | analysis/out/per_target_primary.csv, column `secreted` | rows ANGPTL3, APOC3, PCSK9, PLG, REN | all five read YES |
| Discussion | 16 of 21 locked heart-failure targets never assayed | manuscript/recovery_map_summary.txt [R25] | Heart failure row | same line as the 70-target row above |
| Results 8; Table S22 | overlap upper limits 4.0%, 0.8%, 5.1% of the outcome sample | sample sizes above | 52,363 / 1,320,016 = 3.97% (PCSK9, ANGPTL3); 10,708 / 1,320,016 = 0.81% (APOC3); 52,363 / 1,028,980 = 5.09% (REN) | upper limits, not measured overlap |
| Table S22 | Fenland effects per 1 SD of normalized plasma abundance | PMC9904207 main text (prior_art_search/fulltext/PMC9904207.txt) | phrase '1 s.d. higher normalized plasma abundance' | the scale of the supplementary table itself is not stated in the workbook |
| Results 3 | 11 named targets; PLG and IL1B shown in Table 2 only | manuscript/display/table2_named_targets.md | 13 rows minus the two rows marked with a dagger or section sign |  |
| Results 7 | 11 named launched targets: 7 / 3 / 1 by tier | manuscript/display/table2_named_targets.md | the 13 rows of Table 2 without IL1B (other-tissue tier) and PLG (CV-relevant or plasma-protein tier) | 13 rows give 8 / 4 / 1 ([R26]); removing the two gives 7 / 3 / 1 |
| Results 3 | a currently marketed drug for 10 of the 11 named targets | manuscript/display/table2_named_targets.md, column 'Regulatory status (FDA / EMA)' | count named rows that begin 'Marketed' | AGTR1, ANGPTL3, APOC3 (EMA only), AVPR2, F2, PCSK9, REN, SLC5A1, ITGA2B, NPC1L1; not APOB |
| Results 8 | REN lead pQTL frequency 0.8% | manuscript/rescue_controls/exposure_rows.json | REN_SBP, af_effect | 0.00762883 |
| Results 8; Discussion | APOC3 variant 51 kb from the gene | manuscript/recovery_map.csv, row APOC3, fenland_bp_from_gene_body |  | 51505 |
| Discussion | heart tissue samples in GTEx about 80 times smaller than eQTLGen | gtex_v8_tissue_sample_sizes.csv | 31,684 / 386 = 82; 31,684 / 372 = 85 | rounded |
| Results 5 | kidney cortex is the smallest tissue sample | gtex_v8_tissue_sample_sizes.csv | minimum eqtl_samples among tissues with eGenes | 73 |
| Results 7 | sites of drug action for SLC5A1, ITGA2B, NPC1L1 | manuscript/revision_sources/drug_status.csv, column label_mechanism_text; drug_status/extra_openfda_brand_name_Zetia.json | FDA label mechanism-of-action text | quoted in Table 2 |

Entries traced: 159; UNVERIFIED: 0.
