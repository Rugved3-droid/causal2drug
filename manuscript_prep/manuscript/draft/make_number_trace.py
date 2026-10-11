#!/usr/bin/env python
"""Write NUMBER_TRACE.md: every number in the manuscript, the output file and line that carries
it, and the command that regenerates that file. A number whose pattern is not found in the
stated file is written as UNVERIFIED, so the trace cannot silently drift from the outputs.

    python manuscript_prep/manuscript/draft/make_number_trace.py
"""
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
F = {
    "K": ("analysis/classify_output.txt", r"python manuscript_prep\analysis\classify.py"),
    "R": ("manuscript/recovery_map_summary.txt",
          r"python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py"),
    "V": ("verify/verify_output.txt", r"python manuscript_prep\verify\verify_numbers.py"),
    "T1": ("tables/table1_flow_and_cohort.md", r"python manuscript_prep\tables\make_tables.py"),
    "S3": ("tables/tableS3_sensitivity_cv_set.md", r"python manuscript_prep\tables\make_tables.py"),
    "S4": ("tables/tableS4_sensitivity_cis_window.md",
           r"python manuscript_prep\tables\make_tables.py"),
    "S9": ("tables/tableS9_enrichment.md", r"python manuscript_prep\tables\make_tables.py"),
    "S10": ("tables/tableS10_no_eqtl_data_by_chromosome.md",
            r"python manuscript_prep\tables\make_tables.py"),
    "S11": ("tables/table2_status_by_phase.md", r"python manuscript_prep\tables\make_tables.py"),
    "T2": ("manuscript/display/table2_named_targets.md",
           r"python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py"),
    "T3": ("manuscript/display/table3_recovery_tiers.md",
           r"python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py"),
    "S16": ("manuscript/display/tableS16_recovery_by_resource.md",
            r"python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py"),
    "MAP": ("manuscript/recovery_map.csv",
            r"python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py"),
    "PA": ("manuscript/PRIOR_ART.md", "quotation from the cited paper's full text; see "
           r"manuscript\prior_art_search\fulltext"),
    "PCSK9": ("verify/pcsk9_scan.txt", r"see header of verify\pcsk9_scan.txt"),
    "KN": ("KEY_NUMBERS.md", r"python manuscript_prep\analysis\classify.py (summarised in "
           "KEY_NUMBERS.md)"),
    "LIM": ("LIMITATIONS.md", r"python manuscript_prep\analysis\classify.py (summarised in "
            "LIMITATIONS.md)"),
    "FEN": ("manuscript/DATA_ACCESS.md", "Fenland workbook, sheet ST1, row 'Participants [N]'"),
    "PCI": ("analysis/phase_ci_output.txt", r"python manuscript_prep\analysis\phase_ci.py"),
    "RC": ("manuscript/rescue_controls/run_log.txt",
           r"python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py"),
    "IL6": ("manuscript/il6r_audit/il6r_audit_output.txt",
            r"python manuscript_prep\manuscript\il6r_audit\il6r_audit.py"),
    "GN": ("manuscript/revision_sources/gtex_v8_tissue_sample_sizes.csv",
           "GTEx Portal API, tissueSiteDetail, datasetId=gtex_v8 (cached JSON beside the CSV)"),
    "S18": ("manuscript/display/tableS18_regulatory_status.md",
            r"python manuscript_prep\manuscript\revision_sources\fetch_drug_status.py"),
    "UO": ("manuscript/revision_sources/units_and_overlap_extracts.txt",
           r"python manuscript_prep\manuscript\revision_sources\units_and_overlap_lookup.py"),
    "FT_SUN": ("manuscript/prior_art_search/fulltext/PMC10567551.txt",
               r"python manuscript_prep\manuscript\prior_art_search\papers.py ctx PMC10567551 \"52,363\""),
    "FT_GRAHAM": ("manuscript/prior_art_search/fulltext/PMC8730582.txt",
                  r"python manuscript_prep\manuscript\prior_art_search\papers.py ctx PMC8730582 \"UK Biobank\""),
    "AP": ("manuscript/revision_sources/approved_35_output.txt",
           r"python manuscript_prep\manuscript\revision_sources\check_approved_35.py"),
    "PB": ("manuscript/revision_sources/purple_book/purple_book_check_output.txt",
           r"python manuscript_prep\manuscript\revision_sources\purple_book\purple_book_check.py"),
    "A35": ("manuscript/APPROVED_35_CHECK.csv",
            r"python manuscript_prep\manuscript\revision_sources\check_approved_35.py"),
    "S23": ("manuscript/recovery_out/tableS23_output.txt",
            r"python manuscript_prep\manuscript\recovery_scripts\make_table_s23.py"),
    "AUD": ("manuscript/recovery_out/compare_with_audit_output.txt",
            r"python manuscript_prep\manuscript\recovery_scripts\compare_with_audit.py"),
    "AUDMD": ("../audit/AUDIT.md", "second-implementation audit, D:\\causal2drug\\audit (not part of manuscript_prep)"),
    "SENS": ("manuscript/rescue_controls/sensitivity_pcsk9_discovery_output.txt",
             r"python manuscript_prep\manuscript\rescue_controls\sensitivity_pcsk9_discovery.py"),
    "DEV": ("DEVIATION_CLASS_C_FIGURE.md", r"python manuscript_prep\verify\verify_numbers.py"),
}
# (section, statement, key, regex that must match one line of the file)
N = [
 ("Abstract; Results 1", "726 targets; 443 instrumented (61.0%); 217 no eQTL data (29.9%); 63 no qualifying instrument (8.7%); upper bound 280 (38.6%); 3 unresolved; autosomal-only 188 of 697 (27.0%)", "K",
  r"^PRIMARY\s+N= 726 \| instrumented\s+443 \( 61\.0%\) \| no eQTL data\s+217 \( 29\.9%\).*63 \( 8\.7%\).*280 \( 38\.6%\).*unresolved 3.*188 / 697 \(27\.0%\)"),
 ("Methods; Results 4", "27 X-linked and 2 mitochondrial targets", "K", r"^PRIMARY\s+N= 726.*X-linked 27 \+ mitochondrial 2"),
 ("Methods", "29,476 rows; 29,460 pairs; 29,456 pairs; 2,517 targets; 1,014 indications", "T1", r"Distinct targets / indications in those pairs \| 2517 / 1014"),
 ("Methods", "29,456 pairs after removing 4 programmes", "T1", r"removing 4 programmes.*29456"),
 ("Methods", "93 indications, 2,128 pairs, 726 targets; sensitivity sets 630, 663, 594", "K", r"^PRIMARY\s+terms\s+93\s+pairs\s+2128\s+targets\s+726"),
 ("Methods", "sensitivity sets 630, 663 and 594 targets", "T1", r"Cardiovascular targets \| 726 \| 630 \| 663 \| 594"),
 ("Methods; Results 8", "127,341,798 rows in the full eQTLGen file", "PCSK9", r"127,?341,?798"),
 ("Results 1", "lower bound 27.6% to 29.9%, upper bound 36.9% to 38.6% across set definitions", "T1", r"No eQTL data \(lower bound\) \| 217 \(29\.9%\) \| 180 \(28\.6%\) \| 191 \(28\.8%\) \| 164 \(27\.6%\)"),
 ("Results 1", "upper bounds 38.6%, 37.6%, 37.9%, 36.9%", "T1", r"upper bound \| 280 \(38\.6%\) \| 237 \(37\.6%\) \| 251 \(37\.9%\) \| 219 \(36\.9%\)"),
 ("Results 1", "all 2,517 targets: 26.3% and 35.7%", "KN", r"All 2,517 targets \| 2,517 .*663 \(26\.3 %\).*899 \(35\.7 %\)"),
 ("Results 1", "±1 Mb window moves 1 target, IL1B", "K", r"gain an instrument in a wider window: \[\{'gene': 'IL1B'"),
 ("Results 2", "no eQTL data by phase: 30.0%, 28.6%, 33.3%, 18.7%, 34.8%", "S11", r"^\| Preclinical \| 323 \| 196 \(60\.7%\) \| 97 \(30\.0%\)"),
 ("Results 2", "phase I 28.6%", "S11", r"^\| Phase I \| 91 .* 26 \(28\.6%\)"),
 ("Results 2", "phase II 33.3%", "S11", r"^\| Phase II \| 171 .* 57 \(33\.3%\)"),
 ("Results 2", "phase III 18.7%", "S11", r"^\| Phase III \| 75 .* 14 \(18\.7%\)"),
 ("Results 2", "launched 34.8%", "S11", r"^\| Launched \| 66 .* 23 \(34\.8%\)"),
 ("Abstract; Results 2", "chi-square 6.33, df 4, P=0.18", "KN", r"χ² = 6\.33, df = 4, p = 0\.18"),
 ("Results 3", "35 launched targets pass the ChEMBL strict rule", "S11", r"strict\s*\| 35 \| 22 \(62\.9%\) \| 10 \(28\.6%\).*3 \(8\.6%\) \| 13 \(37\.1%\)"),
 ("Results 3", "33 launched targets with a drug ever approved by FDA or EMA; 9 no eQTL data (27.3%); 3 no qualifying instrument; upper bound 12 (36.4%)", "S11", r"ever approved by FDA or EMA for a CV indication\s*\| 33 \| 21 \(63\.6%\) \| 9 \(27\.3%\).*3 \(9\.1%\) \| 12 \(36\.4%\)"),
 ("Results 3", "33 of 35 have an FDA or EMA approval record", "AP", r"^\[Q02\] targets with an FDA or EMA approval record for >=1 strict-rule drug: 33 of 35"),
 ("Results 3", "MTOR and PLG have no approval record", "AP", r"^\[Q03\] targets that DROP OUT .*MTOR \(INSTRUMENTED; drugs: Perhexiline\); PLG \(NO_EQTL_DATA"),
 ("Results 3", "same count from the classifier", "K", r"verified, APPROVED\s+N=\s+33 \| instrumented\s+21 \( 63\.6%\) \| no eQTL data\s+9 \( 27\.3%\).*upper bound\s+12 \( 36\.4%\)"),
 ("Results 3", "permissive rule: 20 and 4 of 58", "S11", r"permissive\s*\| 58 \| 34 \(58\.6%\) \| 20 \(34\.5%\).*4 \(6\.9%\)"),
 ("Results 3", "names of the 13 targets", "KN", r"\*\*No eQTL data\*\* \| AGTR1, ANGPTL3, APOB, APOC3, AVPR2 \(X-linked\), F2, PCSK9, PLG, REN, SLC5A1"),
 ("Results 3", "IL1B, ITGA2B, NPC1L1 tested without a qualifying instrument", "KN", r"\*\*No qualifying instrument\*\* \| IL1B, ITGA2B, NPC1L1"),
 ("Abstract; Results 4", "188 autosomal no-eQTL-data targets; 182 with a GTEx value; 180 below 1 TPM", "S10", r"^\| Primary set \| 726 \| 217 \| 188 \| 27 \| 2 \| 0 \| 29 \| 12 \| 2 \| 180 / 182"),
 ("Results 4", "12 of 27 X-linked targets expressed in blood at ≥1 TPM", "S10", r"^\| Primary set \| 726 \| 217 \| 188 \| 27 \| 2 \| 0 \| 29 \| 12 \|"),
 ("Results 4", "tissue restriction OR 7.47 (4.12–13.55)", "K", r"^PRIMARY\s+NO_EQTL_DATA\s+vs INSTRUMENTED x tissue_restricted OR\s+7\.47 \[ 4\.12,\s+13\.55\]"),
 ("Results 4", "autosomal only OR 12.18 (5.58–26.62)", "KN", r"autosomal only, OR 12\.18 \(5\.58–26\.62\)"),
 ("Results 4", "no qualifying instrument OR 1.39 (0.77–2.51)", "K", r"^PRIMARY\s+NO_QUALIFYING_INSTRUMENT\s+vs INSTRUMENTED x tissue_restricted OR\s+1\.39 \[ 0\.77,\s+2\.51\]"),
 ("Methods", "GTEx v8: 49 tissues", "R", r"^\[R02\] GTEx v8: tissues with target rows 49"),
 ("Methods", "UKB-PPP threshold P<1.7×10⁻¹¹ (workbook sheet title)", "FEN", r"1\.7 × 10⁻¹¹"),
 ("Methods", "3 targets matched after an identifier change; 1 changed category", "R", r"^\[R31\] amendment 8, before -> after: ADRA2B: locked -> candidate-tissue"),
 ("Abstract; Results 5", "280 in scope; 223 fixable (79.6%); 57 locked", "R", r"all in-scope targets\s+280\s+128\s+24\s+71\s+57\s+223 \(79\.6%\)"),
 ("Abstract; Results 5", "tiers: 151 (53.9%), 72 (25.7%), 35 never assayed, 22 assayed; locked 57 (20.4%)", "R", r"all in-scope targets\s+280\s+151 \(53\.9%\)\s+72 \(25\.7%\)\s+35 \(12\.5%\)\s+22 \(7\.9%\).*57 \(20\.4%\)"),
 ("Results 5", "GTEx non-blood instrument for 199 targets", "R", r"^\[R16\].*non-blood tissue: 199"),
 ("Results 5", "plasma cis-pQTL instrument for 95 targets", "R", r"^\[R20\].*instrument in >=1 study 95"),
 ("Results 5", "candidate: CV-relevant tissue or plasma protein: 56 tissue only, 60 protein only, 35 both", "R", r"^\[R27\] candidate: CV-relevant tissue or plasma protein \(151\): through a CV-relevant tissue only 56; through a plasma cis-pQTL only 60; through both 35"),
 ("Results 5", "tibial artery 44", "R", r"^\s+Artery_Tibial\s+tested\s+\d+\s+instrument\s+44"),
 ("Results 5", "aorta 34", "R", r"^\s+Artery_Aorta\s+tested\s+\d+\s+instrument\s+34"),
 ("Results 5", "atrial appendage 29", "R", r"^\s+Heart_Atrial_Appendage\s+tested\s+\d+\s+instrument\s+29"),
 ("Results 5", "kidney cortex 0", "R", r"^\s+Kidney_Cortex\s+tested\s+\d+\s+instrument\s+0"),
 ("Results 5", "off-target tissues: testis 21, tibial nerve 17, fibroblasts 16, subcutaneous adipose 14", "R", r"^\[R29\].*Testis 21; Nerve_Tibial 17; Cells_Cultured_fibroblasts 16; Adipose_Subcutaneous 14"),
 ("Results 5", "27 of 72 on a single tissue; 7 on testis alone", "R", r"^\[R29\].*single tissue: 27; on testis only: 7"),
 ("Results 5", "never assayed 35 (list); 9 microRNA, 2 mitochondrial", "R", r"^\[R30\] locked: never assayed \(35\)"),
 ("Results 5", "assayed, no instrument 22 (list)", "R", r"^\[R30\] locked: assayed, no instrument \(22\)"),
 ("Results 5", "25 of 35 never-assayed targets tested in ≥1 GTEx tissue (35 minus the 10 in no egenes file)", "R", r"^\[R21\].*not in any GTEx egenes file: 10"),
 ("Results 5", "no qualifying instrument: locked 14 of 63 (22.2%)", "R", r"^\s+no qualifying instrument\s+63\s+31.*total 14 \(22\.2%\)"),
 ("Results 5", "autosomal no data: locked 37 of 188 (19.7%)", "R", r"autosomal, not tested\s+188\s+105.*total 37 \(19\.7%\)"),
 ("Results 5", "X-linked: 23 of 27 with an instrument, 4 locked", "R", r"X-linked, not analysed\s+27\s+15 \(55\.6%\)\s+8 \(29\.6%\).*total 4 \(14\.8%\)"),
 ("Results 6", "locked by phase: 27.0%, 16.7%, 15.2%, 24.0%, 3.7%", "R", r"^\s+Preclinical\s+126\s+65.*total 34 \(27\.0%\)"),
 ("Results 6", "phase I 16.7%", "R", r"^\s+Phase I\s+36\s+18.*total 6 \(16\.7%\)"),
 ("Results 6", "phase II 15.2%", "R", r"^\s+Phase II\s+66\s+35.*total 10 \(15\.2%\)"),
 ("Results 6", "phase III 24.0%", "R", r"^\s+Phase III\s+25\s+16.*total 6 \(24\.0%\)"),
 ("Results 6", "launched 3.7%", "R", r"^\s+Launched\s+27\s+17.*total 1 \(3\.7%\)"),
 ("Abstract; Results 6; Discussion", "heart failure / cardiomyopathy: 70 targets; 25 (35.7%), 24 (34.3%), locked 21 (30.0%), 16 never assayed", "R", r"Heart failure / cardiomyopathy\s+70\s+25 \(35\.7%\)\s+24 \(34\.3%\)\s+16 \(22\.9%\).*total 21 \(30\.0%\)"),
 ("Results 6", "thrombosis: 0 of 22 locked; 17 candidate: CV-relevant tissue or plasma protein", "R", r"Thrombosis\s+22\s+17 \(77\.3%\).*total 0 \(0\.0%\)"),
 ("Results 7", "13 targets of Table 2: 8, 4 and 1 by tier (the text reports the 12 named, without IL1B: 8, 3, 1)", "R", r"the 13 named targets\s+13\s+8 \(61\.5%\)\s+4 \(30\.8%\)\s+1 \(7\.7%\)"),
 ("Results 7", "tissues and pQTL status of each named target", "T2", r"^\| AGTR1 \|"),
 ("Results 7", "AVPR2: locked, never assayed; tested in 45 GTEx tissues", "MAP", r"^AVPR2,ENSG00000126895,.*locked: never assayed.*,45,0,"),
 ("Results 7", "HCN4: locked, never assayed; tested in 25 GTEx tissues", "MAP", r"^HCN4,.*locked: never assayed.*,25,0,"),
 ("Results 7", "AVPR1A, SLC6A2, PDE3A: candidate in a CV-relevant tissue or plasma protein (key 'candidate: CV-relevant tissue or plasma protein' in the CSV)", "MAP", r"^PDE3A,.*candidate: CV-relevant tissue or plasma protein"),
 ("Results 8", "HMGCR β 0.389, P=4.2×10⁻²²⁹", "KN", r"HMGCR → LDL-C \| MR β \+0\.389 \(SE 0\.012\), p = 4\.2 × 10⁻²²⁹"),
 ("Results 8", "IL6R β 0.376, P=1.2×10⁻⁵⁹, fails on direction", "KN", r"IL6R → CRP \| MR β \+0\.376 \(SE 0\.023\), p = 1\.2 × 10⁻⁵⁹"),
 ("Results 8", "PCSK9: 0 of 127,341,798 rows", "KN", r"PCSK9 \| 0 of 127,341,798 rows"),
 ("Limitations", "26 of 63 with a cis-eQTL at FDR<0.05", "K", r"T2 cis-eQTL at FDR<0\.05, none at p<5e-8': 26"),
 ("Limitations", "46 of 93 indications in the 'other' area", "LIM", r"\"Other\" holds 46 of 93 terms"),
 ("Introduction", "Storm 2021: 2,786 of 4,863", "PA", r"Storm 2021.*2786.*4,863"),
 ("Introduction", "Liang 2025: 2,888 of 4,479", "PA", r"Liang 2025.*4,479.*2,888"),
 ("Introduction", "Liu 2024: 2,525 of 4,479", "PA", r"Liu 2024.*2,525.*4,479"),
 ("Introduction", "Wang 2024: 2,531 of 4,317", "PA", r"Wang 2024.*2531 of the 4317"),
 ("Introduction", "Sadler 2023: ~12,550 of ~19,150; HCN4", "PA", r"Sadler 2023.*12,550 of ~19,150"),
 ("Introduction", "Gaziano 2021: 1,016 of 1,263; 49 tissues", "PA", r"Gaziano 2021.*49 tissues.*1,016 of the 1,263"),
 ("Introduction", "eQTLGen: 19,250 autosomal genes expressed in blood and tested", "PA", r"19,250 autosomal genes expressed in blood and tested"),
 ("Supplement", "Class C: 50 paired, 32 higher, 18 lower, medians 0.00821 and 0.01390, Wilcoxon P=0.0262, binomial P=0.0649", "DEV", r"50 paired; PP4 higher in eQTLGen at 32, lower at 18"),
 ("Supplement", "Fenland build check: 3,917 of 3,933 identical to GRCh37.87; 11 to GENCODE v26", "R", r"^\[R08\].*3933 compared.*for 3917.*for 11"),
 ("Supplement", "X coverage: GTEx 1,293 X-linked genes", "R", r"^\[R02\].*chromosome X 1293"),
 ("Supplement", "amendment 8 before: 222 fixable, 58 locked", "R0", r"all in-scope targets\s+280\s+127\s+24\s+71\s+58\s+222 \(79\.3%\)"),
]
N += [
 ("Abstract; Results 2", "no eQTL data, preclinical 30.0% (95% CI 25.3-35.2)", "PCI", r"^\s+Preclinical\s+97 / 323 =\s+30\.0% \(95% CI 25\.3% to 35\.2%\)"),
 ("Results 2", "phase I 28.6% (20.3-38.6)", "PCI", r"^\s+Phase I\s+26 /\s+91 =\s+28\.6% \(95% CI 20\.3% to 38\.6%\)"),
 ("Results 2", "phase II 33.3% (26.7-40.7)", "PCI", r"^\s+Phase II\s+57 / 171 =\s+33\.3% \(95% CI 26\.7% to 40\.7%\)"),
 ("Results 2", "phase III 18.7% (11.5-28.9)", "PCI", r"^\s+Phase III\s+14 /\s+75 =\s+18\.7% \(95% CI 11\.5% to 28\.9%\)"),
 ("Results 2", "launched 34.8% (24.5-46.9)", "PCI", r"^\s+Launched\s+23 /\s+66 =\s+34\.8% \(95% CI 24\.5% to 46\.9%\)"),
 ("Abstract; Results 2", "launched minus preclinical +4.8 points (95% CI -6.8 to +17.8)", "PCI", r"Launched minus Preclinical: \+4\.8 percentage points \(95% CI -6\.8 to \+17\.8\)"),
 ("Results 2", "upper bound difference +1.9 points (-10.4 to +15.0)", "PCI", r"Launched minus Preclinical: \+1\.9 percentage points \(95% CI -10\.4 to \+15\.0\)"),
 ("Abstract; Results 8", "PCSK9 candidate control: rs11591147, MR beta +0.387 (0.377 to 0.396); PASS", "RC", r"MR beta \+0\.3868 \(95% CI \+0\.3774 to \+0\.3963\)"),
 ("Results 8", "PCSK9 control verdict", "RC", r"^PCSK9_LDL .*rs11591147\s+-> PASS"),
 ("Abstract; Results 8", "APOC3 candidate control: rs964184, MR beta +1.171 (1.153 to 1.190); PASS", "RC", r"MR beta \+1\.1711 \(95% CI \+1\.1526 to \+1\.1897\)"),
 ("Results 8", "APOC3 control verdict", "RC", r"^APOC3_TG .*rs964184\s+-> PASS"),
 ("Abstract; Results 8", "AGTR1 candidate control: MR beta -0.037 (-0.268 to 0.194), P=0.75; FAIL", "RC", r"MR beta -0\.0372 \(95% CI -0\.2682 to \+0\.1937\), p 0\.752"),
 ("Results 8", "AGTR1 control verdict", "RC", r"^AGTR1_SBP .*-> FAIL"),
 ("Abstract; Results 8", "ANGPTL3 control not harmonizable (indel)", "RC", r"^ANGPTL3_TG .*-> NOT HARMONIZABLE"),
 ("Results 8", "ANGPTL3 reason: not_snp", "RC", r"harmonize_one status: not_snp"),
 ("Abstract; Results 8", "REN control not harmonizable (absent from outcome file)", "RC", r"^REN_SBP .*-> NOT HARMONIZABLE"),
 ("Results 8", "REN reason", "RC", r"variant position absent from the outcome file"),
 ("Results 8", "IL6R audit: MR variant rs2229238; eQTLGen lead rs10908839", "IL6", r"^\[A01\] control as run: eQTLGen lead rs10908839 .*MR variant rs2229238"),
 ("Results 8", "IL6R audit: r2 = 0.16 between rs2229238 and rs2228145", "IL6", r"rs2229238\s+r = \+0\.401, r2 = 0\.161"),
 ("Results 8", "IL6R audit: 358Ala allele, expression z -3.17", "IL6", r"rs2228145\s+chr1:154426970 assessed C other A z\s+-3\.17"),
 ("Results 8", "IL6R audit: harmonization verified by hand", "IL6", r"^\[A08\] hand check of the MR variant rs2229238"),
 ("Limitations", "GTEx v8 eQTL samples: kidney cortex 73", "GN", r"^Kidney_Cortex,.*,73,"),
 ("Limitations", "liver 208", "GN", r"^Liver,Liver,208,"),
 ("Limitations", "coronary artery 213", "GN", r"^Artery_Coronary,.*,213,"),
 ("Limitations", "adrenal gland 233", "GN", r"^Adrenal_Gland,.*,233,"),
 ("Limitations", "heart atrial appendage 372", "GN", r"^Heart_Atrial_Appendage,.*,372,"),
 ("Limitations", "heart left ventricle 386", "GN", r"^Heart_Left_Ventricle,.*,386,"),
 ("Limitations", "aorta 387", "GN", r"^Artery_Aorta,.*,387,"),
 ("Limitations", "tibial artery 584", "GN", r"^Artery_Tibial,.*,584,"),
 ("Results 3", "APOB: mipomersen FDA discontinued, EMA refused", "S18", r"^\| APOB \| mipomersen \| Discontinued: 1 product\(s\) \| discontinued \| Kynamro: Refused"),
 ("Results 3", "PLG: no record in either source (anistreplase)", "S18", r"^\| PLG \| anistreplase \| no record in Drugs@FDA \(openFDA\) \| no FDA record \| no record in the EMA medicines table"),
 ("Results 3", "PLG: no record in either source (streptokinase)", "S18", r"^\| PLG \| streptokinase \| no record in Drugs@FDA \(openFDA\) \| no FDA record \| no record in the EMA medicines table"),
 ("Results 3; Table 2", "IL1B: rilonacept label says IL-1 alpha and IL-1 beta trap", "T2", r"IL1B.*Soluble IL-1.* and IL-1.*decoy receptor"),
]
N += [
 ("Table S22", "UKB-PPP discovery cohort 34,557, European ancestry (no longer the PCSK9 exposure sample)", "UO", r"participants of European ancestry from the randomly selected baseline cohort \( n = 34,557\)"),
 ("Results 8; Table S22", "UKB-PPP full cohort 52,363, the PCSK9 exposure sample (body text of the paper; the file is one line)", "FT_SUN", r"full UKB-PPP cohort \( n = 52,363\)"),
 ("Results 8; Table S22", "UKB-PPP protein levels inverse-rank normalized", "UO", r"Each protein level was inverse-rank normalized"),
 ("Results 8; Table S22", "Graham 2021: triglycerides log-transformed, residuals inverse-normalized", "UO", r"Triglyceride levels were natural-log transformed before generating residuals\. Inverse normalization"),
 ("Table S22", "Graham 2021: UK Biobank participants in the datasets behind the summary statistics (body text; the file is one line)", "FT_GRAHAM", r"UK Biobank individuals not present in datasets used to generate the summary statistics"),
 ("Results 8; Table S22", "Keaton 2024: effects in mmHg", "UO", r"measured effect in the meta-analysis \(mmHg\)"),
 ("Table S22", "Keaton 2024 cohorts: UKB, ICBP, MVP, BioVU", "UO", r"combining all available genetic data from the UKB, ICBP and MVP"),
 ("Table S22", "Koskeridis 2022: CRP statistics from UK Biobank and CHARGE, MTAG", "UO", r"CRP summary statistics were obtained from the GWAS meta-analysis of UKB and CHARGE"),
 ("Results 8", "F statistics 2,305 (PCSK9) and 115 (APOC3)", "RC", r"exposure beta -1\.1211 \(SE 0\.0234, F 2305\)"),
 ("Results 8", "F 115 (APOC3)", "RC", r"exposure beta -0\.2103 \(SE 0\.0196, F 115\)"),
]
N += [
 ("Results 3", "PLG: no streptokinase or anistreplase row in four Purple Book releases", "PB", r"^\[B04\] .*four Purple Book releases: 0"),
 ("Results 3", "PLG not restored after the Purple Book and Federal Register search", "PB", r"^\[B06\] conclusion: no FDA approval record for streptokinase or anistreplase"),
 ("Results 3; Table 2", "APOB: mipomersen FDA original approval 2013-01-29, NDA 203568; label indication HoFH", "A35", r"^APOB,Mipomersen,FDA,.*original approval 2013-01-29, NDA 203568"),
 ("Tables 1, S3, S11 footnote", "NR3C1 counted through myocarditis", "A35", r"^NR3C1,Dexamethasone,FDA,.*Myocarditis \[the cardiovascular indication relied on\]"),
 ("Tables 1, S3, S11 footnote", "MS4A1 counted through systemic vasculitis", "A35", r"^MS4A1,Rituximab,FDA,.*Systemic Vasculitis \[the cardiovascular indication relied on\]"),
]
N += [
 ("Results 7; Table S23", "11 named targets; 9 with a matching GTEx tissue; 1 with a candidate in it (AGTR1)", "S23", r"^\[S23\] named targets 11; with a matching GTEx tissue 9; with a candidate instrument in a matching tissue 1 \(AGTR1\)"),
 ("Results 7; Table S23", "AGTR1: candidate in tibial artery, lead p 2.97e-14", "S23", r"^\[S23\] AGTR1 .*Artery Tibial \(n=584\): YES, lead P 2\.97"),
]
N += [
 ("Methods; Supplement 11", "audit: category identical for all 280 targets", "AUD", r"^\[X02\] category\s+differences: 0"),
 ("Methods; Supplement 11", "audit: reporting tier identical for all 280 targets", "AUD", r"^\[X02\] reporting tier\s+differences: 0"),
 ("Supplement 11", "audit: UKB-PPP lead variant, sheet and -log10 p identical after the fix", "AUD", r"^\[X02\] UKB-PPP lead variant \(GRCh38 position\)\s+differences: 0"),
 ("Supplement 11", "no category or tier changed by the fixes", "AUD", r"^\[X03\] tier changed by the fixes: 0"),
 ("Supplement 11", "UKB-PPP lead variant changed for 9 targets", "AUD", r"^\[X04\] ukbppp lead variant changed by the fixes: 9"),
 ("Supplement 11", "deCODE lead variant changed for 2 targets", "AUD", r"^\[X04\] decode lead variant changed by the fixes: 2"),
 ("Supplement 11", "UKB-PPP sheet changed for 30 targets", "AUD", r"^\[X05\] UKB-PPP sheet of the lead changed by the fixes: 30"),
 ("Supplement 11", "Fenland assayed changed for CCN1 and GRK2", "AUD", r"^\[X06\] Fenland assayed changed by the fixes: 2 CCN1, GRK2"),
 ("Supplement 11", "deCODE 'not included' reading: 4 targets, 0 tier changes", "R", r"^\[R32\] targets whose deCODE aptamers are all marked not included: 4 .*would change if they were counted as assayed: 0"),
 ("Methods; Supplement 11", "audit: eQTLGen status exact match, 726 of 726", "AUDMD", r"CV set, eQTLGen status: exact match, 726 of 726 targets"),
 ("Supplement 11", "audit: PCSK9 +0.387, F 2,305; HMGCR match", "AUDMD", r"PCSK9 → LDL-C: MISMATCH.*mine \+0\.387 \(SE 0\.0048, F 2,305\)"),
]
N += [
 ("Supplement 7; Table S22", "PCSK9 sensitivity, European discovery cohort: MR beta +0.376 (0.367 to 0.385)", "SENS", r"^\[D03\] sensitivity MR beta \+0\.3761 \(95% CI \+0\.3668 to \+0\.3853\)"),
 ("Supplement 7; Table S22", "PCSK9 sensitivity: F 1,723", "SENS", r"^\[D02\] .*F 1723"),
]
F["MR"] = ("drug_links/manual_review_output.txt",
           r"python manuscript_prep\drug_links\complete_manual_review.py")
N += [
 ("Methods", "327 evidence rows, 32 targets reviewed", "MR", r"^\[MR01\] evidence rows reviewed: 327; targets: 32"),
 ("Supplement 4", "308 rows: drug recorded against a family or group", "MR", r"^\[MR02\] .*not the gene: 308"),
 ("Supplement 4", "19 rows: direct link, indication outside cardiovascular branches", "MR", r"^\[MR03\] .*branches: 19"),
 ("Methods; Supplement 4", "no row accepted; strict rule retained", "MR", r"^\[MR05\] .*primary definition: 0"),
 ("Supplement 4; Results 7", "AVPR2: tolvaptan FDA marketed, label text names heart failure", "A35", r"^AVPR2,Tolvaptan,FDA,.*yes: heart failure"),
 ("Supplement 4", "abciximab: FDA label indication (percutaneous coronary intervention)", "A35", r"^ITGAV,Abciximab,FDA,.*yes: coronary, cardiac, angina"),
 ("Table 2", "ITGA2B: 2 of 3 linked drugs marketed", "T2", r"^\| ITGA2B \|.*Marketed: 2 of 3 drugs"),
]
F["LEAD"] = ("manuscript/recovery_out/lead_findings_output.txt",
             r"python manuscript_prep\manuscript\recovery_scripts\lead_findings.py")
N += [
 ("Abstract; Results 1; Discussion", "all 6 lipid targets lack an eQTLGen instrument: 5 with no eQTL data, NPC1L1 tested without one", "LEAD", r"^\[L02\] of these 6 lipid targets, without an eQTLGen whole-blood instrument: 6; with no eQTL data: 5; tested without a qualifying instrument: 1"),
 ("Results 1", "PCSK9 no eQTL data", "LEAD", r"^\[L01\] PCSK9\s+eQTLGen: NO_EQTL_DATA"),
 ("Results 1", "LPA no eQTL data, phase III", "LEAD", r"^\[L01\] LPA\s+eQTLGen: NO_EQTL_DATA\s+phase: Phase III"),
 ("Results 1", "NPC1L1 tested without a qualifying instrument", "LEAD", r"^\[L01\] NPC1L1\s+eQTLGen: NO_QUALIFYING_INSTRUMENT"),
 ("Results 1", "HMGCR, CETP, LDLR, LPL instrumented", "LEAD", r"^\[L03\] for contrast, HMGCR\s+eQTLGen: INSTRUMENTED"),
 ("Abstract; Results 5; Discussion", "heart failure: 21 locked = 4 microRNA + 1 mitochondrial + 16 other protein-coding; 11 of 16 never assayed", "LEAD", r"^\[L04\] heart failure / cardiomyopathy: uninstrumentable 70; locked 21 = 4 microRNA genes \+ 1 mitochondrial genes \+ 16 other protein-coding genes; of those 16, never assayed by any proteomic study: 11"),
 ("Results 4", "all areas: 57 locked = 9 microRNA + 2 mitochondrial + 46 protein-coding; 24 never assayed, 22 assayed", "LEAD", r"^\[L04\] all areas: uninstrumentable 280; locked 57 = 9 microRNA genes \+ 2 mitochondrial genes \+ 46 other protein-coding genes; of those 46, never assayed by any proteomic study: 24 .*assayed without a candidate: 22"),
 ("Results 5; Discussion", "DMD locked (assayed, no candidate)", "LEAD", r"^\[L05\] DMD\s+phase Preclinical; tier locked: assayed, no instrument"),
 ("Results 5; Discussion", "PLN locked (never assayed)", "LEAD", r"^\[L05\] PLN\s+phase Preclinical; tier locked: never assayed"),
 ("Results 6", "liver: no candidate for PCSK9", "S23", r"^\[S23\] PCSK9\s+Liver \(n=208\): no"),
 ("Results 6", "liver: no candidate for ANGPTL3, APOC3, APOB", "S23", r"^\[S23\] APOB\s+Liver \(n=208\): no"),
]
N += [
 ("Abstract; Results 1; Discussion", "10 lipid targets with landmark genetic evidence, all in the set: 6 uninstrumentable, 4 instrumented", "LEAD", r"^\[L07\] lipid targets with landmark human genetic evidence: 10, all in the set; uninstrumentable 6 \(PCSK9, APOC3, ANGPTL3, APOB, NPC1L1, LPA\); instrumented 4 \(HMGCR, CETP, LDLR, LPL\)"),
 ("Results 1", "GTEx whole blood holds a qualifying variant for PCSK9", "LEAD", r"^\[L08\] PCSK9 in GTEx v8 whole blood: tested True; qualifying variant True"),
 ("Results 1", "Table 2 shows a GTEx whole-blood instrument for PCSK9", "T2", r"^\| PCSK9 \|.*\| Yes \|"),
 ("Table 3 footnote", "never assayed: 9 microRNA and 2 mitochondrial overall; 4 and 1 in heart failure", "LEAD", r"^\[L09\] never-assayed locked targets: 35 overall, of which microRNA genes 9 and mitochondrial genes 2; 16 in heart failure / cardiomyopathy, of which microRNA genes 4 and mitochondrial genes 1"),
]
EXTRA_REV = [
 ("Results 8; Table S22", "overlap upper limits 4.0%, 0.8%, 5.1% of the outcome sample", "sample sizes above", "52,363 / 1,320,016 = 3.97% (PCSK9, ANGPTL3); 10,708 / 1,320,016 = 0.81% (APOC3); 52,363 / 1,028,980 = 5.09% (REN)", "upper limits, not measured overlap"),
 ("Table S22", "Fenland effects per 1 SD of normalized plasma abundance", "PMC9904207 main text (prior_art_search/fulltext/PMC9904207.txt)", "phrase '1 s.d. higher normalized plasma abundance'", "the scale of the supplementary table itself is not stated in the workbook"),
 ("Results 3", "11 named targets; PLG and IL1B shown in Table 2 only", "manuscript/display/table2_named_targets.md", "13 rows minus the two rows marked with a dagger or section sign", ""),

 ("Results 7", "11 named launched targets: 7 / 3 / 1 by tier", "manuscript/display/table2_named_targets.md", "the 13 rows of Table 2 without IL1B (other-tissue tier) and PLG (CV-relevant or plasma-protein tier)", "13 rows give 8 / 4 / 1 ([R26]); removing the two gives 7 / 3 / 1"),
 ("Results 3", "a currently marketed drug for 10 of the 11 named targets", "manuscript/display/table2_named_targets.md, column 'Regulatory status (FDA / EMA)'", "count named rows that begin 'Marketed'", "AGTR1, ANGPTL3, APOC3 (EMA only), AVPR2, F2, PCSK9, REN, SLC5A1, ITGA2B, NPC1L1; not APOB"),
 ("Results 8", "REN lead pQTL frequency 0.8%", "manuscript/rescue_controls/exposure_rows.json", "REN_SBP, af_effect", "0.00762883"),
 ("Results 8; Discussion", "APOC3 variant 51 kb from the gene", "manuscript/recovery_map.csv, row APOC3, fenland_bp_from_gene_body", "", "51505"),
 ("Discussion", "heart tissue samples in GTEx about 80 times smaller than eQTLGen", "gtex_v8_tissue_sample_sizes.csv", "31,684 / 386 = 82; 31,684 / 372 = 85", "rounded"),
 ("Results 5", "kidney cortex is the smallest tissue sample", "gtex_v8_tissue_sample_sizes.csv", "minimum eqtl_samples among tissues with eGenes", "73"),
 ("Results 7", "sites of drug action for SLC5A1, ITGA2B, NPC1L1", "manuscript/revision_sources/drug_status.csv, column label_mechanism_text; drug_status/extra_openfda_brand_name_Zetia.json", "FDA label mechanism-of-action text", "quoted in Table 2"),
]
F["R0"] = ("manuscript/recovery_out/recovery_map_summary_before_amendment8.txt",
           "first run of build_recovery_map.py, kept as a record")
EXTRA = [
 ("Methods", "Fenland 10,708 participants", "manuscript/recovery_data/pietzner2021_supplementary_tables.xlsx, sheet ST1", "python -I manuscript_prep/manuscript/recovery_scripts/inspect_xlsx.py <workbook> 5 ST1", "row 'Participants [N]' reads 10708"),
 ("Introduction; Methods", "eQTLGen tested 19,250 genes", "analysis/classify_output.txt : 185", "python manuscript_prep/analysis/classify.py", "the line gives 19,123 autosomal + 127 not in GRCh37.87 = 19,250; the figure is also quoted from the eQTLGen paper in PRIOR_ART.md"),
 ("Introduction; Methods", "eQTLGen n = 31,684; UKB-PPP 54,219; deCODE 35,559", "abstracts of the cited papers in PubMed", r"python manuscript_prep\manuscript\prior_art_search\papers.py lookup \"<title>[ti]\" (abstract text)", "verified against the PubMed abstract on 2026-10-11; not stored in a project file"),
 ("Results 5", "23 of 72 off-target rescues had been measured in a proteomic study", "manuscript/recovery_map.csv", "count rows with reporting_tier = 'candidate: other tissue only' and non-empty pqtl_studies_assayed", "computed on 2026-10-11 from the CSV; the count is not printed in recovery_map_summary.txt"),
 ("Results 5", "9 microRNA genes among the never-assayed locked targets", "manuscript/recovery_map_summary.txt [R30]", "count of names beginning MIR in the [R30] never-assayed list", "counted from the printed list"),
 ("Results 7", "5 protein-only rescues among the 13 are secreted proteins", "analysis/out/per_target_primary.csv, column `secreted`", "rows ANGPTL3, APOC3, PCSK9, PLG, REN", "all five read YES"),
 ("Discussion", "16 of 21 locked heart-failure targets never assayed", "manuscript/recovery_map_summary.txt [R25]", "Heart failure row", "same line as the 70-target row above"),
]
lines = ["# NUMBER_TRACE", "",
         "Every number in `manuscript.md`, with the output file and line that carries it and "
         "the command that regenerates that file. Built by "
         "`python manuscript_prep/manuscript/draft/make_number_trace.py`, which searches each "
         "file for the stated value; a value it cannot find is written **UNVERIFIED**.", "",
         "Paths are relative to `D:\\causal2drug\\manuscript_prep\\`. Run commands from "
         "`D:\\causal2drug`. Numbers marked in the manuscript depend on the "
         "pending manual review of `drug_links/MANUAL_REVIEW.csv`.", "",
         "| # | where in the manuscript | number(s) | file : line | regenerate with |",
         "|---|---|---|---|---|"]
bad = 0
cache = {}
for i, (sec, what, key, pat) in enumerate(N, 1):
    fn, cmd = F[key]
    if fn not in cache:
        cache[fn] = open(os.path.join(MP, fn), encoding="utf-8", errors="replace").read().split("\n")
    hit = [j + 1 for j, l in enumerate(cache[fn]) if re.search(pat, l)]
    if hit:
        loc = "`%s` : %d" % (fn, hit[0])
    else:
        loc = "**UNVERIFIED** (pattern not found in `%s`)" % fn
        bad += 1
    lines.append("| %d | %s | %s | %s | `%s` |" % (i, sec, what, loc, cmd))
lines += ["", "## Numbers not carried by a single output line", "",
          "| where | number | source | how | note |", "|---|---|---|---|---|"]
lines += ["| %s | %s | %s | %s | %s |" % e for e in EXTRA + EXTRA_REV]
lines += ["", "Entries traced: %d; UNVERIFIED: %d." % (len(N), bad)]
open(os.path.join(MS, "NUMBER_TRACE.md"), "w", encoding="utf-8", newline="\n").write(
    "\n".join(lines) + "\n")
print("entries %d, UNVERIFIED %d" % (len(N), bad))
for l in lines:
    if "UNVERIFIED** (" in l:
        print(l[:200])
