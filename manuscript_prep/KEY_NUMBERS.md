# KEY_NUMBERS — corrected numbers for the manuscript

> **Final, 2026-10-11.** Drug verification is complete and nothing below is provisional: the primary definition of an approved-drug target is the strict ChEMBL rule plus an FDA or EMA approval record (33 launched targets; section 4). Where a passage below still says "provisional", read "verified".

**Built 2026-10-10 on the approved definitions. Replaces `NUMBERS.md` and `HEADLINE.md`
(Phase 1 records for the original 630-target set).**

> **Headline numbers about approved drugs (§4) are PROVISIONAL.** They wait on your review of
> `drug_links/MANUAL_REVIEW.csv`. Sections 1–3 and 5–9 do not depend on drug verification and
> are final.

Everything is printed by one script; tags **[K01]–[K10]** name its output blocks.

```powershell
cd D:\causal2drug
python manuscript_prep\verify\build_chembl_links.py     # drug links, strict and permissive
python manuscript_prep\analysis\classify.py > manuscript_prep\analysis\classify_output.txt
```

The classifier reads only existing project inputs, writes only under
`manuscript_prep\analysis\out\`, and stops unless it reproduces the shipped state of all 2,457
already-mapped targets and the shipped Human Protein Atlas annotation. Machine-readable copy:
`analysis/out/key_numbers.json`.

## Terms

| term | meaning |
|---|---|
| **Instrumented** | ≥ 1 eQTLGen whole-blood cis-eQTL at p < 5 × 10⁻⁸ within gene body ± 100 kb |
| **No eQTL data** | the gene is not among the 19,250 genes eQTLGen tested. **Lower bound.** Sub-types: |
| — autosomal, not tested | an autosomal gene eQTLGen did not test |
| — X-linked (27) and mitochondrial (2), not analysed | eQTLGen analysed no gene on X, Y or the mitochondrial genome. Relabelled 2026-10-11: the group was previously called "X-linked (29)"; 2 of the 29 are mitochondrial genes (MT-ND1, MT-ND2). No total changes |
| **Autosomal-only lower bound** | autosomal no-eQTL-data targets ÷ targets not on X, Y or MT |
| **No qualifying instrument** | tested, but no variant meets the criterion |
| **Upper bound** | no eQTL data + no qualifying instrument. Reported only as a bound |
| **Launched** | highest Minikel phase over the target's indications *in the same set* |
| **Verified, permissive** | an approved drug with any ChEMBL link to the gene has an approved indication under MeSH C14 or a lipid term |
| **Verified, strict** | as permissive, and the link is direct to this gene, and the indication sits only in cardiovascular MeSH branches (§4) |

## 1. Flow and sets  [K01]

29,476 rows → 29,460 distinct pairs → 29,456 pairs after dropping 4 programmes (2,517 targets,
1,014 indications).

| set | rule | terms | pairs | targets |
|---|---|---|---|---|
| **Primary** | every indication Minikel tags cardiovascular, plus 8 lipid tokens | **93** | **2,128** | **726** |
| S1 | original atlas: MeSH branch C14 + lipid tokens | 44 | 1,543 | 630 |
| S2 | S1 + cardiovascular-tagged terms filed under signs/symptoms | 62 | 1,879 | 663 |
| S3 | MeSH branch C14 only | 37 | 1,364 | 594 |

All three sensitivity sets are in the supplement (Table S3).

## 2. Lower and upper bound, overall  [K04, K10]

| set | targets | instrumented | **no eQTL data (lower)** | autosomal / X-linked + mitochondrial | **autosomal-only lower bound** | no qualifying instrument | **upper bound** | unresolved |
|---|---|---|---|---|---|---|---|---|
| **Primary** | 726 | 443 (61.0 %) | **217 (29.9 %)** | 188 / 27 + 2 | **188 of 697 (27.0 %)** | 63 (8.7 %) | **280 (38.6 %)** | 3 |
| S1 | 630 | 391 (62.1 %) | 180 (28.6 %) | 157 / 21 + 2 | 157 of 607 (25.9 %) | 57 (9.0 %) | 237 (37.6 %) | 2 |
| S2 | 663 | 409 (61.7 %) | 191 (28.8 %) | 167 / 22 + 2 | 167 of 639 (26.1 %) | 60 (9.0 %) | 251 (37.9 %) | 3 |
| S3 | 594 | 373 (62.8 %) | 164 (27.6 %) | 144 / 18 + 2 | 144 of 574 (25.1 %) | 55 (9.3 %) | 219 (36.9 %) | 2 |
| All 2,517 targets | 2,517 | 1,603 (63.7 %) | 663 (26.3 %) | 553 / 107 + 3 | 553 of 2,407 (23.0 %) | 236 (9.4 %) | 899 (35.7 %) | 15 |

Across definitions: lower bound 27.6–29.9 %, autosomal-only 25.1–27.0 %, upper bound
36.9–38.6 %.

Sub-types of the 63 "no qualifying instrument" targets (primary): 36 have no cis-eQTL at
FDR < 0.05; 26 have one below the threshold; 1 (IL1B) has significant variants only outside the
window.

## 3. By highest phase, primary set  [K05]

| phase | targets | instrumented | no eQTL data (lower) | autosomal / X-linked + mitochondrial | autosomal-only lower bound | no qualifying instrument | upper bound |
|---|---|---|---|---|---|---|---|
| Preclinical | 323 | 196 (60.7 %) | 97 (30.0 %) | 85 / 10 + 2 | 85 of 311 (27.3 %) | 29 (9.0 %) | 126 (39.0 %) |
| Phase I | 91 | 53 (58.2 %) | 26 (28.6 %) | 20 / 6 + 0 | 20 of 85 (23.5 %) | 10 (11.0 %) | 36 (39.6 %) |
| Phase II | 171 | 105 (61.4 %) | 57 (33.3 %) | 48 / 9 + 0 | 48 of 162 (29.6 %) | 9 (5.3 %) | 66 (38.6 %) |
| Phase III | 75 | 50 (66.7 %) | 14 (18.7 %) | 14 / 0 + 0 | 14 of 75 (18.7 %) | 11 (14.7 %) | 25 (33.3 %) |
| **Launched** | **66** | 39 (59.1 %) | **23 (34.8 %)** | 21 / 2 + 0 | **21 of 64 (32.8 %)** | 4 (6.1 %) | **27 (40.9 %)** |

Unresolved: 1 Preclinical, 2 Phase I. No phase gradient: χ² = 6.33, df = 4, p = 0.18
(descriptive, not pre-registered). Every set: `tables/tableS3_sensitivity_cv_set.csv`.

## 4. Launched targets with an approved cardiovascular drug — PROVISIONAL  [K06]

### The two rules

| | permissive | **strict** |
|---|---|---|
| Drug is approved (ChEMBL `max_phase` 4) | required | required |
| Link from drug to gene | any ChEMBL target containing the gene | **direct only:** single protein, a specific protein complex, or an mRNA target. Links through a protein family, complex group or selectivity group do not count |
| Approved indication (`max_phase_for_ind` 4) | any MeSH descriptor under C14, or a lipid term | **strictly cardiovascular:** under C14 or a lipid term, **and** in no disease branch other than C14, C23 (signs), C10 (nervous system, for stroke), C08 (respiratory, for pulmonary hypertension); lipid terms may also sit in C18 and C16 |

You asked for a rule in which "a linked approved drug carries an approved cardiovascular
indication". Read literally that is the permissive rule, so I tightened the two words that
were loose in it: *linked* and *cardiovascular*. **The strict rule is my reading of your
instruction; change it if you meant something else.** The exact code is in
`verify/build_chembl_links.py` under "strict rule".

### Counts side by side

| set | launched (Minikel) | rule | verified | **no eQTL data (lower)** | autosomal-only lower | no qualifying instrument | **upper bound** | beyond PCSK9 (lower–upper) |
|---|---|---|---|---|---|---|---|---|
| **Primary** | 66 | **strict** | **35** | **10 (28.6 %)** | 9 of 34 (26.5 %) | 3 (8.6 %) | **13 (37.1 %)** | **9–12** |
| | | permissive | 58 | 20 (34.5 %) | 18 of 56 (32.1 %) | 4 (6.9 %) | 24 (41.4 %) | 19–23 |
| S1 | 52 | strict | 32 | 9 (28.1 %) | 9 of 32 (28.1 %) | 3 (9.4 %) | 12 (37.5 %) | 8–11 |
| | | permissive | 46 | 16 (34.8 %) | 16 of 46 (34.8 %) | 4 (8.7 %) | 20 (43.5 %) | 15–19 |
| S2 | 55 | strict | 35 | 10 (28.6 %) | 9 of 34 (26.5 %) | 3 (8.6 %) | 13 (37.1 %) | 9–12 |
| | | permissive | 49 | 17 (34.7 %) | 16 of 48 (33.3 %) | 4 (8.2 %) | 21 (42.9 %) | 16–20 |
| S3 | 44 | strict | 26 | 6 (23.1 %) | 6 of 26 (23.1 %) | 2 (7.7 %) | 8 (30.8 %) | 5–7 |
| | | permissive | 40 | 13 (32.5 %) | 13 of 40 (32.5 %) | 3 (7.5 %) | 16 (40.0 %) | 12–15 |

### Applying the Methods definition of "approved" (added 2026-10-11)

"Approved" = ever approved by the FDA or the EMA for a cardiovascular indication: the strict
rule, and at least one strict-rule drug has an approval record in Drugs@FDA or the EMA
medicines table (`manuscript/APPROVED_35_CHECK.csv`; `classify.py` [K06], rows "APPROVED").
Two of the 35 drop out: **MTOR** (perhexiline; instrumented) and **PLG** (anistreplase,
streptokinase; no eQTL data).

| set | strict | **approved** | no eQTL data | autosomal-only | no qualifying instrument | upper bound |
|---|---|---|---|---|---|---|
| **Primary** | 35 | **33** | **9 (27.3 %)** | 8 of 32 (25.0 %) | 3 (9.1 %) | **12 (36.4 %)** |
| S1 | 32 | 30 | 8 (26.7 %) | 8 of 30 (26.7 %) | 3 (10.0 %) | 11 (36.7 %) |
| S2 | 35 | 33 | 9 (27.3 %) | 8 of 32 (25.0 %) | 3 (9.1 %) | 12 (36.4 %) |
| S3 | 26 | 24 | 5 (20.8 %) | 5 of 24 (20.8 %) | 2 (8.3 %) | 7 (29.2 %) |
| Primary, Phase III or launched | 41 | 39 | 12 (30.8 %) | 11 of 38 (28.9 %) | 4 (10.3 %) | 16 (41.0 %) |

The manuscript text uses the "approved" row. The strict and permissive rows below are kept
for the tables.

### Names, primary set

| | strict | permissive adds |
|---|---|---|
| **No eQTL data** | AGTR1, ANGPTL3, APOB, APOC3, AVPR2 (X-linked), F2, PCSK9, PLG, REN, SLC5A1 | ADRA1A, ADRA1B, ADRA1D, ADRA2B, ADRA2C, EDNRA, EDNRB, GHR, GLA (X-linked), PDE4C |
| **No qualifying instrument** | IL1B, ITGA2B, NPC1L1 | GUCY1B1 |

### Where the two rules disagree

23 launched targets in the primary set are verified by the permissive rule only:

| reason the strict rule rejects them | targets |
|---|---|
| Every approved cardiovascular drug is linked through a family or group target | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, EDNRA, EDNRB, CACNA1C, GUCY1A1, GUCY1B1, PDE4A, PDE4B, PDE4C, PDE4D, TUBB (16) |
| The approved indication also sits outside cardiovascular MeSH branches | GHR, GLA, KLKB1, SERPING1, C5AR1, EPAS1, VEGFA (7) |

The exact reason per target is the `strict_fail_reason` column of
`drug_links/chembl_target_summary.csv`.

**Read this before deciding.** The first group is where the strict rule is probably too
strict. Prazosin and doxazosin do act on α1-adrenoceptors, bosentan on both endothelin
receptors, and nitrates and riociguat on soluble guanylate cyclase. ChEMBL records those drugs
against the receptor family or the enzyme complex group because they are not subtype-selective,
so the strict rule cannot attribute them to one gene. Whether a non-selective drug "verifies"
each subtype gene is the judgement the review file asks you to make. The second group is where
the strict rule does what you asked: GHR and GLA drop out because Turner syndrome and Fabry
disease are not cardiovascular indications.

`drug_links/MANUAL_REVIEW.csv` holds every piece of permissive evidence for the targets where
the rules disagree: **327 rows, 32 targets** (23 launched, the rest Phase III). Each row has the
drug, the link type, the indication, why it fails, and two empty columns for your decision.

Not verified by either rule (primary, 8): AGT, EGFR, ITGB2, LPL, PPARG, RXFP1, SLC12A2, TNF.

Phase III in Minikel with an approved cardiovascular drug in ChEMBL: strict 6 (ADORA2A, AVPR1A,
DPP4, HCN4, PDE3A, SLC6A2); permissive 15. Labels were not edited. Counting those as well:
strict 41 targets, 13 no eQTL data (31.7 %), 17 upper bound (41.5 %); permissive 73, 27
(37.0 %), 33 (45.2 %). All label disagreements: `tables/tableS6_phase_label_vs_chembl.csv`.

### Targets that may be named as examples

Only targets that are strict-verified **and** hold the same phase in S2, so nothing that enters
through a non-cardiac addition. Column `example_eligible` in Table S1 and S7.

| | eligible |
|---|---|
| Launched, no eQTL data (10) | AGTR1, ANGPTL3, APOB, APOC3, AVPR2, F2, PCSK9, PLG, REN, SLC5A1 |
| Launched, no qualifying instrument (3) | IL1B, ITGA2B, NPC1L1 |
| Phase III label with an approved drug, uninstrumentable (4) | AVPR1A, HCN4, SLC6A2, PDE3A |

GHR, GLA, PDE4C and the adrenoceptor and endothelin subtypes are not eligible and appear in no
figure. Launched instrumented targets eligible for contrast include ACE, HMGCR, ITGB3, ADRB1.

## 5. Cis-window sensitivity  [K07]

| window | primary, all 726: upper bound | primary, 66 launched: upper bound |
|---|---|---|
| ± 100 kb (primary) | 280 (38.6 %) | 27 (40.9 %) |
| ± 500 kb | 279 (38.4 %) | 26 (39.4 %) |
| ± 1 Mb | 279 (38.4 %) | 26 (39.4 %) |

One target changes: IL1B (nearest significant variant 119,995 bp from the gene body). The lower
bound and the autosomal-only lower bound cannot change with the window. eQTLGen tested variants
within 1 Mb of the gene centre, so ± 1 Mb from the gene body is not fully covered (largest
distance observed: 812,919 bp).

## 6. Why targets have no eQTL data  [K08, K10]

* **Autosomal, not tested (188):** not expressed in blood. Of the 182 with a GTEx value, 180
  have whole-blood median TPM < 1.
* **X-linked (27) and mitochondrial (2), not analysed:** a property of the resource. None of
  eQTLGen's 19,250 tested genes is on X, Y or MT. 12 of the 27 X-linked targets and both
  mitochondrial targets are expressed in blood at TPM ≥ 1 (`tables/tableS10_no_eqtl_data_by_chromosome.csv`).
* **Tissue restriction**, no eQTL data versus instrumented: OR 7.47 (95 % CI 4.12–13.55),
  p = 1.0 × 10⁻¹⁵; autosomal only, OR 12.18 (5.58–26.62). S1 10.13, S2 8.51, all targets 7.01.
* **No qualifying instrument is a different group:** OR 1.39 (0.77–2.51), p = 0.32; median blood
  TPM 0.54, against 0.016 for no eQTL data and 5.44 for instrumented.

Instrument yield in other tissues, primary set, upper bounds: heart left ventricle 6.6 %, atrial
appendage 6.5 %, aorta 9.1 %, liver 4.4 %, any of the four 17.2 %. Of the 217 no-eQTL-data
targets, 28 have a genome-wide significant top cis-eQTL in one of those tissues.

## 7. Outcome-GWAS file attrition  [K09]

Of 2,128 primary pairs, 771 (22 indications) match a European-only harmonised GWAS Catalog
trait; for 404 (52.4 %) the file cannot be range-queried; 8 were never probed, so this is a
floor.

## 8. Unchanged from Phase 1

`python manuscript_prep\verify\verify_numbers.py`

| | value | block |
|---|---|---|
| HMGCR → LDL-C | MR β +0.389 (SE 0.012), p = 4.2 × 10⁻²²⁹, F = 224; PASS | [V08] |
| IL6R → CRP | MR β +0.376 (SE 0.023), p = 1.2 × 10⁻⁵⁹, F = 134.5; expected negative; FAIL on direction | [V08] |
| PCSK9 | 0 of 127,341,798 rows of the full eQTLGen file | `verify/pcsk9_scan.txt` |
| Class C test | 50 paired loci; PP4 higher with eQTLGen at 32; Wilcoxon p = 0.0262; refuted | [V09] |

## 9. Gene records  [K02]

60 symbols had no unique GRCh37 record; the HGNC rule resolves 45. Unresolved in the primary
set: **MIR29B1, MIR92A1, MIR195**.

## 10. Open conflicts

| | conflict | state |
|---|---|---|
| N1 | "Strict" as you worded it equals the permissive rule | I tightened link and indication; my reading, stated in §4 |
| N2 | The strict rule drops well-known targets (α-adrenoceptors, endothelin receptors, sGC) because ChEMBL links their drugs at family level | in `MANUAL_REVIEW.csv` for your decision |
| N3 | ChEMBL indication rows include weak entries | not filtered; visible in the review file |
| N4 | IL1B is strict-verified through rilonacept for pericarditis, while its Minikel Launched label is for atherosclerosis | kept; column `same_indication_as_minikel_launched` in the summary file shows which targets match on the indication itself |
| N5 | `deliverables/cv/cv_robustness.json` still says 20 | left as instructed; marked superseded in `README.md` |
| N6 | One passing positive control (HMGCR) | unchanged; in LIMITATIONS |
