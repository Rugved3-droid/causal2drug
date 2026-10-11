# RECOVERY_RULE — pre-specified rule for the recovery map

**Written 2026-10-10, before any GTEx or pQTL result for any target was downloaded, opened or
looked at.** The only things known at the time of writing are which resources were named by
the project owner and how those resources are generally organised. This file is not edited
after results are seen; anything that has to change is added under "Amendments" at the end,
dated, with the reason.

## Question

For each cardiovascular target that has no instrument in eQTLGen whole blood, does a cis
instrument exist in another tissue (GTEx v8) or at the protein level (published plasma
cis-pQTL tables)?

## Scope

Every target in the primary set (726 targets) whose status is "no eQTL data" or "no qualifying
instrument" in `manuscript_prep/analysis/out/per_target_primary.csv`, read as it stands today.
The 13 named targets are part of that list and are also reported separately. Targets with an
unresolved gene record (MIR29B1, MIR92A1, MIR195) are out of scope and listed as such.

## Definition of an instrument

The same as the main analysis: **a variant with p < 5 × 10⁻⁸ located within the gene body
± 100,000 bp.** One qualifying variant is enough.

- Gene body coordinates: the start and end of the gene in the annotation and genome build of
  the resource being checked (GTEx v8: GENCODE v26, GRCh38, as given in GTEx's own files). No
  coordinate is lifted between builds to make a variant qualify.
- Genes are matched by Ensembl gene ID with the version suffix removed. If a target's Ensembl
  ID is absent from a resource, the symbol is **not** used as a fallback for GTEx. For pQTL
  tables, which are organised by protein, matching is by the gene symbol or UniProt ID the
  table itself gives; the match used is recorded.
- A resource's own significance flag is necessary but not sufficient. The p < 5 × 10⁻⁸ and
  ± 100 kb conditions are applied on top of it.

## Resources, in this order

### (a) GTEx v8 single-tissue cis-eQTL, all tissues

- File set: the GTEx v8 significant variant–gene pairs, one file per tissue
  (`*.v8.signif_variant_gene_pairs.txt.gz`), with gene coordinates from the matching
  `*.v8.egenes.txt.gz` files.
- A target has a GTEx instrument in a tissue if at least one significant pair for its gene in
  that tissue has `pval_nominal` < 5 × 10⁻⁸ and a variant position within gene body ± 100 kb.
- Recorded per target: every tissue with an instrument; for each, the lead qualifying variant
  (lowest p) and its p-value.
- **CV-relevant tissues**, flagged: Heart – Left Ventricle; Heart – Atrial Appendage; Artery –
  Aorta; Artery – Coronary; Artery – Tibial; Liver; Kidney – Cortex; Adrenal Gland. A target
  is "CV-tissue fixable" if any of these eight has an instrument.
- Also recorded: whether the gene was tested in GTEx at all in each tissue (present in the
  egenes file), so that "not expressed or not tested" can be told apart from "tested, no
  instrument".
- Whole Blood in GTEx is reported but does **not** count toward "fixable-tissue". The question
  is whether a tissue other than blood rescues the target.

### (b) Published plasma cis-pQTL tables

- UK Biobank Pharma Proteomics Project: Sun et al., *Nature* 2023 (10.1038/s41586-023-06592-6),
  the published supplementary table of pQTLs.
- deCODE: Ferkingstad et al., *Nat Genet* 2021 (10.1038/s41588-021-00978-w), the published
  supplementary table of pQTLs.
- Fenland: Pietzner et al., *Science* 2021 (10.1126/science.abj1541), the published
  supplementary table of cis-pQTLs.
- A target has a pQTL instrument in a study if the table lists a pQTL for its protein that the
  study labels **cis**, with p < 5 × 10⁻⁸, **and** whose lead variant lies within gene body
  ± 100 kb in the build the table uses.
- If the table gives no variant position, or the build cannot be established, the distance
  cannot be checked. The cell is then recorded as "cis-pQTL reported; distance UNVERIFIED" and
  does **not** count as an instrument.
- Recorded per target and study: assayed or not assayed (where the table or its companion
  protein list allows this to be read); lead cis variant; p-value; distance from the gene body.
- Only the published tables are used. Full summary statistics are not downloaded, so a
  qualifying variant that is not a reported lead is not seen. This makes the pQTL check
  conservative, and it is stated wherever the result is reported.

## Categories

Assigned mechanically from the cells above.

| category | condition |
|---|---|
| **fixable-tissue** | instrument in ≥ 1 GTEx tissue other than Whole Blood; none in any pQTL table |
| **fixable-protein** | instrument in ≥ 1 pQTL table; none in any GTEx non-blood tissue |
| **fixable-both** | both of the above |
| **locked** | no instrument in any checked resource |

- **"Locked" means locked relative to the resources checked here** (GTEx v8 non-blood tissues;
  UKB-PPP, deCODE and Fenland published cis-pQTL tables). It is never a statement that no
  instrument exists. This qualification accompanies the word wherever it appears.
- X-linked targets keep an **X flag** in addition to their category.
- A target with an instrument only in GTEx Whole Blood is "locked" under these categories and
  is flagged "GTEx whole blood only", because it shows a blood instrument exists outside
  eQTLGen.
- If any cell needed to decide a target's category is UNVERIFIED, the category is given as
  computed from the verified cells and flagged "category depends on UNVERIFIED cell".

## Summaries, fixed in advance

Counts of the four categories, with the X flag:

1. overall;
2. by eQTLGen status (no eQTL data: autosomal / X-linked; no qualifying instrument);
3. by highest clinical phase in the primary set;
4. by indication area (the ten areas of Table 3; a target counted in each area it appears in);
5. for the 13 named targets, one row each;
6. fixable through a CV-relevant tissue versus only through another tissue.

No significance test is planned. These are descriptive counts.

## Data access

Before any classification, the resources that can actually be downloaded and read are listed.
**If any of the four resources cannot be obtained, work stops and the owner is told. No other
resource is substituted.**

## What this rule does not do

- It does not test whether a recovered instrument gives a valid or concordant MR estimate.
- It does not use trans-QTLs, splicing QTLs, rare variants or any resource not listed.
- It does not relax the threshold or the window for any target.

## Amendments

None at the time of writing (2026-10-10).

### 2026-10-11 — clarifications fixed before classification

Logged after reading only sheet names, column headers and the first data rows of each workbook
(none of which is a target in scope), and the column headers of one GTEx tissue file. No
result for any in-scope target had been seen. None of these changes the threshold, the window
or the categories.

1. **deCODE gene body.** The deCODE table gives the transcription start site, not the gene
   end. The gene body is taken from GENCODE v26 (GRCh38), the same build as the table, as
   carried in the GTEx v8 `egenes` files (matched by Ensembl gene ID). If a gene is in no GTEx
   `egenes` file, the distance is UNVERIFIED and the cell does not count.
2. **deCODE p-value.** The table has adjusted and unadjusted −log10(P). The adjusted column is
   used (the study's own reported value; the more conservative of the two).
3. **UKB-PPP sheets.** The workbook has a discovery-cohort sheet (ST9) and a combined-cohort
   sheet (ST10). Both are checked; a qualifying cis variant in either counts, and the sheet is
   recorded. Gene body: the gene start and end in the study's own protein table (ST3, GRCh38).
4. **Fenland build.** The Fenland table does not state its genome build. It is established by
   comparing the gene start/end in the study's own protein table (ST3) with Ensembl GRCh37
   release 87 and with GENCODE v26 (GRCh38); the comparison is printed by the script. Gene
   body: ST3 start and end. If neither build matches, every Fenland distance is UNVERIFIED.
5. **Matching to pQTL tables.** A target matches a protein if the table's own gene symbol
   equals the target symbol, or the table's own Ensembl gene ID equals the target's Ensembl ID
   (UKB-PPP's protein table has no Ensembl ID, so symbol only). The key that matched is
   recorded. Multi-gene entries (protein complexes) are split on the table's separator.
6. **"Lead variant".** deCODE and Fenland list secondary (conditionally independent) signals
   as well as the regional sentinel. Any variant listed for a signal the study labels cis is
   eligible; the one with the lowest p-value that meets both conditions is recorded.
7. **X chromosome.** Each resource is used as published. Whether it covers chromosome X is
   reported, not assumed.

### 2026-10-11 (later) — amendments made AFTER the first results were seen

These two were decided by the project owner after the first run
(`recovery_out/recovery_map_summary_before_amendment8.txt`,
`recovery_out/recovery_map_before_amendment8.csv`). They are post-hoc and are reported as such.

8. **Symbol match in GTEx after a documented Ensembl ID change (post-hoc).** The original rule
   forbids a symbol fallback for GTEx. Exception: a target that is in no GTEx `egenes` file
   under its GRCh37 Ensembl ID is matched through the Ensembl ID that HGNC currently records
   for the same symbol, if and only if (i) that ID differs from the GRCh37 ID (a documented ID
   change), (ii) GTEx v8 carries that ID under the same gene symbol, (iii) the chromosome is
   the same in the GRCh37 record, the GTEx record and the current Ensembl GRCh38 record, and
   (iv) at least 95% of the GTEx gene body lies inside the current Ensembl GRCh38 record of
   that ID. Strand and gene length versus GRCh37 are reported but do not gate the match,
   because coordinates cannot be compared across builds without lifting them.
   Applied to three targets (`recovery_out/gtex_id_change_resolution.csv`,
   `recovery_map_summary.txt` [R02b], [R31]):

   | target | GRCh37 ID → ID in GTEx v8 | check | category before → after |
   |---|---|---|---|
   | ADRA2B | ENSG00000222040 → ENSG00000274286 | chr 2 in all three; 99.97% of GTEx body inside the GRCh38 record; same strand | locked → fixable-tissue (18 non-blood tissues, including aorta and left ventricle) |
   | GDF2 | ENSG00000128802 → ENSG00000263761 | chr 10 in all three; 100% inside; **strand differs between GRCh37 (−) and GRCh38 (+)** | locked → locked |
   | MIR145 | ENSG00000269936 → ENSG00000276365 | chr 5 in all three; 100% inside; **gene model differs: 2,549 bp in GRCh37.87, 88 bp in GTEx** | locked → locked |

   HGNC also documents an ID change for MIR122, MIR143 and MIR21, but the new IDs are not in
   GTEx v8, so nothing changes for them. Net effect: locked 58 → 57; fixable 222 → 223.

9. **Reporting tiers (post-hoc; no change to the instrument rule or to the four categories).**
   - **plausible rescue**: an instrument in any of the eight pre-specified CV-relevant GTEx
     tissues, or a plasma cis-pQTL instrument in any of the three studies;
   - **off-target-tissue rescue**: an instrument only in GTEx tissues outside the eight
     (and no pQTL instrument);
   - **locked: never assayed**: locked, and measured by none of the three proteomic studies;
   - **locked: assayed, no instrument**: locked, and measured by at least one.
   "Off-target" refers to the pre-specified tissue list, not to a per-drug judgement about
   the site of action.

### 2026-10-11 (after a second-implementation audit) — amendment 10 and correction 11

A second implementation of this rule, written in a separate AI-assisted session from the rule
documents alone, without access to the original code or per-target outputs (`audit/AUDIT.md`),
reproduced every category and tier and found one ambiguity and two errors in the original
code. None changes a category or a tier.

10. **deCODE "Included in analysis" (ambiguity resolved; post-hoc).** The deCODE aptamer table
    has a column "Included in analysis". The rule says "assayed or not assayed" without saying
    which column decides. **Chosen reading: an aptamer marked "No" is not counted as assayed**,
    because the study did not analyse it and so could not have reported a pQTL for it. This is
    how the first run was coded; it is now stated. Four in-scope targets have only excluded
    aptamers (ALPP, CDK5, FGA, ITGA2B). Under the other reading (any aptamer listed counts)
    their deCODE cell would read "assayed", and **no tier changes under either reading**: all
    four are assayed by another study, and their excluded aptamers have no cis rows
    (`recovery_map_summary.txt` [R32]).

11. **Corrections to the code (errors, not rule changes).**
    - *P-value underflow.* Candidates carried p as a floating-point number and every p below
      10^-300 was stored as 0. Such values tied, and the tie was broken by sheet or row order,
      so the variant recorded as "lead" was not always the one with the lowest p, contrary to
      clarification 6. Candidates are now compared on -log10(p), from the table's own log-p
      column (UKB-PPP, deCODE) or from effect/SE where the table prints "<1E-300" (Fenland).
      Effect: the recorded UKB-PPP sheet changes for 30 targets and the recorded lead variant
      for 9 (ANGPTL3, CDH17, F11, F9, GHR, GLA, ITGB6, SCG3, TYRO3); the deCODE lead changes
      for 2 (MBL2, SERPINA4); no Fenland lead changes. No instrument yes/no, category or tier
      changes (`recovery_out/compare_with_audit_output.txt` [X03]-[X07]).
    - *Fenland assay coverage.* A protein was matched to the Fenland table through its Ensembl
      ID only on cis rows. CCN1 (listed as CYR61) and GRK2 (listed as ADRBK1) appear only on
      trans rows and were recorded as "not assayed". The match now uses every row. Both become
      "assayed; no cis-pQTL in table"; neither changes tier (CCN1 has a candidate elsewhere;
      GRK2 was already "locked: assayed" through deCODE).
    - The map before these corrections is kept as
      `recovery_out/recovery_map_before_audit_fix.csv`.

### 2026-10-11 — naming only (no rule change)

The category and tier *names* written to `recovery_map.csv` were changed to match the
manuscript, which describes candidate instrument availability and not rescue. The definitions
above are untouched. Old name -> new name: fixable-tissue -> candidate-tissue; fixable-protein ->
candidate-protein; fixable-both -> candidate-both; plausible rescue -> candidate: CV-relevant
tissue or plasma protein; off-target-tissue rescue -> candidate: other tissue only. "Locked",
"locked: never assayed" and "locked: assayed, no instrument" are unchanged. The text above keeps
the names it was written with.
