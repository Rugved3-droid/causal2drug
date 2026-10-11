# Supplementary bundle — cardiovascular target untestability atlas

## What this is

A row-per-target record of **which cardiovascular drug targets can be studied with human
genetics, and which cannot**. It covers **630 targets** carrying a cardiovascular
indication. For each one it records whether a usable genetic instrument exists in whole
blood, whether one exists in heart, vessel or liver tissue, how strongly the gene is
expressed in each of those tissues, and — where no instrument exists — which
measurement would be needed to create one.

**232 targets (36.8%) cannot be instrumented at all in whole blood** (`untestability_class
= CLASS_A`), including 21 with an approved medicine. **386 (61.3%) can be.**

## What this is not

* **Not a causal claim about any target.** No column states that a target does or does not
  affect disease. The atlas records whether the question can be *asked*, not its answer.
* **Not a colocalization resource.** Colocalization is excluded as an evidence gate
  (Gate D, closed as a documented negative result). Colocalization figures elsewhere in
  this project are *method availability*, not evidence.
* **Not a curation of drug mechanism.** `required_modality` is a declared heuristic over
  expression and secretion. It says where a signal could exist, not that one does.
* **Not a complete Class B count.** `untestability_class = CLASS_B` is assigned at N = 1
  (IL6R). The true count is unknown; `classB_upper_bound.csv` bounds it at 62 of 386
  instrumented targets (16.1%) on cardiovascular traits, 196 (50.8%) on any trait.
* **Not tissue-complete.** The heart, vessel and liver instrument columns are **upper
  bounds**: they use a wider window than pre-registered and are unclumped.
* **Not a clinical decision aid.** Counts are genes and target-disease pairs, never
  patients.

## Files

| file | rows | what it holds |
|---|---|---|
| `cardiovascular_atlas.csv` | 630 | the atlas, 25 columns |
| `cardiovascular_atlas.parquet` | 630 | identical content, typed |
| `classB_upper_bound.csv` | 196 | instrumented targets with a protein-altering association |
| `DATA_DICTIONARY.md` | — | every column: type, units, permitted values, source, derivation |

`DATA_DICTIONARY.md` is written to be sufficient on its own. A reader should be able to
use these files without reading the paper.

## Licence position of incorporated sources

See the table at the end of `DATA_DICTIONARY.md`. Positions are stated as commonly
documented and were **not re-verified against each source's terms page in this build**;
treat them as UNVERIFIED and check before redistributing. The Human Protein Atlas
columns (`tissue_restricted`, `secreted`) carry a share-alike condition under CC BY-SA
3.0, which is the strictest condition in the bundle.

## Regenerating

Exactly two commands, from the repository root, in this order:

```bash
python src/cv_atlas.py            # builds cardiovascular_atlas.csv / .parquet
python src/a1_classB_bound.py     # builds classB_upper_bound.csv
python src/make_supplementary.py  # assembles this bundle and this README
```

`src/cv_atlas.py` and `src/a1_classB_bound.py` read only from `data/raw/`, which is
written once and never modified. `src/make_supplementary.py` copies their frozen outputs
and recomputes every statistic in the dictionary from the shipped files, so the
dictionary cannot drift from the data. The build **fails loudly** if a shipped column has
no dictionary entry, if a dictionary entry names a column that is not shipped, or if the
parquet and CSV disagree on row count or column order.

## Checksums

Recorded at build time in `MANIFEST.json`.
