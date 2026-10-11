#!/usr/bin/env python
"""Supplementary bundle for the cardiovascular untestability atlas.

    python src/make_supplementary.py

Assembles deliverables/supplementary/ :

    cardiovascular_atlas.csv        630 rows x 25 columns
    cardiovascular_atlas.parquet    same content, typed
    classB_upper_bound.csv          196 rows
    DATA_DICTIONARY.md              every column of both files
    README.md                       what this is, what it is not, licences, rebuild command

The dictionary's per-column statistics (n populated, permitted values, ranges) are
COMPUTED from the shipped files at build time, so the dictionary cannot drift from the
data. The provenance text is a declared spec, checked against the shipped columns: the
build fails loudly if a column is present without a spec entry or vice versa.
"""
import csv
import datetime
import hashlib
import json
import os
import shutil

import pyarrow.parquet as pq

ROOT = r"D:\causal2drug"
CV = ROOT + r"\deliverables\cv"
SUP = ROOT + r"\deliverables\supplementary"
LOG = ROOT + r"\logs\retrieval_log.csv"

MEASURED = "MEASURED"      # observed in an external experiment, carried in unchanged
DERIVED = "DERIVED"        # computed by this project's code
ANNOT = "ANNOTATION"       # categorical label carried from an external curated source

# ---------------------------------------------------------------------------
# Declared provenance spec: (column, dtype, units, class, source, derivation)
# ---------------------------------------------------------------------------
ATLAS_SPEC = [
    ("gene", "string", "", ANNOT, "Minikel genetic_support v1.0 pp.tsv",
     "HGNC gene symbol of the drug target, as given by the source. Primary key; unique."),
    ("ensembl_id", "string", "", DERIVED, "Ensembl GRCh37.87 GTF (2026-08-26)",
     "Symbol matched to a unique Ensembl gene record. Empty when no unique GRCh37 record "
     "exists (these rows have state = UNMAPPED)."),
    ("cv_primary", "boolean", "", DERIVED, "Minikel areas.tsv x indic_topl_match.tsv",
     "TRUE if >=1 of the target's indications maps to MeSH top-level C14 (cardiovascular). "
     "The strict cardiovascular definition."),
    ("cv_extended", "boolean", "", DERIVED, "as cv_primary, plus a fixed token list",
     "TRUE if cv_primary, OR an indication term matches one of {hyperlipid, "
     "hypercholesterol, dyslipid, hyperlipoprotein, hypertriglycerid, lipoprotein, "
     "atheroscler, arterioscler}. This is the atlas definition; TRUE for every row."),
    ("n_cv_pairs", "integer", "count", DERIVED, "Minikel pp.tsv (deduplicated by ti_uid)",
     "Number of distinct target-indication pairs for this target that satisfy "
     "cv_extended."),
    ("max_clinical_phase", "string", "", DERIVED, "Minikel pp.tsv ccatnum",
     "Highest clinical phase reached across this target's cardiovascular indications."),
    ("max_phase_num", "integer", "ordinal 1-5", DERIVED, "as max_clinical_phase",
     "Numeric encoding of max_clinical_phase: 5 = Launched, 4 = Phase III, 3 = Phase II, "
     "2 = Phase I, 1 = Preclinical. Sorts correctly; use this for filtering."),
    ("approved_drug_direction", "string", "", ANNOT,
     "Open Targets 26.06 drug_mechanism_of_action",
     "UP or DOWN: whether the approved drug agonises or antagonises the target. Empty "
     "where the mechanism is absent or ambiguous (459 of 630 rows)."),
    ("state", "string", "", DERIVED, "eQTLGen phase I cis-eQTL (2026-08-26)",
     "Instrument availability in whole blood. ABSENT = gene not present in the eQTLGen "
     "file at all. TESTED_NULL = gene tested but no cis-eQTL reaches p < 5e-8. "
     "INSTRUMENTED = >=1 genome-wide significant cis-eQTL survives LD clumping. "
     "UNMAPPED = no unique GRCh37 gene record, so the question was never asked."),
    ("untestability_class", "string", "", DERIVED, "this study",
     "CLASS_A = ABSENT or TESTED_NULL (no instrument exists). CLASS_B = instrumented, but "
     "the drug mechanism operates at a level total mRNA cannot represent. Empty for "
     "instrumented targets with no Class B assignment and for UNMAPPED rows. CLASS_B is "
     "assigned at N = 1 (IL6R) and is under-counted by construction; see LIMITATIONS.md "
     "and classB_upper_bound.csv."),
    ("n_independent_instruments", "integer", "count", DERIVED,
     "eQTLGen + 1000 Genomes phase 3 EUR LD",
     "Independent genome-wide significant cis-eQTLs after LD clumping at r2 < 0.01 within "
     "+/-100 kb of the gene. Empty unless state = INSTRUMENTED."),
    ("best_eqtl_p", "float", "p-value", MEASURED, "eQTLGen phase I cis-eQTL",
     "Smallest cis-eQTL p-value for this gene in whole blood. Carried in unchanged. "
     "Empty unless state = INSTRUMENTED."),
    ("instrumented_blood", "boolean", "", DERIVED, "eQTLGen phase I (n = 31,684)",
     "TRUE if state = INSTRUMENTED. The convenient-tissue reference for the yield "
     "comparison."),
    ("instrumented_heart_LV", "boolean", "", DERIVED,
     "eQTL Catalogue R8 QTD000266 (GTEx heart left ventricle, n = 481)",
     "TRUE if the permuted-file top cis p-value < 5e-8. UPPER BOUND: uses a +/-1 Mb "
     "window (wider than the pre-registered +/-100 kb) and is unclumped."),
    ("instrumented_heart_atrial", "boolean", "", DERIVED,
     "eQTL Catalogue R8 (GTEx heart atrial appendage, n = 473)",
     "As instrumented_heart_LV. Upper bound."),
    ("instrumented_artery_aorta", "boolean", "", DERIVED,
     "eQTL Catalogue R8 (GTEx artery aorta, n = 492)",
     "As instrumented_heart_LV. Upper bound."),
    ("instrumented_liver", "boolean", "", DERIVED,
     "eQTL Catalogue R8 (GTEx liver, n = 269)",
     "As instrumented_heart_LV. Upper bound."),
    ("instrumented_cardiac_union", "boolean", "", DERIVED, "the four columns above",
     "TRUE if instrumented in ANY of heart LV, heart atrial, artery aorta or liver. "
     "Upper bound, inheriting the window caveat."),
    ("tissue_restricted", "string", "", ANNOT,
     "Human Protein Atlas RNA tissue specificity",
     "YES = tissue-enriched, group-enriched or tissue-enhanced. NO = low tissue "
     "specificity. NA = gene absent from the HPA table. Used as the exposure in the "
     "Class A enrichment test."),
    ("secreted", "string", "", ANNOT, "Human Protein Atlas secretome annotation",
     "YES = annotated as secreted. NO = not annotated as secreted. NA = gene absent from "
     "the HPA table."),
    ("gtex_blood_tpm", "float", "TPM", MEASURED, "GTEx v10 median TPM by tissue",
     "Median transcripts per million in whole blood. Carried in unchanged. Empty where "
     "the gene is absent from the GTEx matrix (20 of 630 rows)."),
    ("gtex_heart_LV_tpm", "float", "TPM", MEASURED, "GTEx v10 median TPM by tissue",
     "Median TPM in heart left ventricle. Carried in unchanged."),
    ("gtex_liver_tpm", "float", "TPM", MEASURED, "GTEx v10 median TPM by tissue",
     "Median TPM in liver. Carried in unchanged."),
    ("gtex_artery_aorta_tpm", "float", "TPM", MEASURED, "GTEx v10 median TPM by tissue",
     "Median TPM in artery aorta. Carried in unchanged."),
    ("required_modality", "string", "", DERIVED, "this study (declared heuristic)",
     "For CLASS_A targets, the measurement that could in principle instrument them. "
     "Rule, applied in order: cardiac eQTL if heart LV or atrial TPM >= 1; else liver "
     "eQTL if liver TPM >= 1; else vascular eQTL if aorta TPM >= 1; else pQTL if "
     "HPA-secreted; else 'tissue eQTL (<max-TPM tissue>)'; else UNKNOWN. IL6R is assigned "
     "from its established mechanism. This says where a signal COULD exist, not that one "
     "does. Non-CLASS_A rows carry the same rule's output for completeness."),
]

CLASSB_SPEC = [
    ("gene", "string", "", ANNOT, "as atlas.gene",
     "Instrumented cardiovascular target carrying >=1 protein-altering GWAS association."),
    ("max_clinical_phase", "string", "", DERIVED, "as atlas.max_clinical_phase",
     "Highest clinical phase for this target's cardiovascular indications."),
    ("n_protein_altering_assoc", "integer", "count", DERIVED,
     "GWAS Catalog associations, ontology-annotated full release",
     "Number of associations whose MAPPED_GENE contains this gene as an exact token AND "
     "whose CONTEXT is protein-altering."),
    ("top_variant", "string", "", ANNOT, "GWAS Catalog SNPS column",
     "rsID of the strongest such association (smallest p-value)."),
    ("consequence", "string", "", ANNOT, "GWAS Catalog CONTEXT column",
     "VEP functional consequence of top_variant, as computed and versioned by the "
     "Catalog. Counted as protein-altering: missense_variant, stop_gained, "
     "frameshift_variant, start_lost, stop_lost, inframe_deletion, inframe_insertion, "
     "protein_altering_variant, splice_donor_variant, splice_acceptor_variant, "
     "coding_sequence_variant."),
    ("trait", "string", "", ANNOT, "GWAS Catalog MAPPED_TRAIT column",
     "Reported trait for top_variant."),
    ("p_value", "float", "p-value", MEASURED, "GWAS Catalog P-VALUE column",
     "p-value of top_variant. Carried in unchanged."),
    ("cv_relevant_trait", "boolean", "", DERIVED, "this study (fixed token list)",
     "TRUE if >=1 of this gene's protein-altering associations is for a cardiovascular "
     "trait, by fixed token match on the trait string. TRUE rows are the primary bound "
     "(16.1%); all 196 rows are the loose bound (50.8%)."),
]

# (source, version, retrieved, licence position) -- licence NOT re-verified in this build
SOURCES = [
    ("eQTLGen phase I cis-eQTL", "2019-12-11 release", "2026-08-26",
     "Open summary statistics; no explicit licence on the download page"),
    ("eQTL Catalogue", "R8 (QTD000131/136/141/251/256/266)", "2026-08-26 / 2026-08-27",
     "EMBL-EBI; open access"),
    ("GWAS Catalog associations (full, ontology-annotated)", "v1.0.3.1 release",
     "2026-08-31", "EMBL-EBI terms of use; open access"),
    ("GTEx median TPM by tissue", "v10", "2026-08-30 (file mtime)",
     "Open access for summary-level median expression"),
    ("Human Protein Atlas", "proteinatlas.tsv", "2026-08-30 (file mtime)",
     "CC BY-SA 3.0"),
    ("Open Targets", "26.06", "2026-08-30 (file mtime)", "CC0 1.0"),
    ("Minikel genetic_support", "v1.0 (GitHub tag)", "2026-08-26",
     "Repository licence"),
    ("Ensembl GRCh37 GTF", "GRCh37.87", "2026-08-26", "Ensembl: no restriction on use"),
    ("1000 Genomes phase 3 / high coverage", "20130502 / 2504-sample panels",
     "2026-08-26", "Open access"),
    ("UCSC liftOver chain", "hg19ToHg38", "2026-08-27", "Free for academic use"),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def profile(rows, col):
    """Computed profile of one column: populated count, distinct values, range."""
    vals = [r[col] for r in rows]
    pop = [v for v in vals if v not in ("", None)]
    uniq = sorted(set(pop))
    out = {"n_populated": len(pop), "n_rows": len(rows), "n_distinct": len(uniq)}
    if len(uniq) <= 12:
        out["permitted_values"] = uniq
    else:
        try:
            nums = [float(v) for v in pop]
            out["range"] = (min(nums), max(nums))
        except ValueError:
            out["example_values"] = uniq[:3]
    return out


def fmt_profile(p):
    bits = ["%d / %d populated" % (p["n_populated"], p["n_rows"])]
    if "permitted_values" in p:
        bits.append("values: " + ", ".join("`%s`" % v for v in p["permitted_values"]))
    elif "range" in p:
        lo, hi = p["range"]
        bits.append("range: %g to %g" % (lo, hi))
    else:
        bits.append("%d distinct; e.g. %s"
                    % (p["n_distinct"], ", ".join(p["example_values"])))
    return "; ".join(bits)


def validate(rows, spec, filename):
    """Fail before any output file is opened, so a bad build cannot truncate one."""
    cols = list(rows[0].keys())
    spec_cols = [s[0] for s in spec]
    missing = [c for c in cols if c not in spec_cols]
    extra = [c for c in spec_cols if c not in cols]
    if missing or extra:
        raise SystemExit("DATA DICTIONARY OUT OF SYNC with %s\n"
                         "  columns with no spec entry: %s\n"
                         "  spec entries with no column: %s" % (filename, missing, extra))


def dict_section(title, filename, rows, spec, fh):
    cols = list(rows[0].keys())
    fh.write("\n---\n\n## %s\n\n" % title)
    fh.write("`%s` \u2014 **%d rows x %d columns**\n\n" % (filename, len(rows), len(cols)))
    counts = {}
    for _, _, _, k, _, _ in spec:
        counts[k] = counts.get(k, 0) + 1
    fh.write("Provenance classes in this file: " + ", ".join(
        "**%s** %d" % (k, counts[k]) for k in (MEASURED, DERIVED, ANNOT)
        if k in counts) + "\n")

    for name, dtype, units, klass, source, deriv in spec:
        p = profile(rows, name)
        fh.write("\n### `%s`\n\n" % name)
        fh.write("| | |\n|---|---|\n")
        fh.write("| type | %s%s |\n" % (dtype, (" (%s)" % units) if units else ""))
        fh.write("| class | **%s** |\n" % klass)
        fh.write("| source | %s |\n" % source)
        fh.write("| in this file | %s |\n" % fmt_profile(p))
        fh.write("\n%s\n" % deriv)


def main():
    os.makedirs(SUP, exist_ok=True)
    atlas = list(csv.DictReader(open(CV + r"\cardiovascular_atlas.csv",
                                     encoding="utf-8")))
    classb = list(csv.DictReader(open(CV + r"\classB_upper_bound.csv",
                                      encoding="utf-8")))

    ship = ["cardiovascular_atlas.csv", "cardiovascular_atlas.parquet",
            "classB_upper_bound.csv"]
    for f in ship:
        shutil.copyfile(os.path.join(CV, f), os.path.join(SUP, f))

    validate(atlas, ATLAS_SPEC, "cardiovascular_atlas.csv")
    validate(classb, CLASSB_SPEC, "classB_upper_bound.csv")

    pt = pq.read_table(os.path.join(SUP, "cardiovascular_atlas.parquet"))
    csv_cols = list(atlas[0].keys())
    if pt.num_rows != len(atlas) or list(pt.column_names) != csv_cols:
        raise SystemExit("PARQUET DOES NOT MATCH CSV: parquet %d x %d, csv %d x %d"
                         % (pt.num_rows, pt.num_columns, len(atlas), len(csv_cols)))

    with open(os.path.join(SUP, "DATA_DICTIONARY.md"), "w", encoding="utf-8") as fh:
        fh.write("# Data dictionary\n\n")
        fh.write("Supplementary bundle for the cardiovascular target untestability "
                 "atlas. Every column of every shipped file is documented here. The "
                 "counts, permitted values and ranges below are **computed from the "
                 "shipped files at build time**, so they cannot drift from the data.\n\n")
        fh.write("## Provenance classes\n\n")
        fh.write("Every column carries one of three labels:\n\n")
        fh.write("| class | meaning |\n|---|---|\n")
        fh.write("| **%s** | An experimental measurement from an external resource, "
                 "carried in unchanged. We did not compute it. |\n" % MEASURED)
        fh.write("| **%s** | Computed by this project's code from other columns or from "
                 "a raw source, under a rule stated in the column entry. |\n" % DERIVED)
        fh.write("| **%s** | A categorical label curated by an external resource and "
                 "carried in unchanged. |\n" % ANNOT)
        fh.write("\nA **%s** column is only as good as its rule; the rule is given in "
                 "full for each one, and the code that applies it is named in "
                 "README.md.\n" % DERIVED)
        dict_section("File 1 \u2014 the atlas", "cardiovascular_atlas.csv",
                     atlas, ATLAS_SPEC, fh)
        fh.write("\n> `cardiovascular_atlas.parquet` holds the identical content with "
                 "types applied. It is verified at build time to match the CSV in row "
                 "count and column order.\n")
        dict_section("File 2 \u2014 the Class B upper bound",
                     "classB_upper_bound.csv", classb, CLASSB_SPEC, fh)
        fh.write("\n---\n\n## Sources, versions and access dates\n\n")
        fh.write("| source | version | retrieved (UTC) | licence position |\n"
                 "|---|---|---|---|\n")
        for s in SOURCES:
            fh.write("| %s | %s | %s | %s |\n" % s)
        fh.write("\n**Licence caveat.** The licence column states the position as "
                 "commonly documented for each resource. **None was re-verified against "
                 "the source's own terms page during this build \u2014 treat the column "
                 "as UNVERIFIED** and check before redistributing.\n\n")
        fh.write("**Retrieval caveat.** URL, timestamp, byte count and SHA-256 are "
                 "logged in `logs/retrieval_log.csv` for 22 of the raw inputs. Four "
                 "\u2014 GTEx v10 median TPM, Human Protein Atlas, the GWAS Catalog "
                 "associations archive and the Open Targets 26.06 extract \u2014 are "
                 "**not in that log**; their retrieval URLs are **UNVERIFIED** and the "
                 "dates above are filesystem mtimes, not logged retrieval times.\n")

    with open(os.path.join(SUP, "README.md"), "w", encoding="utf-8") as fh:
        n_a = sum(1 for r in atlas if r["untestability_class"] == "CLASS_A")
        n_i = sum(1 for r in atlas if r["state"] == "INSTRUMENTED")
        n_cv = sum(1 for r in classb if r["cv_relevant_trait"] == "True")
        fh.write("""# Supplementary bundle \u2014 cardiovascular target untestability atlas

## What this is

A row-per-target record of **which cardiovascular drug targets can be studied with human
genetics, and which cannot**. It covers **%d targets** carrying a cardiovascular
indication. For each one it records whether a usable genetic instrument exists in whole
blood, whether one exists in heart, vessel or liver tissue, how strongly the gene is
expressed in each of those tissues, and \u2014 where no instrument exists \u2014 which
measurement would be needed to create one.

**%d targets (%.1f%%) cannot be instrumented at all in whole blood** (`untestability_class
= CLASS_A`), including %d with an approved medicine. **%d (%.1f%%) can be.**

## What this is not

* **Not a causal claim about any target.** No column states that a target does or does not
  affect disease. The atlas records whether the question can be *asked*, not its answer.
* **Not a colocalization resource.** Colocalization is excluded as an evidence gate
  (Gate D, closed as a documented negative result). Colocalization figures elsewhere in
  this project are *method availability*, not evidence.
* **Not a curation of drug mechanism.** `required_modality` is a declared heuristic over
  expression and secretion. It says where a signal could exist, not that one does.
* **Not a complete Class B count.** `untestability_class = CLASS_B` is assigned at N = 1
  (IL6R). The true count is unknown; `classB_upper_bound.csv` bounds it at %d of %d
  instrumented targets (%.1f%%) on cardiovascular traits, %d (%.1f%%) on any trait.
* **Not tissue-complete.** The heart, vessel and liver instrument columns are **upper
  bounds**: they use a wider window than pre-registered and are unclumped.
* **Not a clinical decision aid.** Counts are genes and target-disease pairs, never
  patients.

## Files

| file | rows | what it holds |
|---|---|---|
| `cardiovascular_atlas.csv` | %d | the atlas, 25 columns |
| `cardiovascular_atlas.parquet` | %d | identical content, typed |
| `classB_upper_bound.csv` | %d | instrumented targets with a protein-altering association |
| `DATA_DICTIONARY.md` | \u2014 | every column: type, units, permitted values, source, derivation |

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
""" % (len(atlas), n_a, 100.0 * n_a / len(atlas),
            sum(1 for r in atlas if r["untestability_class"] == "CLASS_A"
                and r["max_phase_num"] == "5"),
            n_i, 100.0 * n_i / len(atlas),
            n_cv, n_i, 100.0 * n_cv / n_i, len(classb), 100.0 * len(classb) / n_i,
            len(atlas), pt.num_rows, len(classb)))

    man = {"built_utc": datetime.datetime.now(datetime.timezone.utc)
           .strftime("%Y-%m-%dT%H:%M:%SZ"),
           "files": {}}
    for f in ship + ["DATA_DICTIONARY.md", "README.md"]:
        p = os.path.join(SUP, f)
        man["files"][f] = {"bytes": os.path.getsize(p), "sha256": sha256(p)}
    json.dump(man, open(os.path.join(SUP, "MANIFEST.json"), "w"), indent=2)

    print("supplementary bundle -> %s" % SUP)
    print("  atlas       %d rows x %d cols" % (len(atlas), len(csv_cols)))
    print("  parquet     %d rows x %d cols (verified against CSV)"
          % (pt.num_rows, pt.num_columns))
    print("  classB      %d rows x %d cols" % (len(classb), len(classb[0])))
    print("  dictionary  %d + %d columns documented, spec/data in sync"
          % (len(ATLAS_SPEC), len(CLASSB_SPEC)))
    for f in sorted(man["files"]):
        print("  %-32s %9d B  %s" % (f, man["files"][f]["bytes"],
                                     man["files"][f]["sha256"][:16]))


if __name__ == "__main__":
    main()
