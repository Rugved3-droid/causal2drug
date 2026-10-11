# DEFINITIONS — "tested-null" and "launched"

> **Update 2026-10-10: reviewed and decided.** Tested-null is now reported as **"no qualifying instrument"** and absent as **"no eQTL data"**, separately everywhere; their sum is only an explicitly labelled upper bound (option C, with sub-types T1–T3 in Supplementary Table S5). The per-set launched definition below is adopted. Window sensitivity at ± 500 kb and ± 1 Mb was added. Corrected numbers on the approved sets are in `KEY_NUMBERS.md`; the counts in this file describe the original 630-target set.

**Status: for review, 2026-10-10. Nothing has been reclassified.**

Reproduce:

```powershell
cd D:\causal2drug
python manuscript_prep\verify\describe_tested_null.py     # Part 1; rows -> tested_null_detail.csv
python manuscript_prep\verify\launched_reconciliation.py  # Part 2
```

---

# Part 1 — What "tested-null" means

## The code that assigns it

Every target gets exactly one of four states. The assignment is in
`src/partb_classA.py:50-58`; `src/cv_atlas.py:139` copies the result into the atlas.

```python
# src/partb_classA.py:41-58   (FLANK = 100_000)
hits = set()
with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c = coords.get(r["MinikelSymbol"])
        if not c or r["SNPChr"] != c["chrom"]:
            continue
        if int(c["start"]) - FLANK <= int(r["SNPPos"]) <= int(c["end"]) + FLANK:
            hits.add(r["MinikelSymbol"])
st = {}
for s in targets:
    if s not in coords:
        st[s] = ("UNMAPPED", None)
    elif coords[s]["ensembl_id"] not in tested:
        st[s] = ("ABSENT", coords[s]["ensembl_id"])
    elif s in hits:
        st[s] = ("INSTRUMENTED", coords[s]["ensembl_id"])
    else:
        st[s] = ("TESTED_NULL", coords[s]["ensembl_id"])
```

The three inputs:

| name in code | file | built by | content |
|---|---|---|---|
| `coords` | `d1b/gene_coords_grch37.tsv` | `src/d1b_gene_coords.py` | one GRCh37 gene record per Minikel symbol (Ensembl release 87 GTF) |
| `tested` | `d1b/eqtlgen_tested_universe.tsv` | `src/d1b_tested_universe.py` | every Ensembl gene ID in the **full** eQTLGen cis-eQTL file: 19,250 genes |
| `eqtlgen_hits_by_ensg.tsv` | same | `src/d1b_hits_by_ensg.py:19` | every eQTLGen row for a Minikel target with **p < 5e-8** |

## The definition in words

A target is **tested-null** when all three hold:

1. its symbol maps to one GRCh37 gene record;
2. eQTLGen tested that gene (its Ensembl ID appears in at least one row of the full file);
3. **no variant with eQTLGen p < 5 × 10⁻⁸ lies within the gene body ± 100 kb.**

I re-applied this rule to all 2,517 targets and it reproduces the shipped states with 0
disagreements.

## What it does not mean

* **Not "no eQTL".** 22 of the 55 CV tested-null targets have a cis-eQTL at eQTLGen's own
  FDR < 0.05.
* **Not "not expressed in blood".** eQTLGen only tests genes it detects in blood. 26 of the 55
  have GTEx whole-blood median TPM ≥ 1, and 17 have TPM ≥ 10.
* **Not a statement about clumping.** State depends on a hit existing in the window. All 1,576
  instrumented targets have at least one independent instrument after clumping, so the two
  never disagree, but clumping is not what assigns the state.
* **Not a null MR result.** No MR is run for these targets.

## What is inside it

Three sub-types, measured from files Deliverable 1b already produced:

| sub-type | meaning | CV set (55) | all targets (229) |
|---|---|---|---|
| **T1** | a genome-wide significant cis-eQTL exists, but every one lies **outside** gene body ± 100 kb | 1 (1.8 %) | 14 (6.1 %) |
| **T2** | a cis-eQTL exists at FDR < 0.05, none reaches p < 5e-8 | 21 (38.2 %) | 92 (40.2 %) |
| **T3** | tested, no cis-eQTL at FDR < 0.05 | 33 (60.0 %) | 123 (53.7 %) |

So tested-null mixes three different situations:

* **T1 is an artefact of the window.** The one CV case is **IL1B**, a launched target: it has 50
  eQTLGen variants at p < 5e-8, the best at p = 4.5 × 10⁻¹⁶, and the nearest is 119,995 bp from
  the gene body, 20 kb outside the window. The ± 100 kb window is pre-registered (R1) and I am
  not proposing to move it.
* **T2 is an artefact of the threshold.** A signal exists and is below the frozen p < 5e-8.
* **T3 is the only group where the resource found nothing** at n = 31,684.

Launched and Phase III CV targets that are tested-null:

| gene | phase | blood TPM | sub-type |
|---|---|---|---|
| IL1B | Launched | 29.5 | T1 — 50 significant variants, nearest 119,995 bp out |
| ITGA2B | Launched | 36.3 | T2 |
| NPC1L1 | Launched | 0.03 | T3 |
| ENPEP, HCN2, KIT | Phase III | 0.01, 0.54, 0.22 | T2 |
| APOA1, BRD4, GRIN1, PDE3A, PDGFRA, SLC6A4 | Phase III | 0.71, 21.1, 0.04, 0.04, 0.10, 0.14 | T3 |

All 55 rows are in `tested_null_detail.csv`.

## Why this matters for the paper

The existing text describes Class A as targets that "cannot be instrumented at all in whole
blood" and explains them by tissue-restricted expression. That fits the **absent** group
(genome-wide OR 6.90 for tissue restriction). It fits tested-null poorly (OR 1.46; secreted OR
0.82, p = 0.27; `d1e/classA_tests.json`), and the enrichment analysis and Figure 4 in fact use
the absent group only.

Options for you to choose between. I have applied none of them.

| option | what changes | effect on the launched count (original set) |
|---|---|---|
| **A** keep Class A = absent + tested-null | wording only: say "no instrument at p < 5e-8 within ± 100 kb" instead of "not measurable in blood" | 21 of 52 |
| **B** report two named categories, "absent" and "tested, no instrument", and never sum them in a headline | headline becomes the absent count | 18 of 52 absent; 3 tested-no-instrument |
| **C** as B, and report T1/T2/T3 in the supplement | adds the sub-types as a column | same as B |

My recommendation is C. The absent group is the finding the biology supports, and separating
the two removes the easiest reviewer objection (ITGA2B and IL1B are abundantly expressed in
blood).

---

# Part 2 — One definition of "launched"

## The definition

> **phase(target, S)** is the highest Minikel `ccatnum` over the target's pairs whose indication
> belongs to set S. A target is **launched in S** when phase(target, S) = 5.

Three points fix the details:

1. **The set that defines membership also defines phase.** A target is launched in a CV set
   only if a CV indication *in that set* carries the Launched label.
2. **`ccat` is Minikel's combined label.** In this dataset the historical/active distinction
   does not arise for launched pairs: all 1,814 Launched pairs have `hcat` = `acat` = `ccat` =
   Launched.
3. **The word is "launched", not "approved".** It is Minikel's label for a target–indication
   pair. Figure 1 and the supplementary README say "approved medicine"; whether that wording is
   justified per target is what the ChEMBL review (`drug_links/`) is for.

`src/cv_atlas.py:110-113` already implements this for the atlas set.

## Reconciling 21/52, 17/44 and 20

| set | launched targets | Class A among them | correct under the definition? |
|---|---|---|---|
| Original 630 (C14 + lipid tokens) | **52** | **21** (40.4 %) | yes |
| CV-PRIMARY (C14 only) | **44** | **17** (38.6 %) | yes |
| CV-PRIMARY as reported in the existing documents | — | 20 | **no** |

Both 21/52 and 17/44 are right; they describe different sets. The **20** is the error. It comes
from `src/a2_robustness.py:84`, which selects CV-PRIMARY targets but reads a phase computed over
the extended set. Three targets are counted as launched under CV-PRIMARY although none of their
CV-PRIMARY indications is launched: ANGPTL3, APOB and NPC1L1, each Preclinical in its C14
indications and launched only in a lipid indication. The fourth lipid-only target, APOC3, is
not in CV-PRIMARY at all. 21 − 4 = 17.

The 20 appears in `ANALYSES_1_2_3.md:88`, `CARDIOVASCULAR_ATLAS.md:170`, `LIMITATIONS.md:274`
and `cv_robustness.json`. Those files stay as they are; the paper will use 17 of 44.

## Effect of the gene-record fixes (provisional)

Four launched targets were unmapped. With the fixes in `UNMAPPED_RESOLUTION.md`, ACE, ITGB3 and
GUCY1A1 are instrumented and GUCY1B1 is tested-null:

| set | launched | absent | tested-null | Class A | instrumented | unmapped |
|---|---|---|---|---|---|---|
| Original 630, as shipped | 52 | 18 | 3 | 21 (40.4 %) | 27 | 4 |
| Original 630, after fixes | 52 | 18 | 4 | 22 (42.3 %) | 30 | 0 |
| CV-PRIMARY, as shipped | 44 | 15 | 2 | 17 (38.6 %) | 23 | 4 |
| CV-PRIMARY, after fixes | 44 | 15 | 3 | 18 (40.9 %) | 26 | 0 |

These are provisional until you approve the mapping fixes. The extended set has not been
computed.

## What I need from you

1. Tested-null: option A, B or C.
2. Confirm the launched definition.
3. Confirm the paper says "launched" throughout, with "approved" reserved for targets the ChEMBL
   review supports.
