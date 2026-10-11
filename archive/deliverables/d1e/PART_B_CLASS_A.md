# Part B — Class A quantified: absence is strongly non-random

```bash
cd /d/causal2drug && python src/partb_classA.py
```

States recomputed from the Deliverable 1b artefacts, not carried from its summary:
**INSTRUMENTED 1,576 · ABSENT 652 · TESTED-BUT-NULL 229 · UNMAPPED 60** (sums to 2,517).

**Reference: Human Protein Atlas `proteinatlas.tsv`** (7,460,135 B, sha256 `ad401f8519ecdee5…`,
retrieved 2026-08-30T21:52:32Z). Chosen over a GTEx-derived tau statistic because it supplies
**both** a curated tissue-specificity class **and** a curated secretome annotation from one
consistent source, so the two hypotheses are tested against the same reference rather than two
ad-hoc thresholds of my own choosing. **2,401 of 2,457 mapped targets (97.7 %) join to HPA.**

**Independent cross-check: GTEx v10 median TPM by tissue** (8,846,936 B, sha256
`b951fd32aae5cad4…`), used for the one quantity that is mechanistically direct — whole-blood
expression, which is what eQTLGen can actually see.

Definitions, fixed before testing: *tissue-restricted* = HPA `RNA tissue specificity` ∈
{Tissue enriched, Group enriched, Tissue enhanced}; *not restricted* = {Low tissue specificity};
anything else excluded as NA. *Secreted* = non-empty HPA `Secretome location` or "secreted" in
`Protein class`.

---

## H1 — tissue-restricted expression

| | tissue-restricted | low-specificity |
|---|---|---|
| **ABSENT** | 533 | 53 |
| **INSTRUMENTED** | 926 | 635 |

**OR = 6.90, 95 % CI [5.11, 9.31], Fisher exact p = 2.89e-51**

| | tissue-restricted | low-specificity |
|---|---|---|
| **TESTED-BUT-NULL** | 153 | 72 |
| **INSTRUMENTED** | 926 | 635 |

OR = 1.46, 95 % CI [1.08, 1.96], p = 0.013

## H2 — secreted proteins

| | secreted | not secreted |
|---|---|---|
| **ABSENT** | 215 | 397 |
| **INSTRUMENTED** | 374 | 1,190 |

**OR = 1.72, 95 % CI [1.41, 2.11], Fisher exact p = 1.84e-07**

| | secreted | not secreted |
|---|---|---|
| **TESTED-BUT-NULL** | 46 | 179 |
| **INSTRUMENTED** | 374 | 1,190 |

**OR = 0.82, 95 % CI [0.58, 1.15], p = 0.275 — NULL, and reported as such.**

## H3 — whole-blood expression (GTEx v10, independent of HPA)

| state | n | median TPM | Q1 | Q3 | % below 1 TPM |
|---|---|---|---|---|---|
| **ABSENT** | 634 | **0.013** | 0.000 | 0.076 | **90.4 %** |
| TESTED-BUT-NULL | 227 | 1.157 | 0.112 | 14.625 | 48.9 % |
| INSTRUMENTED | 1,566 | 5.491 | 1.041 | 26.785 | 24.6 % |

Mann-Whitney ABSENT vs INSTRUMENTED: **p = 4.5e-207**. TESTED-BUT-NULL vs INSTRUMENTED:
p = 8.9e-12.

---

## Conclusion

**The hypothesis is confirmed, with a large effect size. It is no longer an anecdote.**

Targets absent from eQTLGen are **~6.9× more likely to be tissue-restricted** and **~1.7× more
likely to be secreted** than instrumented targets. The mechanism is visible directly in the
expression data: **90.4 % of ABSENT targets have GTEx whole-blood median TPM below 1** — they are
simply not transcribed in the tissue the exposure resource measures.

Three things worth separating:

1. **Attrition is not random with respect to biology.** It removes tissue-restricted and secreted
   proteins preferentially. Any benchmark of "does causal-genetic evidence predict clinical
   success" computed on the surviving set is computed on a target list that has been filtered
   along an axis correlated with target class — and, plausibly, with therapeutic area and
   mechanism.

2. **The two failure states are genuinely different, and the data says so.** ABSENT is strongly
   enriched for both tissue-restriction (OR 6.90) and secretion (OR 1.72). TESTED-BUT-NULL is only
   weakly enriched for tissue-restriction (OR 1.46) and **shows no secretion effect at all
   (OR 0.82, p = 0.275)**. That null is reported because it was asked for and because it is
   informative: TESTED-BUT-NULL targets are expressed in blood but lack a strong local regulatory
   variant, which is a statistical-power property, not a tissue-biology property. Collapsing the
   two states would have hidden this.

3. **This quantifies Class A only.** Class B (well-instrumented, wrong exposure modality — the IL6R
   case) is invisible to every test in this file, because Class B targets sit inside the
   INSTRUMENTED group and look entirely normal there. Nothing measured here bounds how large
   Class B is.

Per-target annotations: `deliverables/d1e/classA_annotation.tsv`
(symbol, state, ENSG, HPA specificity class, secretome location, GTEx blood TPM, max tissue TPM).
Test outputs: `classA_tests.json`.
