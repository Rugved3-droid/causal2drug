# Instrument availability — a primary result

Per DECISION 9, instrument availability is a finding about which drug targets the causal-genetics
toolkit can and cannot see, not a data-cleaning statistic. Accumulated deliberately; updated as
each lane produces evidence.

Three states, never merged:

* **ABSENT** — the gene is not tested by the exposure resource at all
* **TESTED-BUT-NULL** — tested, but no cis-eQTL reaches the frozen threshold p < 5e-8
* **INSTRUMENTED** — ≥1 independent cis instrument after clumping

---

## Lane A — eQTLGen phase I whole blood (n ≈ 31,684), Minikel target set

`python src/d1b_census.py` · full detail in `deliverables/d1b/DELIVERABLE_1B.md`

| state | N | % of 2,517 |
|---|---|---|
| UNMAPPED (no unique GRCh37 gene record) | 60 | 2.4 % |
| **ABSENT from eQTLGen** | **652** | **25.9 %** |
| **TESTED-BUT-NULL** | **229** | **9.1 %** |
| **INSTRUMENTED** | **1,576** | **62.6 %** |

Partition verified in code to sum to 2,517.

Per-target states: `deliverables/d1b/mapping_methodB_absent.txt`,
`mapping_methodB_tested_no_instrument.txt`, `clumped_instruments_final.tsv`.

---

## The anchor example: PCSK9

**PCSK9 is ABSENT.** Verified by direct scan of the full 127,341,798-row eQTLGen file: **zero
rows** mention `ENSG00000169174` or the symbol `PCSK9`. PCSK4, PCSK5, PCSK6 and PCSK7 are all
present, so this is not nomenclature.

Why it matters:

* PCSK9 is **liver-expressed and secreted** — a hepatocyte protein acting on the LDL receptor.
* It is the target of one of the most successful cardiovascular drug classes of the last two
  decades (evolocumab, alirocumab, inclisiran).
* It is **completely invisible** to a whole-blood exposure resource.

A pipeline that silently dropped PCSK9 would report a clean result on a target list truncated in a
way **correlated with target biology** — secreted and tissue-restricted proteins removed, blood-
expressed ones retained. That is not random attrition, and it would bias any benchmark of
"does causal-genetic evidence predict clinical success" toward the kinds of target that happen to
be visible in blood.

PCSK9 was retired as a positive control by DECISION 8 and deliberately **not replaced**, because
picking a replacement after seeing eQTLGen's coverage is selection on a known result.

---

## Does absence track tissue-restricted or secreted status?

**ANSWERED 2026-08-30 — measured, confirmed, large effect.** This was previously a hypothesis
resting on the single PCSK9 case; it has now been tested on all 2,517 targets. See
"Class A is now measured" below and `deliverables/d1e/PART_B_CLASS_A.md`.

Confirmed individual cases:

| target | state in eQTLGen | expression class | note |
|---|---|---|---|
| PCSK9 | ABSENT | hepatocyte-secreted | canonical case above |
| HMGCR | INSTRUMENTED | broadly expressed | 277 hits, best p 1.21e-50 |
| IL6R | INSTRUMENTED | broadly expressed, incl. leukocytes | 226 hits, best p 4.09e-34 |

Therapeutic-area gradient (D1b §6), consistent with blood-expression bias but much milder than the
per-target effect: yield highest for infection (75.6 %), haematology (73.5 %), immune (68.2 %);
lowest for endocrine (60.0 %), psychiatry (60.6 %), neurology (61.2 %) — a spread of only 15.6
percentage points, because targets are shared across indications and area-specific target sets
overlap heavily.

---

## Lane B — cardiac/vascular expression (pre-registered, not yet run)

`python src/tissue_power.py` · rule in `PREREGISTRATION.md` §R3b

Instrument availability is the binding constraint on Lane B before a single analysis runs. Minikel
targets with top cis p < 5e-8 (upper bounds — eQTL Catalogue ±1 Mb window, unclumped):

| tissue | n | targets instrumented | % of 2,517 |
|---|---|---|---|
| GTEx v10 liver | 269 | 108 | 4.3 % |
| GTEx v10 artery (coronary) | 279 | 100 | 4.0 % |
| GTEx v10 heart (left ventricle) — **Lane B primary** | 481 | 222 | 8.8 % |
| GTEx v10 heart (atrial appendage) | 473 | 232 | 9.2 % |
| GTEx v10 artery (aorta) | 492 | 268 | 10.6 % |
| GTEx v10 artery (tibial) | 736 | 359 | 14.3 % |
| union of LV + atrial + aorta | — | 422 | 16.8 % |
| **eQTLGen whole blood (measured, clumped, ±100 kb)** | **31,684** | **1,576** | **62.6 %** |

**Cardiac expression data is instrument-limited by roughly four-fold relative to whole blood.**
A Lane B negative result at expression level will usually mean *no instrument*, not *no effect* —
which is exactly the failure mode this document exists to make visible.

---

## Consequence for interpretation, applying to every lane

A target with no instrument produces **no estimate**, never a null. "Not instrumented" and "tested
and null" must remain distinct in every downstream table, model and figure. Collapsing them turns a
coverage property of the exposure resource into apparent biological evidence.

---

## Failure classes (DECISION 12) — classify every target against these

* **CLASS A — target absent or uninstrumented in the exposure resource.** Anchor **PCSK9**
  (0 of 127,341,798 eQTLGen rows). The analysis visibly cannot run.
* **CLASS B — target well-instrumented, but the exposure modality cannot represent the mechanism
  the drug acts on** (splicing, shedding, protein activity, localisation). Anchor **IL6R**
  (F = 134.5, p = 1.16e-59 — a textbook-strong instrument, answer still uninterpretable for the
  drug mechanism).

**Class B is worse than Class A and is invisible on every existing platform, because the analysis
appears to succeed.** Class A announces itself as a missing row; Class B returns a confident,
well-powered, correctly computed estimate of the wrong quantity.

### Class A is now measured, not anecdotal (Part B, 2026-08-30)

| test | effect | p |
|---|---|---|
| ABSENT vs INSTRUMENTED, tissue-restricted (HPA) | **OR 6.90** [5.11, 9.31] | 2.9e-51 |
| ABSENT vs INSTRUMENTED, secreted (HPA) | **OR 1.72** [1.41, 2.11] | 1.8e-07 |
| TESTED-BUT-NULL vs INSTRUMENTED, tissue-restricted | OR 1.46 [1.08, 1.96] | 0.013 |
| TESTED-BUT-NULL vs INSTRUMENTED, secreted | OR 0.82 [0.58, 1.15] | 0.275 — **NULL** |
| GTEx whole-blood median TPM, ABSENT vs INSTRUMENTED | 0.013 vs 5.491; **90.4 % of ABSENT below 1 TPM** | 4.5e-207 |

Attrition is **not random with respect to biology**: it preferentially removes tissue-restricted
and secreted proteins. Full detail: `deliverables/d1e/PART_B_CLASS_A.md`.

### Class B remains unquantified — UNVERIFIED

Part C produced 284 discordant target-indication pairs, 29 with MR p < 0.05. Several are strong
Class B candidates on mechanistic grounds (TYK2 — the protective allele is the coding P1104A
loss-of-function variant, altering kinase activity not transcript level; TNF, CTLA4, IL4R — biologic
targets acting on protein). **No systematic Class B classification was possible**, because no
resource labels "this drug's mechanism operates at a level total mRNA cannot represent". Assigning
Class B currently requires per-target expert judgement.

**Nothing measured so far bounds how large Class B is.** Class B targets sit inside the
INSTRUMENTED group and look entirely normal there — including in every Part B test above.

---

## Infrastructure-level untestability: the `.tbi` filter

Not every untestable target is untestable for a biological or statistical reason. Some are
untestable because the data, though public and complete, cannot be read.

**Measured cost (Part C, `deliverables/d1f/partc_skipped_accessions.tsv`):** of 1,310 qualifying
target-indication pairs, **578 (44.1 %) were lost across 70 GWAS accessions** whose harmonised
summary statistics are **plain gzip with no `.tbi` index**. Those files have complete summary
statistics; they simply cannot be region-queried, and using them would mean downloading 300-500 MB
per accession.

This is the ADDITION 3 category from Deliverable 0 made concrete. The Deliverable 0 sampling
measured **86.50 %** of harmonised GWAS Catalog files as indexed (95 % CI [81.8 %, 91.2 %], n=200),
and **74.49 %** on the LDL-trait census - so the loss is trait-dependent and cannot be applied as a
single global rate.

**Why it belongs in the atlas.** A target-indication pair dropped here is indistinguishable, in any
downstream table, from one that was tested and found null - unless the distinction is carried
explicitly. The attrition is a property of **file packaging**, not of biology, statistics, or the
target. It is the least interesting reason for a target to be untestable and one of the largest.

*(Note: case-control files reporting `odds_ratio` rather than `beta` were NOT lost - beta = ln(OR)
is applied, with the accompanying standard error already on the log-OR scale.)*

---

## CLASS C — PROPOSED, TESTED, **REFUTED**, and DROPPED (2026-08-31)

Class C was to be a third untestability class: *evidence exists but the standard analytic method
fails to detect it, in a way that worsens as data quality improves.*

**The power-dependence half of that definition was pre-registered as falsifiable and has been
falsified.** In the paired head-to-head (`python src/part3_headtohead.py`, 50 targets, same outcome
GWAS and window, only the exposure resource differing):

| | GTEx whole blood n=853 | eQTLGen n=31,684 |
|---|---|---|
| median coloc.abf PP4 | 0.00821 | **0.0139** |
| PP4 higher in the more powerful resource | — | **32 / 50** |
| Wilcoxon signed-rank on log10(PP4) | — | **p = 0.0262** |

**PP4 INCREASES with exposure power.** Per the pre-registered criterion in `PROTOCOL.md`, the
Class C power-dependence claim is **REFUTED and Class C is dropped from the paper.** This matches
**Rosen et al. 2025** (bioRxiv 2025.08.05.668745), who reported the same direction over
n = 420 → 2,256. We predicted the opposite and were wrong.

**Classes A and B are unaffected and remain the core of the atlas.**

### What survives, as a Gate D finding rather than a target class

Two method-failure mechanisms remain documented (see `PROTOCOL.md`, Gate D closed). They describe
when *colocalization as a method* fails; they are **not** properties of a target, so they do not
constitute an untestability class:

1. **ABF peakiness under lead displacement** — published as "class II false negatives",
   **Hukku et al. 2021, AJHG 108(1):25–35**.
2. **LD misspecification at high z** — DENTIST residual-z RMSE (**Chen et al. 2021, Nat Commun
   12:7117**) median **4.46**, **93 % of loci > 2** across 28 measured loci; and a **640×** larger
   ancestry-matched reference (UK Biobank n = 337,000) moved the HMGCR GWAS-side RMSE only
   **7.80 → 7.12**, so panel size is not the lever. The regime described by **Kanai et al. 2022,
   Cell Genomics 2(12):100210 (SLALOM)**.

A practical consequence worth carrying: credible-set colocalization returned an answer at only
**2 of 28 loci (7 %)** with a public LD reference — at the rest it produced no credible-set pair at
all.
