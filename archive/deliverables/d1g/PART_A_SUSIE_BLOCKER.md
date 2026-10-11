# Part A — credible-set colocalization — **BLOCKER**

```bash
cd /d/causal2drug
python src/test_susie_coloc.py    # method validation, synthetic ground truth
python src/run_susie_real.py      # HMGCR-LDL-C and IL6R-CRP, real data
```

## Verdict against the pre-fixed criterion

> **PASS CRITERION (fixed before running): PP4 for HMGCR–LDL-C must recover to > 0.8
> using eQTLGen as the exposure.**

**Measured: coloc.susie best PP4 = 0.1482. The criterion is NOT met. BLOCKER.**

Nothing was tuned to move that number, no resource was swapped to find a better one, and
the criterion was not relaxed. **Gate D remains NOT TRUSTWORTHY.**

| locus | coloc.abf PP4 | coloc.susie best PP4 | criterion |
|---|---|---|---|
| **HMGCR – LDL-C** | 1.005e-11 | **0.1482** | **> 0.8 → FAIL** |
| IL6R – CRP | 7.959e-17 | ~~1.000~~ **EXCLUDED** | none (observational) |

**The IL6R PP4 = 1.000 is EXCLUDED from every table and figure (user directive, 2026-08-30).**
It is not merely untrustworthy: it is *the same over-fitting artefact* as the withdrawn HMGCR
PP4 = 1.000, appearing on a locus that has **no independent sanity check** to expose it. At HMGCR
the artefact was catchable because rs12916 is a known, externally-anchored lead; at IL6R there is
no such anchor, so an artefactual PP4 = 1.000 would have passed unnoticed. The value is struck
through above only so the record shows it was produced and withdrawn; it must not be cited,
tabulated, or plotted.

---

## Method: what was implemented, and why

**SuSiE-RSS (Zou et al. 2022) + coloc.susie (Wallace 2021)**, written from scratch in
`src/susie_coloc.py`. Chosen over conditional/stepwise coloc because SuSiE-RSS never
inverts the LD matrix — it only multiplies by it — and the LD matrices here are singular
(271 variants from 525 samples; measured min eigenvalue −1.5e-14, condition number
**infinite**). Any method requiring an LD inverse would fail outright at this n_snps/n_samples
ratio.

LD: **1000 Genomes GRCh38 high-coverage panel, EUR unrelated, n = 525** — the same panel
used for clumping throughout. Region for HMGCR: chr5:75,236,329–75,464,241 (gene body
±100 kb); 276 harmonized variants, **271 with LD**. Dosages are re-oriented to the
harmonized effect allele before the correlation matrix is formed, so LD signs match z signs.

### Method validation: 4/4 on synthetic ground truth

| case | want | abf PP4 | susie PP4 | result |
|---|---|---|---|---|
| B1 shared single causal variant | PP4 | 1 | 1 | PASS |
| B2 genuinely distinct causal variants | PP3 | 8.2e-42 | 2.3e-36 | PASS |
| **M1 two signals vs one shared** (susie's actual purpose) | PP4 | 1 | 1 | PASS |
| M2 two signals, neither shared | PP3 | 4.6e-64 | 1.4e-47 | PASS |

The implementation does what it should on data where the truth is known.

## Two corrections to my own earlier work, recorded

**1. The Part A remedy recommendation was imprecise.** I previously wrote that coloc.susie
was the remedy for lead-variant displacement. For any single-effect Bayes factor,
`lbf_j = −½log(1+v/s²) + ½(z_j²/s²)·(v/s²)/(1+v/s²)`, so the lbf difference between
variants scales with `(z_j² − z_k²)/s²`. What flattens it is the **residual variance s²** —
not the credible-set structure. coloc.susie's actual purpose is separating **multiple**
causal signals (case M1), which is a different problem from the one at HMGCR, where
Part A already established exactly **one independent signal per trait**.

**2. A first draft of my SuSiE fixed s² = 1 and produced a spurious PP4 = 1 at HMGCR.**
With z up to 49.6, the residual after removing a component is itself large; with s² pinned
at 1 every subsequent single effect became a delta function and IBSS carved the region into
five fake "signals". One of them topped at **rs6453133 — the exposure's lead — with
PIP 1.000, while the true GWAS lead rs12916 had PIP 5e-87.** That produced PP4 = 1.
**I did not report it.** Adding residual-variance estimation removed the artefact and the
honest answer appeared: 0.1482. Recording this because a PP4 = 1 obtained that way would
have looked exactly like a pass.

## Why it fails: the LD reference cannot support fine-mapping a 1.3 M-sample GWAS

The fit is still not sane on the outcome side: **no credible set tops at rs12916**, the
GWAS lead with z = 49.59. An `L = 1` fit does correctly identify rs12916, so the region's
dominant signal is found — the extra components are absorbing something else.

That something is LD misspecification. If one causal variant drives a region, then
`z_j ≈ z_lead · r(lead, j)`. Measuring the departure (the check `susie_rss` performs
internally and warns about):

| trait | n | lead | z | corr(observed, predicted) | **RMSE of residual z** | max abs residual |
|---|---|---|---|---|---|---|
| **GWAS (LDL-C)** | 1,320,016 | rs12916 | 49.59 | 0.9076 | **7.80** | 14.5 |
| eQTL (eQTLGen) | 31,684 | rs6453133 | 14.97 | 0.8854 | **1.33** | 2.9 |

Residual RMSE should be ≈ 1 if the LD reference matched the study's LD. **On the eQTL side
it is 1.33 — acceptable. On the GWAS side it is 7.80 — a factor of ~8 too large.** A
525-sample reference panel simply cannot reproduce the LD of a 1.32 M-sample meta-analysis
finely enough for a summary-statistic fine-mapper operating at z ≈ 50. SuSiE-RSS attributes
that mismatch to additional components, which is why spurious credible sets appear and why
PP4 does not recover.

This is a **data-infrastructure limit, not a coding defect**: the same code passes 4/4 on
synthetic data where the LD is exact.

## What was tried

1. coloc.susie with empirical-Bayes prior variance — **PP4 = 0.148**.
2. Residual-variance estimation added (fixing the artefactual PP4 = 1) — **PP4 = 0.148**.
3. `L = 1` sanity fit — recovers rs12916 correctly, confirming the region's signal is found.
4. LD-consistency diagnostic — identifies the GWAS-side residual RMSE of 7.80 as the cause.

**Not tried, deliberately**, because each would violate the standing instruction: tuning
p12 or the prior grid; using a larger LD reference to obtain a better number; relaxing the
0.8 threshold; substituting a smaller GWAS whose z-scores the panel could support.

## What would actually be needed

An LD reference matched to the outcome GWAS — in-sample LD, or a reference of comparable
ancestry and scale (e.g. UK Biobank-derived LD matrices for UKB-based GWAS). Whether such a
reference is obtainable for the GWAS Catalog studies this project uses is **UNVERIFIED**;
I did not investigate it, because doing so is the "switch resources" move the criterion
forbids without your decision.

Until then: **single-variant coloc is unusable (Part A, established), and credible-set coloc
is unusable with a 525-sample LD panel against large GWAS (established here). Gate D has no
trustworthy implementation.**

**Stopping as instructed. Parts B and C not started.**
