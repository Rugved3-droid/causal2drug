# Causal2Drug - gate status

| gate | what it is | status | as of |
|---|---|---|---|
| A - instrument retrieval | eQTLGen cis-eQTL, p<5e-8, +-100 kb | OK | 2026-08-30 |
| B - LD clumping | 1000G phase 3 EUR, r2<0.01 | OK | 2026-08-30 |
| C - cis-MR | lead-SNP Wald ratio (PREREGISTRATION R12.2) | OK | 2026-08-30 |
| **D - colocalization** | coloc ABF / coloc.susie | **CLOSED - documented negative result** | 2026-08-31 |
| E - tractability / safety | Open Targets | not started | - |

## Gate D - CLOSED as a documented negative result

Colocalization returned **PP3 ~ 1.000, PP4 ~ 1e-11 at HMGCR-LDL-C**, one of the most robust
colocalizations in human genetics. Near-zero PP4 there indicates an implementation or specification
defect, not a biological result.

**No colocalization output from this project may be used, cited, or carried into any downstream
gate until this is resolved.** MR verdicts do not depend on coloc and are unaffected.

**CAUSE FOUND 2026-08-30** - `deliverables/d1d/PART_A_COLOC_DIAGNOSIS.md`.

Refuted: Z-score conversion (native vs reconstructed beta/SE give PP4 0.0271 vs 0.0403 on the same
data); N imbalance (synthetic case F: GWAS 1.3M vs eQTL 32k with a shared causal variant gives
PP4 = 1); multiple causal variants (exactly 1 independent signal per trait); prior choice (PP4
stays negligible across four orders of magnitude of p12). Implementation validated 6/6 against
synthetic ground truth.

**Actual cause:** coloc-ABF compares single variants and its discrimination scales as (1/2)z^2. The
eQTL lead (rs6453133, z=14.97) and GWAS lead (rs12916, z=49.59) sit 36,237 bp apart at r^2 = 0.464
-- the same LD block, one signal each. The eQTL log-ABF difference between them is 32.23, which
equals the observed s1+s2-s12 exactly and yields PP4 = 1e-11.

**Consequence: this affects every locus in the project and worsens as the exposure resource grows.**

**Remedy ATTEMPTED and FAILED 2026-08-30** - `deliverables/d1g/PART_A_SUSIE_BLOCKER.md`.

SuSiE-RSS + coloc.susie was implemented from scratch and validated 4/4 against synthetic ground
truth, including the multi-signal case it is designed for. On the real HMGCR-LDL-C locus it gives
**PP4 = 0.148**, against a pre-fixed pass criterion of > 0.8. **BLOCKER.**

Cause: the LD reference cannot support fine-mapping a 1.3M-sample GWAS. If one causal variant drove
the region then z_j = z_lead * r(lead,j); the residual RMSE of that relation is **7.80 on the GWAS
side** (should be ~1) versus 1.33 on the eQTL side. A 525-sample 1000G EUR panel cannot reproduce
the LD of a 1,320,016-sample meta-analysis at z ~ 50, so SuSiE-RSS attributes the mismatch to
spurious components and no credible set tops at the true GWAS lead rs12916.

IL6R-CRP coloc.susie PP4 = 1.000 is EXCLUDED from all tables and figures: it is the same
over-fitting artefact as the withdrawn HMGCR PP4 = 1.000, on a locus with no independent sanity
check to expose it.

Also corrected: the earlier recommendation of coloc.susie as the remedy for lead displacement was
imprecise. Its purpose is separating MULTIPLE causal signals; HMGCR has exactly one per trait.

**Gate D therefore has NO trustworthy implementation.** Single-variant coloc is unusable (peakiness)
and credible-set coloc is unusable with this LD panel (misspecification). What would be needed is an
LD reference matched to the outcome GWAS; whether one is obtainable for GWAS Catalog studies is
UNVERIFIED and was not investigated, because that is the "switch resources" move the criterion
forbids without a decision from the user.

---

## PRE-REGISTERED REFUTATION CRITERION — written 2026-08-30, BEFORE Part 3 was computed

Rosen JD, Broadaway KA, Brotman SM, Mohlke KL, Love MI, "Higher eQTL power reveals signals that
boost GWAS colocalization", bioRxiv 2025.08.05.668745, reports the OPPOSITE of the Class C
power-dependence claim over n = 420 -> 2,256: higher eQTL power INCREASES colocalization detection
(eQTL detection 15.6% at n=500 -> 36.8% at n=2,000; colocalizing signals roughly double).

**If PP4 INCREASES with exposure power in our paired head-to-head, that REFUTES the Class C
power-dependence claim. It will be reported as a refutation - not reframed, not narrowed, not
explained away. Class C is then dropped from the paper. Classes A and B are unaffected and remain
the core.**

The measurement is paired, per target, across the declared power points. A null is reported as a
null.

## FRAMING - FIXED 2026-08-30

The contribution is the FOOTPRINT, not the mechanism. All mechanisms are already published and are
cited as such, never as our discovery:

- Hukku et al. 2021 AJHG 108(1):25-35 - lead/causal-variant mismatch causes colocalization false
  negatives ("class II false negatives")
- Kanai et al. 2022 Cell Genomics 2(12):100210 - SLALOM; meta-analysis fine-mapping miscalibration;
  LD reference mismatch "all but guaranteed" in GWAS meta-analyses; 67% of GBMI loci suspicious
- Chen et al. 2021 Nat Commun 12:7117 - DENTIST; observed-vs-LD-predicted residual statistic as the
  standard QC for GWAS/LD-reference heterogeneity
- Wallace 2021 PLOS Genet 17(9):e1009440 - coloc+SuSiE for multiple causal variants

**The residual-z diagnostic used throughout this project is DENTIST's method (Chen et al. 2021),
independently reimplemented here. It is cited as theirs wherever it appears.**


---

## GATE D CLOSED - 2026-08-31

Gate D is closed as a **documented negative result**. It is **not a deliverable**. Colocalization is
**excluded from the atlas as an evidence gate** and appears only as a finding. No further attempts
to fix it will be made.

### Two independent failure mechanisms, both established, both already published

**Mechanism 1 - ABF peakiness under lead displacement.** Single-variant coloc compares one variant
at a time; lBF ~ (1/2)z^2, so the lBF gap between two variants scales as (1/2)(z_j^2 - z_k^2). At
HMGCR the eQTL lead (rs6453133, z=14.97) and GWAS lead (rs12916, z=49.59) sit 36,237 bp apart at
r2=0.464 - one signal per trait - and the gap of 32.23 log units forces PP4 to 1.005e-11. Priors do
not rescue it (PP4 negligible across four orders of magnitude of p12).
Published as "class II false negatives": **Hukku et al. 2021, AJHG 108(1):25-35**.

**Mechanism 2 - LD misspecification at high z.** The DENTIST statistic (**Chen et al. 2021, Nat
Commun 12:7117**, reimplemented here and cited as theirs) gives a GWAS-side residual-z RMSE of
**7.80** with 1000G EUR n=525, where ~1 is expected. Part 2 tested a **640x larger** ancestry-matched
reference (UK Biobank LD, n=337,000, PolyFun/Weissbrod): RMSE moved only to **7.119**. **Panel size
is not the binding constraint.** The outcome GWAS is a heterogeneous meta-analysis (Graham 2021
GLGC, n=1,320,016), whose effective LD is not that of any single reference population - the regime
**Kanai et al. 2022, Cell Genomics 2(12):100210 (SLALOM)** describes as making LD reference mismatch
"all but guaranteed" in GWAS meta-analyses (67% of GBMI loci suspicious).

### Three artefactual PP4 values produced and withdrawn

All three are EXCLUDED from every table and figure. Each has the same signature: a SuSiE component
from each trait converging onto the SAME variant, which colocalizes with itself by construction.

| locus | reference | PP4 | why withdrawn |
|---|---|---|---|
| HMGCR-LDL-C | 1000G n=525, s2 fixed at 1 | 1.000 | outcome CS topped at the EXPOSURE lead rs6453133 (PIP 1.000) while true GWAS lead rs12916 had PIP 5e-87 |
| IL6R-CRP | 1000G n=525 | 1.000 | same over-fitting artefact, on a locus with no independent sanity check to expose it |
| HMGCR-LDL-C | UKB n=337,000 | 0.9957 | winning pair was rs1544755 x rs1544755 - neither trait's lead |

**In no configuration did any credible set top at rs12916, the true GWAS lead.**

### Standing rule

No colocalization output from this project may be used, cited, tabulated or plotted. MR verdicts do
not depend on coloc and are unaffected.

---

## PRE-REGISTERED REFUTATION — OUTCOME, 2026-08-31

**The Class C power-dependence claim was REFUTED.** Paired head-to-head, 50 targets, same outcome
GWAS / window / harmonization, only the exposure resource differing:

- PP4 higher in the more powerful resource: **32 / 50**; lower: 18 / 50
- exact binomial p = 0.0649 ; **Wilcoxon signed-rank on log10(PP4) p = 0.0262**
- median PP4: GTEx blood (n=853) 0.00821 -> eQTLGen (n=31,684) **0.0139**

PP4 INCREASES with exposure power, the direction reported by Rosen et al. 2025
(bioRxiv 2025.08.05.668745). **Class C is dropped from the paper. Classes A and B are unaffected.**

Reproduce: `python src/part3_headtohead.py`

**Non-monotonicity is UNVERIFIABLE**: no eQTL Catalogue R8 `ge` dataset with n in [2,000, 10,000]
exists on the FTP (both INTERVAL sets are QTS000053/054, unreleased; FTP holds QTS000001-042;
largest available n = 954). Two power points give a slope, not a shape. No curve is claimed.

**Known weakness in our own selection rule (R14.1):** it selected on `n_overlap` and did not
require the outcome GWAS to be genome-wide significant, so the LD-dependent arms ran mostly on loci
with no GWAS signal (median GWAS lead |z| = 2.6; 1 of 28 genome-wide significant). This was NOT
re-filtered post hoc. It does not affect the refutation, which is a paired within-locus comparison
where outcome weakness cancels between arms.
