# Part 3 — degradation measurement: **CLASS C POWER-DEPENDENCE IS REFUTED**

```bash
cd /d/causal2drug
python src/part3_headtohead.py    # R14.4, the pre-registered refutation test (no LD needed)
python src/part3_degradation.py   # R14.2/R14.3/R14.5, LD-dependent quantities
```

Rule frozen in `PREREGISTRATION.md` R14 before any of these numbers existed. Selection:
571 Part C pairs → 402 with `n_overlap ≥ 100` → **106 one-per-target (max overlap)** → cap of
**60** at seed 20260830. 22 distinct outcome GWAS.

---

## R14.4 — the head-to-head, and the pre-registered verdict

`coloc.abf` with the **same** outcome GWAS, window and harmonization; only the exposure resource
differs. 60 loci attempted, **50 paired** (10 lost: GTEx blood had < 30 overlapping cis variants).

| | GTEx whole blood n = 853 | eQTLGen n = 31,684 |
|---|---|---|
| median PP4 | 0.00821 | **0.0139** |
| median lead \|z\| | 3.6 (IQR 3.0–6.2) | **24.4 (IQR 16.4–35.6, max 83.8)** |

| test | result |
|---|---|
| PP4 **higher** in the more powerful resource | **32 / 50** |
| PP4 **lower** in the more powerful resource | 18 / 50 |
| exact binomial (lower vs higher) | p = 0.0649 |
| **Wilcoxon signed-rank on log10(PP4)** | **W = 408.0, p = 0.0262** |
| median log10 difference (eQTLGen − GTEx) | **+0.031** |

**PP4 INCREASES with exposure power.** The pre-registered Wilcoxon is significant (p = 0.026) and
the binomial count trends the same way (32 vs 18, p = 0.065). Both point in the same direction.

> **PROTOCOL.md, pre-registered before computing:** *"If PP4 INCREASES with exposure power in our
> paired head-to-head, that REFUTES the Class C power-dependence claim. It will be reported as a
> refutation — not reframed, not narrowed, not explained away. Class C is then dropped."*

**The claim is refuted. Class C is dropped. Classes A and B are unaffected and remain the core.**

This is the direction **Rosen et al. 2025** (bioRxiv 2025.08.05.668745) reported over
n = 420 → 2,256. Our data extend the same direction to n = 853 → 31,684. We predicted the
opposite and were wrong.

### One residual observation, reported but explicitly NOT a rescue

Within the head-to-head, the eQTLGen-minus-GTEx log10 PP4 advantage *shrinks* as the eQTLGen lead
|z| grows (Spearman rho = −0.455, p = 0.00090). That is directionally consistent with the
mechanism we proposed biting at the highest z. **It does not rescue the claim**: the overall
effect is still that more power gives higher PP4, the analysis is correlational and confounded
(a higher eQTLGen z generally means a stronger true eQTL, which also raises GTEx's PP4), and
converting a refuted primary into a surviving subgroup claim is exactly the move the
pre-registration forbids. Recorded as an observation for whoever follows this up.

---

## R14.5 — LD misspecification at scale (DENTIST; Chen et al. 2021, Nat Commun 12:7117)

**28 of 60 loci completed** before the run was stopped for time — the 1000 Genomes LD fetches
stalled repeatedly (single loci blocking 15–20 min in retry chains). **Reported as partial;
the remaining 32 loci are UNVERIFIED.** Checkpoints in `deliverables/d1h/ckpt/`.

| trait side | n loci | median RMSE | IQR | max | % > 2 |
|---|---|---|---|---|---|
| **eQTLGen (median lead \|z\| = 25.6)** | 28 | **4.46** | [3.77, 6.12] | 13.76 | **93 %** |
| GWAS (median lead \|z\| = 2.6) | 27 | 0.71 | [0.51, 0.88] | 1.29 | 0 % |

Expected ≈ 1 under a well-specified reference.

**The HMGCR finding generalises: LD misspecification against 1000G EUR is pervasive on the
high-z side, at 93 % of loci measured.** The GWAS side looks clean here only because these
particular loci have no GWAS signal to misspecify (see the limitation below).

**SuSiE returned a credible-set pair at only 2 of 28 loci (7 %).** At **93 % of loci the LD-aware
method could not return an answer at all** — not a low PP4, no answer. Combined with Part 2, that
is the practical verdict on credible-set colocalization with a public LD reference at this scale.

---

## R14.3 — abf vs susie discrepancy: **NOT TESTABLE**

The test requires loci where *both* methods return a PP4. Only **2 of 28** qualify. With n = 2 no
correlation with eQTL lead z or with r² can be estimated. **UNVERIFIED — not run.**

---

## A design weakness in my own pre-registered rule, stated plainly

**R14.1 selected on `n_overlap` and did not require the outcome GWAS to be genome-wide significant
in the window.** The consequence is visible in the numbers: median GWAS lead |z| = 2.6, and only
**1 of 28** loci has a genome-wide significant GWAS signal.

Colocalization is only meaningful when *both* traits are associated. So the LD-dependent arms
(R14.3, R14.5) largely ran on loci where there is no outcome signal — which is why the GWAS-side
RMSE looks well-behaved and why so few credible-set pairs exist. **This is a flaw in the rule I
wrote, not in the data.**

I have **not** re-filtered on GWAS significance and re-run. Doing so after seeing these results
would be selecting on the quantity of interest. The correct fix is to pre-register a
both-traits-significant criterion and re-run from scratch, which is new work.

**Does this weaken the refutation?** No. R14.4 does not depend on GWAS significance: it is a
*paired* within-locus comparison of two exposures against the identical outcome, so any weakness of
the outcome cancels between the arms. The refutation stands on its own.

---

## R14.5b — intermediate power point: UNAVAILABLE

Amendment 2 required a third power point at n ≈ 2,000–10,000, chosen by sample size alone. eQTL
Catalogue R8 lists exactly two `ge` datasets in range — INTERVAL_RNA `QTD001000` (n = 4,729) and
INTERVAL_RNA_WGS `QTD001003` (n = 2,892) — and **both are absent from the FTP** (`QTS000053`,
`QTS000054` → HTTP 404; the FTP holds QTS000001–QTS000042 only). Cross-checking all 466 R8 `ge`
datasets against the FTP listing gives **0 datasets with n in [2,000, 10,000]**; the largest
FTP-available `ge` dataset is n = 954.

**Non-monotonicity is therefore UNVERIFIABLE with available data.** Two power points give a slope,
not a shape. No curve is claimed. As it happens the slope points the opposite way to our
hypothesis, so the non-monotone question is moot for this project.

---

## Summary

| claim | verdict |
|---|---|
| PP4 falls as exposure power rises (Class C power-dependence) | **REFUTED** (Wilcoxon p = 0.026, direction opposite) |
| LD misspecification is pervasive at high z | **SUPPORTED** — DENTIST RMSE median 4.46, 93 % > 2 (partial, 28/60) |
| Credible-set coloc is usable with a public LD reference at this scale | **NO** — returns an answer at 7 % of loci |
| abf-vs-susie discrepancy grows with z | **NOT TESTABLE** (n = 2) |
| Non-monotone power relationship | **UNVERIFIABLE** — no intermediate resource exists |
