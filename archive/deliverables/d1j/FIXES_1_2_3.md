# Fixes 1–3 — complete, no BLOCKER

```bash
cd /d/causal2drug
python src/test_coloc_assertion.py   # FIX 1
python src/test_susie_coloc.py       # FIX 1 regression (4/4 still pass)
python src/fix2_r15.py               # FIX 2 + FIX 3
```

---

## FIX 1 — hard assertion against the artefact signature

`ColocArtefact` in `src/susie_coloc.py`, active by default (`strict=True`). It **raises**; it
never returns a number silently. A credible-set pair is rejected when either:

* **(a)** both credible sets top at the **same variant** — a component colocalizing with itself; or
* **(b)** the winning pair's tops match **neither trait's lead**.

### It catches all four withdrawn artefacts

| artefact | configuration | PP4 | caught by |
|---|---|---|---|
| HMGCR–LDL-C | 1000G n=525, s²=1 | 0.998 | (a) — both top at index 108 |
| IL6R–CRP | 1000G n=525 | 1.000 | (a) — both top at index 93 |
| HMGCR–LDL-C | UKB LD n=337,000 | 0.9957 | (a) — both top at index 62 (rs1544755, neither lead) |
| HMGCR–LDL-C | 1000G, first draft | 1.000 | (a) — the originally withdrawn value |

**Both conditions are needed.** In the HMGCR/1000G case the shared top *was* the exposure's lead,
so (b) alone would have missed it; in the UKB case the shared top was neither lead, and (a) and (b)
both fire.

### The cost, measured rather than assumed

**Condition (a) also rejects a genuine colocalization.** Synthetic case B1 — ground truth is a
single shared causal variant — is rejected with PP4 = 1. Condition (a) cannot distinguish a real
shared causal variant from a component colocalizing with itself, because both produce the same
pattern.

The assertion is therefore **deliberately conservative**: it will refuse true positives. With Gate
D closed and colocalization excluded from the atlas as an evidence gate, that is the safe
direction. Anyone reusing this code for evidence would need a better discriminator. The synthetic
validation suite calls `strict=False`, with a comment saying why, and still passes **4/4**.

### Retained analyses re-run through it

Both surviving Part 3 credible-set results **pass**: `ADORA1` (susie PP4 = 0.9587) and `F10`
(0.0904). Under FIX 2's corrected locus set, **one further real locus was rejected** — `TNF` —
so the assertion is doing work on live data, not only on the known cases.

---

## FIX 2 — R15: the measurement re-specified and re-run

R14.1 never required the outcome GWAS to be genome-wide significant. R15 adds that criterion
(declared in `PREREGISTRATION.md` before running; descriptive, not confirmatory, since Class C is
dead). Rule: `n_overlap ≥ 100`, **outcome GWAS p < 5e-8 in the window**, one locus per target by
largest overlap, **no cap**.

**106 candidate loci screened. Only 6 (6 %) have a genome-wide significant outcome signal in the
cis-window.** That is itself a finding: for 94 % of instrumented target loci in this set,
colocalization is not interpretable at all, because the outcome trait has no signal there.

| quantity | **R14.1 — flawed rule** | **R15 — corrected** |
|---|---|---|
| N loci | 28 (of 60; 32 superseded) | **6** |
| median GWAS lead \|z\| | 2.6 | **11.3** |
| loci genome-wide significant on GWAS side | 1 / 28 | 6 / 6 (100 %) |
| DENTIST RMSE eQTLGen, median | 4.46 | **3.15** |
| DENTIST RMSE eQTLGen, % > 2 | 93 % | **100 %** |
| DENTIST RMSE GWAS, median | 0.71 | **1.70** |
| DENTIST RMSE GWAS, % > 2 | 0 % | **33 %** |
| credible-set pair returned | 2 / 28 (7 %) | **4 / 6 (67 %)** |
| rejected by the FIX-1 assertion | n/a | **1** |

**The R14.1 row is measured under a flawed selection rule and is labelled as such wherever it
appears.** Its GWAS-side RMSE of 0.71 was an artefact of selecting loci with no GWAS signal —
there was nothing to misspecify. Under R15 the GWAS-side RMSE rises to 1.70 with a third of loci
above 2.

**What survives the correction:** LD misspecification on the high-z eQTLGen side is if anything
*more* pervasive — **100 % of R15 loci exceed RMSE 2** (median 3.15). The finding is robust to the
rule change.

**What does not survive:** the claim that credible-set coloc "returns an answer at only 7 % of
loci". Under a proper locus set it returns an answer at **67 %** (4/6). The 7 % figure was an
artefact of the flawed rule and is withdrawn.

### The six R15 loci

| gene | GWAS | eQTL lead z | GWAS lead z | r² leads | coloc.abf PP4 | coloc.susie PP4 | RMSE GWAS | RMSE eQTL |
|---|---|---|---|---|---|---|---|---|
| ADORA1 | GCST90624412 | 65.68 | 17.56 | 0.0034 | 1.24e-53 | **0.9587** | 1.29 | 4.95 |
| KCNH2 | GCST90624412 | 18.38 | −12.29 | 0.0321 | 0.00206 | **0.9988** | 3.05 | 3.96 |
| TSHR | GCST90319319 | 13.28 | 12.81 | 0.4976 | 0.2448 | — no credible-set pair | 1.55 | 2.96 |
| IL4R | GCST90399686 | 29.34 | 9.214 | 0.3382 | 4.17e-07 | 0.00794 | 1.86 | 3.34 |
| PPARG | GCST90399686 | 31.81 | 5.554 | 0.0002 | 2.29e-05 | 0.00306 | 0.85 | 2.89 |
| TNF | GCST90399686 | 24.48 | 10.36 | 0.1118 | 1.17e-17 | **REJECTED by assertion** | 3.55 | 2.94 |

Two loci show the abf→susie gap dramatically: **ADORA1 1.24e-53 → 0.959** and
**KCNH2 0.00206 → 0.999**. With n = 4 usable pairs this is an observation, not a test —
**R14.3 remains NOT TESTABLE (UNVERIFIED)**.

---

## FIX 3 — the 32 incomplete loci

**Excluded explicitly, with reason: superseded.** They belong to the R14.1 set, whose selection
rule is now known to be inadequate (it admitted loci with no outcome signal). Completing them would
have added ~30 more loci that mostly fail R15's criterion anyway. **Final N for the flawed-rule
comparison is 28**, labelled as such; **final N for the corrected measurement is 6**, complete with
no cap and nothing pending.

Nothing is left UNVERIFIED in a table except where explicitly marked: R14.3 (n = 4, not testable).

---

## No BLOCKER

All three fixes completed. Nothing else was found to be artefactual: of the retained
colocalizations, ADORA1 and F10 pass the assertion, and one new rejection (TNF) was caught on live
data. Proceeding to the cardiovascular subset.
