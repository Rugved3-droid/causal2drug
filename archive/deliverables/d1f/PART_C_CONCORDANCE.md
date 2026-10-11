# Part C — rule-based distributional validation: **NULL RESULT**

Rule frozen in `PREREGISTRATION.md` §R13 before any concordance was computed.

```bash
cd /d/causal2drug && python src/partc_concordance.py
```

## Headline

**Concordance between cis-MR direction and approved-drug mechanism direction is 50.3 %
(287 / 571), binomial p = 0.933, 95 % CI [0.461, 0.544] — indistinguishable from chance.**

This does **not** validate the pipeline. It also does not, on its own, invalidate it — the test as
constructed is dominated by a scaffold limitation I declared in advance and can now quantify. Both
statements are developed below.

---

## 1. Mechanism extraction, and how reliable it is

Open Targets Platform **26.06** (pinned), `drug_mechanism_of_action` + `drug_molecule`. Both dataset
directories carry `_SUCCESS` (verified). OPEN ITEM 4 caveat carried forward: the release-level
manifest reports `"result": "failure"`.

| step | N | note |
|---|---|---|
| mechanism rows | 6,500 | |
| ...with ≥1 **approved** drug (`maximumClinicalStage == APPROVAL`) | 3,428 | 3,072 excluded |
| ...`single protein` **and** exactly one Ensembl target | 2,131 | 1,021 not single protein; 276 multi/zero target |
| ...with an unambiguous action direction | 1,954 | **177 excluded as ambiguous actionType** |
| distinct targets with ≥1 usable mechanism | 421 | |
| **targets with internally CONFLICTING direction (excluded)** | **44 (10.5 %)** | both UP and DOWN mechanisms |
| targets with a consistent direction | **377** | 293 DOWN, 84 UP |

**Reliability of the extraction itself, as asked:** 177 / 2,131 = **8.3 % of otherwise-usable
mechanism rows have an action type too ambiguous to direct**, and **44 / 421 = 10.5 % of targets
carry mutually contradictory directions** across their own approved drugs. Both were excluded
rather than resolved by majority vote. So roughly one target in ten cannot be assigned a direction
even in principle from this resource — the direction annotation is good but not clean, and a
concordance rate computed on it inherits that noise.

## 2. Cascade to analysable pairs

| step | N |
|---|---|
| targets: instrumented **and** consistent drug direction | **198** |
| target-indication pairs covering them | 5,810 |
| ...whose indication exact-matches a EUR-only harmonised GWAS trait (R13.4) | **1,310** (177 targets, 155 GWAS) |
| ...outcome GWAS actually usable | 732 |
| — lost: **70 accessions skipped**, mostly plain-gzip with **no `.tbi`** (not range-queryable) | **−578** |
| ...with variants surviving harmonization | **571** |
| — no GWAS chunks in window 26; no overlap after harmonization 135 | −161 |

The 578-pair loss is the Deliverable 0 / ADDITION 3 category made concrete: those files have full
summary statistics but are **not queryable**, and using them would mean downloading 300–500 MB
each. Case-control files reporting `odds_ratio` instead of `beta` were **not** lost — β = ln(OR) is
applied, with the accompanying `standard_error` already on the log-OR scale.

## 3. Concordance

**Pair-level** — pairs are **not independent**; many share a target, so this test overstates
precision:

> **287 / 571 = 50.3 %**, binomial vs 50 %: **p = 0.933**, 95 % CI **[0.461, 0.544]**

**Target-level** (one vote per target, majority across its indications; 30 ties excluded):

> **64 / 121 = 52.9 %**, binomial p = 0.586, 95 % CI [0.436, 0.620]

**Stratified by instrument count:**

| instruments | concordant | % | p |
|---|---|---|---|
| 1 | 74/139 | 53.2 % | 0.498 |
| 2 | 66/139 | 47.5 % | 0.611 |
| ≥3 | 147/293 | 50.2 % | 1.000 |

**Stratified by F (lead SNP):**

| F quartile | range | concordant | % |
|---|---|---|---|
| 1 | 31–82 | 68/142 | 47.9 % |
| 2 | 82–354 | 69/143 | 48.3 % |
| 3 | 361–857 | 80/143 | 55.9 % |
| 4 | 861–37,530 | 70/143 | 49.0 % |

**Stratified by MR significance:**

| MR p | concordant | % | p |
|---|---|---|---|
| < 5e-8 | 4/7 | 57.1 % | 1.000 |
| < 0.05 | 35/64 | 54.7 % | 0.532 |
| ≥ 0.05 | 252/507 | 49.7 % | 0.929 |

**No stratum shows concordance above chance.** Stronger instruments do not help; significant MR
estimates do not help.

## 4. Why — the dominant limitation, and it is not subtle

**89 % of the MR estimates (507/571) are non-significant.** The direction of a null estimate is a
coin flip, so a test dominated by nulls returns 50 % almost by construction. There is *some* real
signal — the fraction with p < 0.05 is **0.112**, about twice the 0.05 expected under a global null
— but it is thin.

**The deeper problem is which outcome trait each target got paired with.** R13.4 declared the
exact-string indication→trait rule as an unvalidated Part-C-only expedient. Its cost is now
measurable, and HMGCR is the clean illustration:

| HMGCR pair in this analysis | MR p |
|---|---|
| HMGCR → Diabetes Mellitus | 0.137 |
| HMGCR → Aneurysm | 0.042 |
| HMGCR → Pneumonia | 0.511 |
| HMGCR → Alzheimer Disease | 0.137 |

**HMGCR is never paired with LDL cholesterol** — the phenotype through which statins actually act,
and the pair that returned β = +0.389, p = 4.2e-229 in the accepted Part C control run. Minikel's
HMGCR indications are clinical endpoints, and exact string matching to GWAS trait names supplies
"Pneumonia" and "Alzheimer Disease" instead. The most-used outcome traits overall are Psoriasis
(62 pairs), Asthma (42), Multiple Sclerosis (41), Alzheimer Disease (41) — diseases, not the
biomarker phenotypes through which most of these drugs are understood to work.

**So the test is not measuring what it was intended to measure.** It asks "does target expression
affect this disease in the direction the drug implies", when the drug's mechanism often runs
through an intermediate phenotype that is not the trait being tested.

## 5. Discordant pairs, and Class B

284 discordant pairs (`partc_discordant.tsv`); 29 of them have MR p < 0.05. The strongest:

| target | indication | got | expected | β | p |
|---|---|---|---|---|---|
| TYK2 | Psoriasis | NEGATIVE | POSITIVE | −0.609 | 7.5e-49 |
| ADORA1 | Atrial Fibrillation | POSITIVE | NEGATIVE | +0.048 | 1.2e-12 |
| IL4R | Asthma | NEGATIVE | POSITIVE | −0.136 | 2.8e-11 |
| DNMT3A | Psoriasis | NEGATIVE | POSITIVE | −0.447 | 1.4e-06 |
| PTGIR | Asthma | POSITIVE | NEGATIVE | +0.072 | 2.5e-05 |
| TNF | Asthma | NEGATIVE | POSITIVE | −0.065 | 9.3e-05 |
| CTLA4 | Psoriasis | NEGATIVE | POSITIVE | −0.239 | 2.6e-04 |

**Class B classification is not systematically possible, and that is itself the finding.** Several
of these are strong candidates on mechanistic grounds — TYK2's protective allele is the coding
P1104A loss-of-function variant, which reduces kinase *activity* rather than transcript level;
TNF, CTLA4 and IL4R are targeted by biologics acting on protein, not on mRNA abundance. But there
is no resource that labels "the drug's mechanism operates at a level total mRNA cannot represent",
so assigning Class B requires per-target expert judgement.

**Marked UNVERIFIED:** no systematic Class B/unexplained split of the 284 discordant pairs was
produced, because the annotation required to do it mechanically does not exist. Class B remains
what DECISION 12 says it is — invisible, because the analysis appears to succeed.

## 6. What this does and does not establish

**Establishes:**
* End to end, on 571 pairs, the pipeline's MR direction carries **no detectable information** about
  approved-drug mechanism direction (50.3 %, CI [0.461, 0.544]).
* The mechanism-direction annotation is itself ~10 % unreliable, and that is now quantified rather
  than assumed.
* 578 of 1,310 qualifying pairs are unreachable because their GWAS files are not range-queryable —
  the ADDITION 3 category, now with a concrete cost.

**Does not establish:**
* That the MR machinery is wrong. It produces strong, plausible estimates where signal exists
  (TSHR→Hyperthyroidism, KCNH2→Atrial Fibrillation, JAK2→Crohn Disease all concordant at p < 1e-4),
  and it reproduced HMGCR→LDL-C at p = 4.2e-229 when given the right outcome trait.
* That the null is biological. It is at least as consistent with the outcome traits being the wrong
  phenotypes, which is a property of the R13.4 scaffold and not of the pipeline.

**Consequence:** Part C cannot serve as the pipeline's validation until Deliverable 3 supplies a
real indication→trait mapping. Re-running it afterwards is the correct test, and the 50.3 % here is
the pre-mapping baseline against which that re-run should be compared. **The concordance rate is
reportable as it stands, and it is 50.3 %.**
