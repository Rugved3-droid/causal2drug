# Cardiovascular Target Untestability Atlas

```bash
cd /d/causal2drug
python src/cv_atlas.py          # set definition, classes, tissue yield, enrichment
python src/cv_atlas_report.py   # named examples, modality, method availability
```

Resource: `cardiovascular_atlas.csv` and `cardiovascular_atlas.parquet` — **one row per CV
target, 630 rows, 25 columns**. Summaries: `cv_summary.json`, `cv_report.json`.

---

## The CV set rule, declared before applying

**Classification source:** Minikel's own `areas.tsv` × `indic_topl_match.tsv` (MeSH top-level
tree code → therapeutic area). All 1,014 indications map, zero misses.

* **CV-PRIMARY** — therapeutic area == `cardiovascular` (MeSH C14): **1,364 pairs, 594 targets**
* **CV-EXTENDED** — CV-PRIMARY ∪ indications whose MeSH term matches a fixed token list
  {`hyperlipid`, `hypercholesterol`, `dyslipid`, `hyperlipoprotein`, `hypertriglycerid`,
  `lipoprotein`, `atheroscler`, `arterioscler`}: **1,543 pairs, 630 targets** (+36 targets)

**The lipid/metabolic decision, made explicitly.** Minikel classifies the atherogenic lipid
disorders as **metabolic (C18), not cardiovascular** — seven indications: Dyslipidemias,
Hypercholesterolemia, Hyperlipidemias, Hyperlipoproteinemia Type I, **Hyperlipoproteinemia
Type II**, Hypertriglyceridemia, Hypoalphalipoproteinemias.

**CV-EXTENDED is the atlas definition.** Hypercholesterolemia and Hyperlipoproteinemia Type II are
the PCSK9 and statin indications; a cardiovascular atlas that excluded them would exclude its own
anchor example and the entire lipid-lowering class. The rule is a fixed token match on the MeSH
term — mechanical and reproducible, not case-by-case judgement. CV-PRIMARY is reported alongside so
the effect of the choice is visible.

---

## 1. Class breakdown — 630 CV targets

| state | N | % |
|---|---|---|
| **ABSENT from eQTLGen** | **177** | **28.1 %** |
| **TESTED-BUT-NULL** | **55** | **8.7 %** |
| **CLASS A total (absent + tested-null)** | **232** | **36.8 %** |
| instrumented in whole blood | 386 | 61.3 % |
| unmapped | 12 | 1.9 % |
| **CLASS B (modality mismatch)** | **1** (IL6R) | see §4 |

**More than one in three cardiovascular drug targets cannot be instrumented at all in the
standard whole-blood exposure resource.**

## 2. Tissue yield — CV-specific, not genome-wide

| exposure resource | CV targets instrumented | % of 630 | genome-wide |
|---|---|---|---|
| **whole blood (eQTLGen, n = 31,684)** | 386 | **61.3 %** | 62.6 % |
| GTEx heart left ventricle (n = 481, Lane B primary) | 45 | **7.1 %** | 8.8 % |
| GTEx heart atrial appendage (n = 473) | 43 | 6.8 % | 9.2 % |
| GTEx artery aorta (n = 492) | 55 | 8.7 % | 10.6 % |
| GTEx liver (n = 269) | 31 | 4.9 % | 4.3 % |
| **cardiac union (4 tissues)** | **110** | **17.5 %** | 16.8 % |

CV targets behave almost exactly like the genome-wide set: blood 61.3 % vs 62.6 %, cardiac union
17.5 % vs 16.8 %. **The cardiac tissue set is ~3.5× worse than blood even for cardiovascular
targets** — the tissue that is mechanistically right is the one that cannot instrument them.

## 3. Class A enrichment within CV targets — **stronger than genome-wide**

| | tissue-restricted | low-specificity |
|---|---|---|
| **ABSENT** | 156 | 8 |
| **INSTRUMENTED** | 253 | 130 |

**OR = 10.02, 95 % CI [4.77, 21.03], Fisher exact p = 4.05 × 10⁻¹⁵**

Genome-wide primary result: OR 6.90, 95 % CI [5.11, 9.31], p = 2.9 × 10⁻⁵¹.

**CV is not underpowered here and the effect is if anything larger** (OR 10.0 vs 6.9, confidence
intervals overlapping). The genome-wide result remains the primary on precision; the CV-restricted
result is a fully-powered replication in the target population of interest.

## 4. Named examples — the atlas core

**CLASS A anchor — PCSK9.** ABSENT from eQTLGen: **0 of 127,341,798 rows** (direct scan; PCSK4/5/6/7
present, so not nomenclature). GTEx whole-blood median **0.29 TPM** vs **liver 28.3 TPM**. Launched
drug class (evolocumab, alirocumab, inclisiran). **Required modality: liver eQTL.**

**CLASS B anchor — IL6R.** INSTRUMENTED and strongly so: 2 independent instruments, F = 134.5,
lead p = 1.16 × 10⁻⁵⁹, blood 60.1 TPM. cis-MR runs and returns a confident POSITIVE estimate. It is
still uninterpretable for the drug mechanism: tocilizumab's effect derives from rs2228145
(Asp358Ala) **shedding**, which acts at the splicing/protein level and total mRNA cannot represent.
**Required modality: sQTL or pQTL.** The analysis appears to succeed — that is what makes Class B
worse than Class A.

**CLASS A targets with LAUNCHED drugs (21).** These are approved cardiovascular medicines whose
targets cannot be instrumented in blood:

| gene | state | blood TPM | heart LV TPM | liver TPM | required modality |
|---|---|---|---|---|---|
| **PCSK9** | ABSENT | 0.29 | 0.0 | 28.3 | liver eQTL |
| **APOC3** | ABSENT | 0.79 | 0.1 | **8,122** | liver eQTL |
| **ANGPTL3** | ABSENT | 0.02 | 0.0 | 127.7 | liver eQTL |
| **APOB** | ABSENT | 0.02 | 1.8 | 350.4 | cardiac eQTL |
| **F2** (thrombin) | ABSENT | 0.09 | 0.1 | 598.4 | liver eQTL |
| **PLG** | ABSENT | 0.02 | 0.0 | 288.6 | liver eQTL |
| **REN** | ABSENT | 0.00 | 0.0 | 0.4 | pQTL |
| **AGTR1** | ABSENT | 0.03 | 1.4 | 28.2 | cardiac eQTL |
| **EDNRA / EDNRB** | ABSENT | 0.02 / 0.53 | 6.2 / 10.4 | 1.1 / 7.8 | cardiac eQTL |
| **ADRA1A / 1B / 1D / 2B / 2C** | ABSENT | ≤0.05 | 0.1–3.4 | 0.0–13.3 | cardiac / vascular eQTL |
| **NPC1L1** | TESTED-NULL | 0.03 | 0.1 | 31.1 | liver eQTL |
| **ITGA2B** | TESTED-NULL | 36.3 | 0.3 | 0.4 | cardiac eQTL |
| **IL1B** | TESTED-NULL | 29.5 | 0.1 | 0.5 | vascular eQTL |
| **RXFP1** | ABSENT | 0.00 | 0.0 | 0.6 | tissue eQTL (brain cortex) |

**CLASS A at Phase III (18)** includes **LPA** (blood 0.00 TPM, liver 14.6 — the lipoprotein(a)
target of pelacarsen), **APOA1**, **HCN4** (ivabradine target), **PDE3A**, **AVPR2**, **IL17A**,
**ABCC8**.

The full 630-row list with every field is the CSV/Parquet resource; `untestability_class` and
`max_clinical_phase` are the columns to filter on.

## 5. Required modality for the 232 untestable CV targets

Rule, declared and mechanical: cardiac eQTL if GTEx heart LV or atrial TPM ≥ 1; else liver eQTL if
liver TPM ≥ 1; else vascular eQTL if aorta/coronary TPM ≥ 1; else pQTL if HPA-secreted; else tissue
eQTL naming the max-TPM tissue; else **UNKNOWN**. IL6R is assigned from its established mechanism.

| modality that would be required | N | % |
|---|---|---|
| **cardiac eQTL** | 101 | 43.5 % |
| tissue eQTL (non-cardiac, non-hepatic) | 50 | 21.6 % |
| **liver eQTL** | 35 | 15.1 % |
| vascular eQTL | 19 | 8.2 % |
| **pQTL** | 15 | 6.5 % |
| **UNKNOWN** | 12 | 5.2 % |

**Nearly half the untestable CV targets would need cardiac eQTL data** — the resource measured in
§2 as instrumenting only 7.1 % of CV targets. The requirement and the shortfall point at the same
gap.

## 6. Method availability, CV-restricted

**Infrastructure attrition — CORRECTED 2026-08-31 (Analysis 3).** Of **448** CV pairs whose
indication maps to a EUR-only harmonised GWAS — every CV pair for which an outcome GWAS could in
principle be used — **265 (59.2 %) are lost** because the GWAS summary statistics are plain gzip
with no `.tbi` index and cannot be region-queried. Under CV-PRIMARY the rate is 218/401 = 54.4 %.
Both exceed the 44.1 % overall rate.

*Superseded:* an earlier draft quoted **57/103 = 55.3 %**. That denominator was conditioned on two
Part C eligibility rules (approved-drug mechanism direction; exact trait match) which have nothing
to do with file indexing, and so mixed two unrelated filters. 57/103 is retained only as the
Part C-conditioned sub-rate. The 57 losses were verified to come from 4 accessions, all plain gzip
with no `.tbi`.

**Colocalization.** Of 65 CV target loci screened, only **3 (4.6 %)** have a genome-wide
significant outcome signal in the cis-window, and only **2 (3.1 %)** yield a usable credible-set
colocalization.

> **Colocalization is unavailable for 96.9 % of cardiovascular target loci.**

This is method availability, not evidence: Gate D is closed as a documented negative result
(`PROTOCOL.md`) and colocalization is excluded from the atlas as an evidence gate.

---

## Robustness of the CV set definition (Analysis 2)

Every headline result was recomputed under **CV-PRIMARY** (594 targets, Minikel C14 only) and
**CV-EXTENDED** (630, the atlas definition). **Every conclusion survives.** Class A 36.2 % vs
36.8 %; heart LV 7.4 % vs 7.1 % against blood 62.0 % vs 61.3 %; Class A enrichment OR 10.60
[4.81, 23.32] vs 10.02 [4.77, 21.03]; launched-drug Class A 17 of 44 vs 21 of 52. [CORRECTED 2026-10-10: was 20. The 20 took CV-PRIMARY membership but a phase computed over the extended set; with phase taken over CV-PRIMARY indications it is 17 of 44 launched targets. See manuscript_prep/DEFINITIONS.md.] **PCSK9 is in the set under
both definitions**, so the anchor does not depend on the lipid extension. Full table in
`ANALYSES_1_2_3.md`.

## Class B upper bound (Analysis 1)

Class B remains **assigned** at N = 1 (IL6R). It is now also **bounded**.

**Primary bound — 62 of 386 instrumented CV targets (16.1 %)** carry a protein-altering GWAS
Catalog association *for a cardiovascular trait*. This is the bound to cite: it requires the coding
signal to be relevant to the disease area the drug treats.

**Loose bound — 196 (50.8 %)** if any trait is allowed (59.3 % of launched-drug targets, 60.0 % of
Phase III-or-launched). Reported for completeness and **discounted**: a protein-altering
association for an unrelated trait says little about the cardiovascular drug mechanism.

The claim is *at least 1, plausibly ~62* — an upper bound, not an assignment.
See `classB_upper_bound.csv`.

## Caveats

* **Class B is under-counted at N = 1.** Only IL6R is assigned, because no resource labels "this
  drug's mechanism operates at a level total mRNA cannot represent". The true Class B count among
  the 386 instrumented CV targets is **UNVERIFIED** and is probably much larger — every one of them
  looks normal in every table here.
* **Cardiac instrumentation figures are upper bounds**: they use the eQTL Catalogue permuted files'
  top cis p-value over a ±1 Mb window, wider than the pre-registered ±100 kb, and unclumped.
* The modality assignment is a declared heuristic over expression and secretion, not a curation of
  drug mechanism. It says where a signal could exist, not that one does.
* 12 CV targets (1.9 %) are UNMAPPED — no unique GRCh37 gene record — and are excluded from the
  class statistics rather than assigned.
