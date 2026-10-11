# Closing analyses 1–3

```bash
cd /d/causal2drug
python src/a1_classB_bound.py    # Analysis 1
python src/a2_robustness.py      # Analysis 2
python src/a3_cascade.py         # Analysis 3
```

---

# ANALYSIS 1 — Class B bounded, not assigned

**Source: GWAS Catalog associations, ontology-annotated full release** (73,544,409 B, retrieved
2026-08-31), 1,191,572 association rows scanned. It carries, per association, the mapped gene and
the variant's **VEP functional consequence** in the `CONTEXT` column — already computed and
versioned by the Catalog, so no VEP calls and no consequence-calling of my own.

*Rejected alternatives, and why:* **ClinVar** is dominated by rare pathogenic variants rather than
the common association signals drug-target genetics rests on. **Open Targets genetics** L2G scores
are model outputs, and putting a trained model inside the denominator of a bound would make the
bound uninterpretable.

**Rule (mechanical):** an instrumented CV target is flagged if the GWAS Catalog reports ≥1
association whose `MAPPED_GENE` contains that gene as an exact token **and** whose `CONTEXT` is
protein-altering (missense, stop_gained, frameshift, start/stop_lost, inframe indel,
protein_altering, splice donor/acceptor, coding_sequence).

## The bound

| set | flagged | total | % |
|---|---|---|---|
| **all instrumented CV targets** | **196** | **386** | **50.8 %** |
| ...with a **launched** drug | 16 | 27 | 59.3 % |
| ...at **Phase III or launched** | 36 | 60 | 60.0 % |
| tightened to **CV-relevant traits only** | **62** | 386 | **16.1 %** |
| ...launched-drug targets, CV traits only | 6 | 27 | 22.2 % |

**The claim is: Class B is at least 1, and plausibly as many as 196 of 386 instrumented CV targets
(50.8 %).** Under the tighter requirement that the protein-altering association is itself for a
cardiovascular trait, the bound is **62 (16.1 %)**.

**This is an upper bound, not an assignment.** A protein-altering association in a gene shows the
mechanism *could* act at protein level; it does not show the drug's mechanism does. Every flagged
target is a candidate, not a case. The atlas keeps `untestability_class = CLASS_B` at N = 1
(IL6R) and carries the bound separately in `classB_upper_bound.csv` (196 rows).

## Worked examples

**IL6R — recovered correctly.** 42 protein-altering associations; the strongest is
**rs2228145, missense_variant**, for blood protein amount and interleukin-6 receptor subunit alpha
measurement. This is exactly the Asp358Ala shedding variant that makes IL6R Class B, so the method
finds the known case unprompted.

**TYK2 — the method works, but TYK2 is not a CV target.** TYK2's Minikel indications are psoriasis
and multiple sclerosis, so it is **not in the 630-target CV set** and could not be flagged here.
Queried directly against the same file, the method finds **142 protein-altering associations across
6 distinct variants, led by rs34536443 (P1104A) with 94 associations** — rheumatoid arthritis,
type 1 diabetes, sarcoidosis, platelet count. The negative in the CV table is a scope boundary, not
a method failure, and is reported as such rather than quietly omitted.

**Flagged launched-drug CV targets** include ADRB1 (rs1801253), ADRB2 (rs1042714, stop_gained),
LPL (rs328, stop_gained), F10, PPARA, PPARG, SERPINC1, PLAT, PLAU, PTGIR, SLC5A2, ITGB2, MTOR,
PDE5A, VEGFA, SLC12A1. At Phase III: CETP (rs1800777), NOS3 (rs1799983, stroke), PLA2G7, GLP1R,
GIPR, DPP4, HGF (stop_gained), IL6R, MPO, GP6, GSK3B (splice_acceptor).

---

# ANALYSIS 2 — every headline result under both CV definitions

| | **CV-PRIMARY** (C14 only) | **CV-EXTENDED** (atlas) |
|---|---|---|
| target-indication pairs | 1,364 | 1,543 |
| **distinct targets** | **594** | **630** |
| ABSENT from eQTLGen | 162 (27.3 %) | 177 (28.1 %) |
| TESTED-BUT-NULL | 53 (8.9 %) | 55 (8.7 %) |
| **CLASS A total** | **215 (36.2 %)** | **232 (36.8 %)** |
| instrumented in blood | 368 (62.0 %) | 386 (61.3 %) |
| yield: whole blood | 62.0 % | 61.3 % |
| **yield: heart LV** | **7.4 %** | **7.1 %** |
| yield: heart atrial | 7.1 % | 6.8 % |
| yield: artery aorta | 9.1 % | 8.7 % |
| yield: liver | 5.2 % | 4.9 % |
| **yield: cardiac union** | **18.0 %** | **17.5 %** |
| **Class A enrichment OR** | **10.60 [4.81, 23.32]** | **10.02 [4.77, 21.03]** |
| Fisher p | 1.83e-14 | 4.05e-15 |
| 2×2 table | [[143, 7], [241, 125]] | [[156, 8], [253, 130]] |
| CLASS A with launched drug | 17 (corrected 2026-10-10; was 20) | 21 |
| CLASS A at Phase III | 18 | 18 |
| CLASS B (assigned) | 1 | 1 |
| modality: cardiac eQTL | 98 (45.6 %) | 101 (43.5 %) |
| modality: tissue eQTL | 46 (21.4 %) | 50 (21.6 %) |
| modality: liver eQTL | 27 (12.6 %) | 35 (15.1 %) |
| modality: vascular eQTL | 18 (8.4 %) | 19 (8.2 %) |
| modality: pQTL | 15 (7.0 %) | 15 (6.5 %) |
| modality: UNKNOWN | 11 (5.1 %) | 12 (5.2 %) |
| **infrastructure attrition** | **218/401 = 54.4 %** | **265/448 = 59.2 %** |

**Stated plainly: every conclusion survives under the stricter definition.** Class A is 36.2 % vs
36.8 %; the cardiac-vs-blood gap is 7.4 % vs 62.0 % (a factor of 8.4) under CV-PRIMARY against
7.1 % vs 61.3 % (8.6) under CV-EXTENDED; the Class A enrichment is *slightly stronger* under
CV-PRIMARY (OR 10.60 vs 10.02) with overlapping intervals; the launched-drug Class A count is 17 of 44 versus 21 of 52 (38.6 % vs 40.4 %). [CORRECTED 2026-10-10: was 20. The 20 took CV-PRIMARY membership but a phase computed over the extended set; with phase taken over CV-PRIMARY indications it is 17 of 44 launched targets. See manuscript_prep/DEFINITIONS.md.]

**The lipid extension is not doing the work.** It adds 36 targets and shifts no headline figure by
more than about one percentage point. **PCSK9 is in the CV set under both definitions** — it
carries cardiovascular (C14) indications as well as hypercholesterolemia — so the anchor example
does not depend on the extension either. The extension is defensible on its merits and, as it
turns out, immaterial to the result.

---

# ANALYSIS 3 — the attrition denominator, audited and **corrected**

## The cascade, step by step

| step | N pairs | N targets | reason for the drop |
|---|---|---|---|
| CV-EXTENDED pairs (frozen set) | **1,543** | 630 | — |
| target instrumented in eQTLGen | 932 | 386 | Class A targets have no instrument |
| ...AND target has a consistent approved-drug direction | 356 | 111 | Part C rule R13.2 (mechanism annotation) |
| ...AND indication exact-matches a EUR harmonised GWAS trait | **103** | 65 | Part C rule R13.4 (trait mapping expedient) |

Cross-checked directly against `partc_pairs.tsv`: 103 CV pairs. Of these, 46 produced a result and
**57 were lost**.

## The cause of the 57 losses — and a labelling error of mine

The skip file records those 57 as `missing_effect_or_se_columns`. **That label is wrong**, because
my Part C runner checks the header before probing for a `.tbi`, and a plain-gzip file yields an
empty header that fails the column check first.

Probing the underlying files directly: all 57 pairs come from just **4 accessions** —
`GCST90029021` (10 pairs), `GCST90043777` (9), `GCST90044006` (23), `GCST90080708` (15) — and
**every one is plain gzip with no `.tbi` (BGZF = False, `.tbi` HTTP 404)**. The substance of the
claim was right; the recorded reason was misleading, and is corrected here.

## 103 is the WRONG denominator — corrected figure

103 is conditioned on two Part C eligibility rules that have nothing to do with file indexing:
the target must carry an unambiguous approved-drug mechanism direction, and the indication must
exact-match a GWAS trait. Quoting an *infrastructure* attrition rate over that denominator mixes
two unrelated filters and understates the base.

**Correct denominator: every CV pair whose indication maps to a EUR-only harmonised GWAS at all** —
i.e. every CV pair for which an outcome GWAS could in principle be used.

| | lost | denominator | rate |
|---|---|---|---|
| ~~Part C-conditioned (as previously reported)~~ | ~~57~~ | ~~103~~ | ~~55.3 %~~ |
| **corrected, CV-EXTENDED** | **265** | **448** | **59.2 %** |
| corrected, CV-PRIMARY | 218 | 401 | 54.4 % |

**The headline infrastructure-attrition figure for CV changes from 55.3 % to 59.2 %.**
`CARDIOVASCULAR_ATLAS.md` has been updated; the 55.3 % figure is retained there only as the
Part C-conditioned sub-rate, labelled as such.
