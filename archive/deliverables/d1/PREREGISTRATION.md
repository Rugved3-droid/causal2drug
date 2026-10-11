# Deliverable 1 — pre-registration

**Written 2026-08-26, before any MR or colocalization result was computed.**
Nothing below may be changed in response to a result. If a control fails under these rules, the
failure is reported as a BLOCKER and the deliverable stops.

---

## R1. Cis-window definition

Gene body ± 100,000 bp on GRCh38. Gene bounds are taken from the Ensembl REST
`lookup/symbol/homo_sapiens/<SYMBOL>` endpoint, recording the Ensembl gene ID, assembly name and
retrieval timestamp. The same window is used for the GWAS extract, the eQTL extract and the LD
extract.

## R2. Outcome GWAS selection rule

For a given trait EFO ID, select the study that satisfies **all** of:

1. `MAPPED_TRAIT_URI` in the GWAS Catalog studies release `v1.0.3.1` contains the EFO ID;
2. the accession appears in `harmonised_list.txt`;
3. a sibling `.h.tsv.gz.tbi` exists (HTTP HEAD → 200);

and among those, has the **largest** sample size, parsed as the maximum integer appearing in the
`INITIAL SAMPLE SIZE` field. Ties broken by lowest accession number.

Traits: LDL cholesterol = `EFO_0004611`; C-reactive protein = `EFO_0004458`.

## R3. eQTL tissue selection rule — FIXED IN ADVANCE, ONE TISSUE FOR ALL THREE PAIRS

**Rule: GTEx_v10 liver, eQTL Catalogue R8 dataset `QTD000266`, `quant_method = ge`, for all three
control pairs. No per-pair tissue selection. No maximum taken over tissues.**

Justification, fixed before running: all three traits are hepatically determined. LDL-C is set by
hepatic cholesterol synthesis (HMGCR, the statin target) and hepatic LDL-receptor recycling (PCSK9,
secreted by hepatocytes). CRP is an acute-phase protein synthesised by hepatocytes downstream of
IL-6 binding IL6R on the hepatocyte surface. Liver is therefore the a priori site of action for
every pair.

Liver is additionally the *only* liver `ge` dataset in eQTL Catalogue R8, so the rule has a unique
solution and no discretion:

```bash
awk -F'\t' 'NR>1 && $9=="ge" && tolower($6) ~ /liver/ {print $1,$2,$3,$6,$8}' \
  data/raw/dataset_metadata_r8.tsv     # -> QTS000015 QTD000266 GTEx_v10 liver 269
```

Fixing one tissue in advance is also what keeps this control honest with respect to the
Deliverable 6 leakage concern about taking a maximum over ~49 correlated GTEx tissues. GTEx_v10
contributes 49 `ge` datasets to R8; none of the other 48 may be substituted if liver fails.

## R4. eGene requirement

The gene must be a significant eGene in the selected tissue, judged by the dataset's **own**
permutation output `QTD000266.permuted.tsv.gz`, using `p_beta` with Benjamini–Hochberg FDR
across all genes in that file, threshold **FDR < 0.05**.

## R5. Instrument selection and clumping

* Candidate instruments: variants inside the cis-window with eQTL `pvalue` < **5 × 10⁻⁸**.
* LD from 1000 Genomes GRCh38 phased panel, **EUR superpopulation, unrelated only**
  (`FatherID == 0 and MotherID == 0`), retrieved by range query.
* Clumping: greedy, lowest eQTL p-value first, remove any remaining variant with
  **r² ≥ 0.01** against an already-selected instrument.
* Instrument strength: F = (β_eQTL / SE_eQTL)² per instrument. Reported as min, mean, and the
  N-weighted overall F. **Strength criterion: F ≥ 10.**

## R6. Harmonization rules

Both sources and the LD panel are GRCh38; this is asserted by checking the eQTL variant ID
(`chr_pos_ref_alt`) against the reported position, and by confirming the GWAS Catalog file is the
`.h.` harmonised (GRCh38) product. Build alignment is verified, not assumed.

Merge key is (chromosome, position). For each merged variant:

1. If GWAS {effect, other} == eQTL {alt, ref} → aligned, sign of β_GWAS kept.
2. If GWAS {effect, other} == eQTL {ref, alt} → flipped, sign of β_GWAS negated.
3. Otherwise, attempt reverse-complement of the GWAS pair and retry 1–2; if it then matches, the
   variant is recorded as strand-flipped and counted.
4. If none match → dropped, counted as `allele_mismatch`.

**Strand-ambiguous (palindromic) SNPs** — allele pairs A/T and C/G — are handled explicitly rather
than silently: they are resolved by allele frequency when GWAS `effect_allele_frequency` and eQTL
alt-allele frequency (`ac/an`) are both available and both < 0.42 or both > 0.58, requiring the two
frequencies to agree in direction; otherwise the variant is **dropped** and counted as
`ambiguous_unresolved`. Counts for every category are reported.

Indels and multi-allelic records are dropped and counted.

## R7. MR estimator

* Exactly 1 instrument → **Wald ratio**: β = β_out/β_exp, SE = |SE_out / β_exp| (first-order delta
  method), z = β/SE, two-sided p.
* ≥ 2 instruments → **inverse-variance-weighted fixed effects** over the clumped, near-independent
  instruments.

## R8. Colocalization

`coloc` approximate Bayes factor method (Giambartolomei et al. 2014), implemented directly.
Priors **p1 = 1e-4, p2 = 1e-4, p12 = 1e-5**. Both traits treated as quantitative with prior effect
SD **W = 0.15** (the coloc default for quantitative traits). Tested over all variants surviving
harmonization in the cis-window. PP0–PP4 reported, along with the number of variants tested.

## R9. Expected directions — committed before running

MR estimates are expressed as change in trait per **+1 unit (SD) of gene expression in liver**.

| pair | expected sign | basis |
|---|---|---|
| HMGCR → LDL-C | **POSITIVE** | HMGCR is rate-limiting for cholesterol synthesis; statins inhibit HMGCR and lower LDL-C, so more HMGCR ⇒ more LDL-C |
| PCSK9 → LDL-C | **POSITIVE** | PCSK9 degrades the LDL receptor; PCSK9 inhibitors lower LDL-C, so more PCSK9 ⇒ more LDL-C |
| IL6R → CRP | **NEGATIVE** (declared ambiguous in advance) | see below |

**IL6R direction is mechanistically ambiguous at the mRNA level, and I am committing to a
prediction now rather than after seeing the answer.** Pharmacological IL6R blockade (tocilizumab)
lowers CRP, which naively implies more IL6R signalling ⇒ more CRP ⇒ positive. But the dominant cis
signal at *IL6R* is rs2228145 (Asp358Ala), which *increases* IL6R mRNA and increases receptor
shedding; the resulting soluble IL6R sequesters IL-6 and **reduces** classical hepatic signalling,
and the higher-expression allele is associated with **lower** CRP. Because the lead cis-eQTL is
expected to tag rs2228145, the committed prediction is **NEGATIVE**.

If IL6R→CRP returns a strong POSITIVE association, that is recorded as a **direction failure under
this pre-registration** and reported as such — I will not retroactively reinterpret the mechanism
to convert it into a pass.

## R10. Pass criteria — committed before running

A control **PASSES** only if all of:

1. gene is an eGene at FDR < 0.05 in liver (R4);
2. ≥ 1 instrument at p < 5e-8 after clumping (R5);
3. F ≥ 10;
4. MR p < 0.05 (non-null);
5. MR sign matches R9.

Any failure ⇒ **BLOCKER**, reported with the failing number, and the deliverable stops. No
threshold, tissue, window, or instrument criterion is changed to convert a failure into a pass.

## R11. Failure-mode diagnostic, defined in advance

HMGCR and PCSK9 are run through the identical code path in the identical tissue as IL6R. If those
two pass and IL6R alone fails, the pipeline machinery is demonstrably intact and the failure is
attributable to the fixed tissue choice or to the biology, not to broken code. This distinction is
recorded but does **not** convert a fail into a pass, and does not license a tissue swap.

---

## AMENDMENT R2a — 2026-08-26, made BEFORE any MR or colocalization result was computed

**What happened.** Applying R2 to `EFO_0004611` produced a tie at "max parsed N = 1,320,016"
between `GCST90239655` and `GCST90239658`, and the tie-break (lowest accession) selected
`GCST90239655`. Inspecting the two:

* `GCST90239658` — `1,320,016 European ancestry individuals` (single-ancestry European)
* `GCST90239655` — `40,963 South Asian, 48,057 Hispanic, 99,432 Admixed African or African,
  146,492 East Asian, 1,320,016 European` (multi-ancestry meta-analysis)

The tie is an artefact of "maximum integer in the string": for the multi-ancestry study the maximum
component happens to be its European sub-cohort.

**Why this is a defect and not a preference.** R5 fixes the LD reference panel to 1000 Genomes
**EUR unrelated**. Clumping a multi-ancestry meta-analysis against a European-only LD panel is an
LD mismatch — the LD structure used to select and de-correlate instruments would not be the LD
structure of the analysed sample. That is a broken chain, not a valid positive control, and it
would be broken regardless of which direction the result came out.

**Amendment.** R2 gains a criterion, applied before the size ranking:

> (d) the study's `INITIAL SAMPLE SIZE` must describe **European ancestry only** — it must mention
> European ancestry and must not mention any other ancestry group — so that the outcome GWAS is
> ancestry-matched to the EUR LD reference panel fixed in R5.

**Disclosure.** This amendment was made after seeing study *metadata* and before computing any
harmonization, MR, F-statistic or colocalization result. It is motivated purely by LD-panel
matching. It is not permitted to be revisited after results are seen.

---

## AMENDMENT R3a — 2026-08-26 — exposure resource changed BEFORE any re-run

**Old R3.** GTEx_v10 liver, eQTL Catalogue R8 dataset `QTD000266`, `quant_method = ge`, fixed for
all three control pairs; no per-pair tissue selection, no maximum over tissues.

**New R3.** **eQTLGen phase I cis-eQTL, whole blood, n ≈ 31,684.** Fixed for all pairs. No per-pair
selection. No maximum over resources or tissues.

**Reason.** Instrument power, and nothing else. Under old R3 all three positive controls returned
zero instruments at p < 5e-8, and none of the three genes was an eGene at FDR < 0.05 in that tissue
(HMGCR p_beta 0.675; PCSK9 0.0248; IL6R 0.448). GTEx v10 liver detects 2,603 eGenes of 20,540 genes
tested (~13 %), which is a limitation of the resource rather than of the controls.

**What did NOT change.** The instrument threshold remains **p < 5 × 10⁻⁸** (frozen by DECISION 6),
the clumping rule remains r² < 0.01 against 1000 Genomes EUR unrelated, the F ≥ 10 criterion
remains, the cis-window remains gene body ± 100 kb, and the pass criteria in R10 remain as written.
Only the source of the exposure statistics changed.

**Ordering.** This amendment is recorded before any re-run of the controls. Deliverable 1b is a
census/lookup against the new resource; it computes no MR and re-runs no control.

**Open issue created by this amendment, not resolved here.** eQTLGen phase I is distributed on
**GRCh37**; the GWAS Catalog harmonised files and the Deliverable 1 chain are **GRCh38**. R6
currently asserts both sides are GRCh38, which will no longer hold. A build-reconciliation rule
(liftover direction, tool, version, and the handling of variants that fail to lift) must be
pre-registered before the controls are re-run.

---

## R3b. LANE B EXPOSURE TISSUE SET — declared 2026-08-27, before any cardiology result exists

### Options measured, not assumed

Every sample size below is read from `data/raw/dataset_metadata_r8.tsv`; every yield figure is
computed from each dataset's own permutation output:

```bash
python src/tissue_power.py
```

| dataset | tissue | n | eGenes (BH FDR<0.05) | genes with top cis p<5e-8 | Minikel targets p<5e-8 | % of 2,517 |
|---|---|---|---|---|---|---|
| QTD000266 | GTEx v10 liver | 269 | 2,603 / 20,540 (12.7 %) | 1,711 (8.3 %) | 108 | 4.3 % |
| QTD000136 | GTEx v10 artery (coronary) | 279 | 2,632 / 22,904 (11.5 %) | 1,762 (7.7 %) | 100 | 4.0 % |
| QTD000251 | GTEx v10 heart (atrial appendage) | 473 | 4,997 / 22,258 (22.5 %) | 3,390 (15.2 %) | 232 | 9.2 % |
| QTD000256 | GTEx v10 heart (left ventricle) | 481 | 4,592 / 21,873 (21.0 %) | 3,197 (14.6 %) | 222 | 8.8 % |
| QTD000131 | GTEx v10 artery (aorta) | 492 | 5,696 / 21,888 (26.0 %) | 3,807 (17.4 %) | 268 | 10.6 % |
| QTD000141 | GTEx v10 artery (tibial) | 736 | 7,130 / 21,101 (33.8 %) | 4,941 (23.4 %) | 359 | 14.3 % |

Union of heart LV + atrial appendage + aorta = **422 targets (16.8 %)**; adding coronary artery and
liver = **476 (18.9 %)**. Compare **eQTLGen whole blood: 1,576 targets (62.6 %)**, measured under
the stricter ±100 kb window after clumping.

**These tissue figures are UPPER BOUNDS.** They use the permuted file top cis p-value over eQTL
Catalogue's ±1 Mb window, wider than our pre-registered ±100 kb, and they are not clumped. The true
Lane B yield will be lower.

GTEx v10 is the **only** cardiac source in eQTL Catalogue R8 — verified across all quantification
methods and against the upcoming-dataset table. No larger cardiac eQTL resource was found.

### The rule, fixed now

**Lane B primary exposure tissue: GTEx v10 heart left ventricle, eQTL Catalogue R8 `QTD000256`,
`quant_method = ge`, n = 481.** One tissue. No maximum over tissues.

**Lane B secondary arm, pre-corrected:** exactly four tissues — heart left ventricle `QTD000256`,
heart atrial appendage `QTD000251`, artery aorta `QTD000131`, liver `QTD000266`. Each analysed
**separately** and reported **separately**, with a Bonferroni factor of **4** applied to the
secondary arm. A target may not be promoted into the primary result because it reached significance
in a secondary tissue. This exists so a hepatic target such as PCSK9 is not structurally invisible
to a cardiac lane, while avoiding the uncorrected max-over-tissues problem Deliverable 6 exists to
catch.

Left ventricle is primary rather than aorta (which has marginally higher yield) because the Lane B
question is cardiac, and choosing a tissue by measured yield rather than by mechanism would be
selection on the very quantity the lane is powered by. Aorta remains in the corrected secondary arm.

### Power verdict, stated plainly

**Expression-level cardiac data is underpowered for Lane B.** At n ≈ 480 the cardiac tissues
instrument roughly one target in nine, against nearly two in three for whole blood. Lane B run on
GTEx cardiac expression alone will be instrument-limited, and its negative results will be
uninterpretable — absence of an instrument, not absence of an effect. This is a finding about
resource availability, not a failure of the design.

**The technically indicated route for Lane B is protein-level exposures (pQTL), not expression.**
What was verified about that route:

| resource | scale | status |
|---|---|---|
| **UKB-PPP** (Sun et al. 2023, n ≈ 54,000) | largest | **NOT in GWAS Catalog** — PMID 37794186 returns 0 studies. Synapse project `syn51364943` resolves (HTTP 200, name "UKB-PPP") but file access requires credentials. **CREDENTIAL GATE — stopped, not worked around.** |
| **deCODE** (Ferkingstad 2021, plasma proteome) | large | Listed on decode.com/summarydata/, but the download site is a JavaScript single-page app; the HTTP 200 responses are the app shell, not files. **Actual file accessibility UNVERIFIED** — not confirmed by download. |
| **Fenland** (Pietzner 2021, n = 10,708, 1,542 proteins) | medium | In GWAS Catalog, but only **3 of 1,542** studies have harmonised summary statistics (0.2 %) — effectively unavailable by the harmonised+tabix route. |
| **INTERVAL** (Sun BB 2018, n = 3,301) | small | 126 studies, **100 % harmonised**. Credential-free and range-queryable, but n only ~7x GTEx heart LV. |
| **CSF proteome** (Western 2024, n = 2,326–3,506) | small | 72.4 % harmonised; cerebrospinal fluid, not plasma — wrong compartment for cardiology. |

**Recommendation, for the user's decision, not acted on:** Lane B should use UKB-PPP plasma pQTLs
as the primary exposure, which requires supplying Synapse credentials and accepting the UK Biobank
data-use terms. If credentials are unavailable, the fallback is INTERVAL (n = 3,301, credential-
free, harmonised) — ~7x better powered than GTEx heart but still two orders of magnitude below
eQTLGen. The GTEx cardiac tissue set fixed above stands as the pre-registered expression arm either
way, so Lane B has a declared rule regardless of how the credential question resolves.

---

## R12. CIS-MR ESTIMATOR RULE — FROZEN 2026-08-27, before any MR output exists

This section replaces R7 for all cis-MR from this point forward. Written before any estimate was
computed; not revised afterwards.

### R12.1 Why the previous rule was inadequate

Deliverable 1b showed `r² < 0.01` returning a **median of 3 and a maximum of 30** "independent"
instruments at a single cis locus. The computation was verified correct (for CBS, max pairwise r²
among kept SNPs = 0.00967 while the lead SNP had r² = 1.000 with six neighbours). The problem is
interpretive, not computational: variants within one ~200 kb cis window share a single regulatory
mechanism and one ancestral haplotype structure. Thirty mutually near-uncorrelated SNPs at one
locus are **not** thirty independent natural experiments, and treating them as such in a standard
IVW makes the standard error too small — a false-precision failure, and exactly what Deliverable 6
exists to catch.

### R12.2 Primary estimator — FIXED

**The primary cis-MR estimator is the lead-SNP Wald ratio.**

* Instrument: the single variant with the smallest eQTL p-value inside the cis-window that survives
  harmonization and liftover.
* β_MR = β_outcome / β_exposure
* SE_MR = |SE_outcome / β_exposure| (first-order delta method; the exposure SE is not propagated,
  which is conservative in the usual direction for cis-MR with a strong instrument)
* z = β_MR / SE_MR, two-sided p from the normal distribution.

### R12.3 Why lead-SNP Wald ratio and not LD-aware IVW

Both were considered. LD-aware IVW (generalised IVW using the full instrument correlation matrix Ω,
β = (G'Ω⁻¹G)⁻¹G'Ω⁻¹β_out) is statistically more efficient **when Ω is accurately known**. It is
rejected as primary for three reasons that are properties of this project, not general objections:

1. **Ω is estimated from 503 EUR samples.** With 30 instruments, a 30x30 correlation matrix
   estimated from 503 individuals is near-singular; inverting it is numerically unstable and
   produces exactly the false-precision failure we are trying to avoid. The instability grows with
   instrument count — worst precisely where LD-aware IVW would appear to gain most.
2. **The LD reference is not the analysis sample.** eQTLGen is a 31,684-sample meta-analysis over
   dozens of cohorts; 1000G EUR is a proxy. Errors in Ω propagate into the point estimate under
   LD-aware IVW, whereas the lead-SNP Wald ratio never uses Ω.
3. **One cis locus supports one instrument's worth of evidence.** The lead-SNP Wald ratio makes
   that assumption explicit rather than burying it in a correlation matrix.

Cost accepted and stated: the lead-SNP Wald ratio is less efficient and discards information when a
locus genuinely contains multiple independent regulatory signals.

**LD-aware IVW is retained as a pre-declared secondary sensitivity analysis** where ≥2 instruments
exist, using the same 1000G EUR panel, with the condition number of Ω reported alongside every
estimate and the estimate flagged unstable when that condition number exceeds 1,000. It never
replaces the primary.

### R12.4 Behaviour by instrument count — FIXED

| instruments after clumping | primary | secondary / sensitivity |
|---|---|---|
| **1** | lead-SNP Wald ratio | none available; reported as such |
| **2** | lead-SNP Wald ratio | LD-aware IVW; no MR-Egger, no weighted median |
| **>=3** | lead-SNP Wald ratio | LD-aware IVW; MR-Egger; weighted median |
| **0** | no estimate; target reported as *not instrumented*, never as null | — |

### R12.5 Sensitivity analyses and their explicit unavailability

Where ≥3 instruments exist: **MR-Egger** (intercept as directional-pleiotropy test) and **weighted
median**, reported with the LD-aware IVW estimate and the condition number of Ω.

**Explicit statement, as required:** in Deliverable 1b, **460 of 1,576 instrumented targets (29.2 %)
have exactly one instrument; 1,116 (70.8 %) have ≥2; 818 (51.9 %) have ≥3.** Therefore **MR-Egger
and weighted median are unavailable for 758 of 1,576 targets — 48.1 %.** For nearly half the
analysable targets there is no pleiotropy-robust sensitivity analysis at all, and no downstream
modelling changes that. Any claim that the pipeline is "pleiotropy-checked" must be restricted to
the 51.9 % where the check exists; the 48.1 % must be reported as unchecked, never as passing.

Applicable at every instrument count: **colocalization** (R8, unchanged) is reported for every pair,
since it uses all variants in the region and does not depend on instrument count.

### R12.6 Cis-window and clumping — restated unchanged

* cis-window: **gene body ± 100,000 bp** (R1, unchanged)
* instrument threshold: **p < 5 x 10^-8** (DECISION 6, frozen)
* clumping: greedy, lowest eQTL p first, drop r² ≥ **0.01**, LD from 1000 Genomes **EUR** (R5,
  unchanged). Clumping still runs: it defines the instrument count that selects the row of the
  R12.4 table and supplies the instrument set for the secondary analyses.
* instrument strength: F = (β_eQTL/SE_eQTL)² for the lead SNP; criterion **F ≥ 10** (R5, unchanged)

---

## R9a. Direction expectations under the changed exposure resource — recorded 2026-08-27, before Part C ran

R9 fixed the expected signs while the exposure was **liver** expression. Under AMENDMENT R3a the
exposure is **eQTLGen whole blood**. The units of the MR estimate therefore change: it is now
"change in trait per +1 SD of gene expression **in whole blood**".

**The sign expectations are unchanged**, because they are anchored in pharmacology rather than in
the measured tissue:

| pair | expected sign | anchor |
|---|---|---|
| HMGCR → LDL-C | **POSITIVE** | statins inhibit HMGCR and lower LDL-C |
| IL6R → CRP | **NEGATIVE** | the lead cis signal at *IL6R* is expected to tag rs2228145, whose higher-expression allele increases receptor shedding, reduces classical hepatic signalling, and lowers CRP |

**Risk recorded in advance, not after the fact:** a cis-eQTL measured in whole blood can be
*discordant in sign* with the same variant's effect in the tissue where the drug acts. If a control
returns the wrong sign, tissue discordance is a candidate explanation — but under R10 it is still
recorded as a **direction failure and a BLOCKER**, and it will not be reinterpreted post hoc into a
pass. Stating the possibility now is what makes it not an excuse later.

PCSK9 is retired as a control by DECISION 8 and is not replaced.

---

## R13. PART C — RULE-BASED DISTRIBUTIONAL VALIDATION — FROZEN 2026-08-30, before any concordance was computed

Named positive controls are retired (DECISION 13). Validation is now the concordance between
cis-MR direction and the direction implied by an approved drug's mechanism, across every
qualifying target-indication pair.

### R13.1 Mechanism resource

**Open Targets Platform release `26.06`**, pinned explicitly rather than `latest`, datasets
`drug_mechanism_of_action` (6,500 rows) and `drug_molecule` (22,407 rows).

*Carrying forward OPEN ITEM 4:* release 26.06's top-level `manifest.json` reports
`"result": "failure"`. Both dataset directories used here carry `_SUCCESS` markers, which is the
per-dataset completion signal. This is recorded as a known, unresolved caveat on the mechanism
annotations, not silently accepted.

### R13.2 Extraction rule

1. **Approved drug** = `drug_molecule.maximumClinicalStage == "APPROVAL"`.
2. Keep a mechanism row only if: at least one of its `chemblIds` is an approved drug; and
   `targetType == "single protein"`; and `targets` contains **exactly one** Ensembl gene ID.
   Multi-target mechanisms are **excluded as ambiguous**, not apportioned.
3. **Direction of the drug on target activity:**
   * **DOWN** — INHIBITOR, ANTAGONIST, BLOCKER, NEGATIVE MODULATOR, NEGATIVE ALLOSTERIC
     MODULATOR, ALLOSTERIC ANTAGONIST, INVERSE AGONIST, ANTISENSE INHIBITOR, RNAI INHIBITOR,
     DEGRADER, DISRUPTING AGENT, GENE EDITING NEGATIVE MODULATOR
   * **UP** — AGONIST, PARTIAL AGONIST, ACTIVATOR, POSITIVE MODULATOR, POSITIVE ALLOSTERIC
     MODULATOR, OPENER, STABILISER
   * **EXCLUDED as ambiguous** — MODULATOR, BINDING AGENT, OTHER, SUBSTRATE, EXOGENOUS PROTEIN,
     EXOGENOUS GENE, VACCINE ANTIGEN, RELEASING AGENT, CROSS-LINKING AGENT, HYDROLYTIC ENZYME,
     PROTEOLYTIC ENZYME
4. **Target-level consistency:** a target qualifies only if **every** approved-drug mechanism for
   it agrees on direction. A target carrying both UP and DOWN mechanisms is **excluded as
   ambiguous**. No majority vote, no tie-breaking.

### R13.3 Concordance rule

An approved drug that acts on target T and treats disease D implies a direction for T's causal
effect on D:

* drug direction **DOWN** ⇒ lowering T improves D ⇒ **raising T worsens D** ⇒ expected **MR β > 0**
* drug direction **UP** ⇒ raising T improves D ⇒ expected **MR β < 0**

A pair is **CONCORDANT** when the sign of the cis-MR estimate matches the expectation above.

### R13.4 Indication → outcome GWAS — a Part-C-only expedient, explicitly NOT the Deliverable 3 mapping

Part C needs an outcome GWAS per indication, and the Deliverable 3 indication-mapping does not yet
exist. Rather than pre-empt D3, Part C uses the most conservative mechanical rule available:

> exact case-insensitive string match between Minikel `indication_mesh_term` and a GWAS Catalog
> `MAPPED_TRAIT`; among matching studies require harmonised summary statistics and a `.tbi`;
> select the European-ancestry-only study (R2a) with the largest parsed sample size.

This is **high precision, low recall, and unvalidated**. It is a scaffold for Part C only, is not
reused downstream, and will be superseded by the Deliverable 3 mapping. Its coverage is reported
as a limitation, not as a result.

### R13.5 Analysis and reporting, fixed in advance

* cis-MR by the frozen R12.2 rule: lead-SNP Wald ratio, p < 5e-8, ±100 kb, r² < 0.01, F ≥ 10.
* Report: N qualifying pairs; concordant fraction; **two-sided exact binomial test against 50 %**.
* Stratify the concordant fraction by instrument count (1, 2, ≥3) and by F.
* List every discordant pair, classified where possible as **CLASS B modality mismatch**
  (splicing, shedding, protein activity, localisation) versus **unexplained**.
* **Quantify the reliability of the mechanism extraction itself**: the fraction of mechanism rows
  excluded as ambiguous, and the fraction of targets excluded for internal direction conflict.
  If that fraction is large, the concordance rate is correspondingly weak evidence and must be
  reported that way.

No threshold, direction map, or exclusion in R13 may be changed after concordance is computed.

---

## R14. PART 3 — DEGRADATION MEASUREMENT: LOCUS SELECTION — FROZEN 2026-08-30, before any PP4 was computed for these loci

### R14.1 Selection rule (mechanical, no hand-picking)

Start from `deliverables/d1f/partc_results.tsv` — the 571 target-indication pairs that already
produced a cis-MR estimate under the frozen Part C rules. Keep a pair only if:

1. `n_overlap >= 100` — enough harmonized variants for a stable LD matrix and a meaningful
   colocalization;
2. the target has >= 1 eQTLGen cis-eQTL at p < 5e-8 inside gene body +-100 kb (true by
   construction of the Part C set);
3. **one locus per target**: where a target appears with several indications, keep the pair with
   the **largest `n_overlap`**. This is deterministic and involves no colocalization or MR result,
   so it cannot select on the quantity being measured. It also makes the paired tests in R14.3/
   R14.4 independent across targets.

If more than **60** loci survive, take a random sample of 60 with **seed 20260830**. The cap exists
because each locus costs one outcome-GWAS region query, one GTEx-blood region query and one 1000G
LD region query; it is a compute budget, not a filter on results.

### R14.2 Per-locus quantities reported

eQTL lead z; GWAS lead z; r2 between the two leads (1000G EUR); distance between leads (bp);
lBF gap `s1+s2-s12`; coloc.abf PP3 and PP4; coloc.susie PP4; and the **residual-z RMSE for both
traits** — the DENTIST statistic (Chen et al. 2021, Nat Commun 12:7117), reimplemented here and
cited as theirs.

### R14.3 Degradation test

Does the abf-vs-susie PP4 discrepancy increase with eQTL lead z, and with lead-to-lead r2 below 1?
Reported as Spearman correlation with p-value, plus a regression of
`log10(PP4_susie) - log10(PP4_abf)` on eQTL lead z and on r2. Effect size and test, not narrative.

### R14.4 Head-to-head across exposure power — the primary comparison

For every selected target that is also present in **GTEx whole blood (QTD000356, n = 853)**, run
**coloc.abf** separately with each exposure resource against the same outcome GWAS, same window,
same harmonization:

* GTEx whole blood, n = 853
* eQTLGen, n = 31,684

Paired per target; **Wilcoxon signed-rank** on `log10(PP4)` plus the count of targets where PP4 is
lower in the more powerful resource, with an exact binomial test.

**PRE-REGISTERED REFUTATION (Amendment 1, already in PROTOCOL.md): if PP4 INCREASES with exposure
power, the Class C power-dependence claim is REFUTED and reported as a refutation. Class C is
dropped from the paper. Classes A and B are unaffected.**

### R14.5 Intermediate power point — UNAVAILABLE, determined before running

Amendment 2 requires a third power point at n ~ 2,000-10,000, chosen by sample size alone.
**No such dataset is obtainable.** eQTL Catalogue R8 lists exactly two `ge` datasets in that range —
INTERVAL_RNA `QTD001000` (n = 4,729) and INTERVAL_RNA_WGS `QTD001003` (n = 2,892) — and **both are
absent from the FTP** (`QTS000053`, `QTS000054` return HTTP 404; the FTP holds QTS000001-QTS000042
only). Cross-checking all 466 R8 `ge` datasets against the FTP listing gives **0 datasets with
n in [2,000, 10,000]**; the largest FTP-available `ge` dataset is n = 954 (OneK1K CD4+ TCM).

**Therefore the non-monotone hypothesis (PP4 rising with power to ~n=2,000 then falling by
n~30,000) is UNVERIFIABLE with available data.** Two power points give a slope, not a shape, and no
curve will be claimed from them. The head-to-head tests the slope only.

---

## R15. RE-SPECIFIED DESCRIPTIVE LD/COLOC MEASUREMENT — declared 2026-08-31, before running

**Status: DESCRIPTIVE, not a hypothesis test.** Class C was refuted on its pre-registered
criterion (PROTOCOL.md) and is dropped, so nothing here is confirmatory and re-specifying the
locus rule cannot contaminate a live hypothesis. R14.1's omission of an outcome-significance
criterion is corrected here.

**Rule (mechanical, declared before applying):** from the 571 Part C pairs that produced a cis-MR
estimate, keep a pair if ALL hold:

1. `n_overlap >= 100`;
2. **the outcome GWAS has >= 1 variant with p < 5e-8 (|z| > 5.4513) inside the merged cis-window** —
   the criterion R14.1 omitted. Colocalization is only interpretable when BOTH traits are
   associated; the exposure side is genome-wide significant by construction;
3. one locus per target, keeping the largest `n_overlap` (deterministic; no coloc or MR quantity
   involved).

**No cap.** Every qualifying locus is run to completion, so nothing is left UNVERIFIED (FIX 3).

**Reported:** DENTIST residual-z RMSE (Chen et al. 2021, Nat Commun 12:7117) for both traits, and
credible-set availability, under this rule — **side by side with the R14.1 numbers, which are
labelled as measured under a flawed selection rule.**

**The R14.1 60-locus set is superseded.** Of its 60 loci, 28 completed and 32 did not (the 1000
Genomes LD endpoint stalled 15-20 min on individual loci). The 32 are **excluded as superseded**,
not left pending: they were selected by a rule now known to be inadequate, so completing them would
add loci that mostly lack an outcome signal. Final N for the flawed-rule comparison is **28**.

**Colocalization remains excluded from the atlas as an evidence gate (Gate D closed).** These
numbers describe method availability, not evidence for any target.
