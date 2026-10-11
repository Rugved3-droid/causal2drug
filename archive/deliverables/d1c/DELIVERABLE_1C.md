# Part A/B/C — estimator freeze, build alignment, control re-run — **BLOCKER**

Run 2026-08-27 UTC. Rules frozen in `deliverables/d1/PREREGISTRATION.md` §R12 (estimator), §R3b
(Lane B tissues), §R9a (direction expectations) — all written before any output existed.

## VERDICT

| control | MR β | SE | p | direction | expected | instruments | F | verdict |
|---|---|---|---|---|---|---|---|---|
| **HMGCR → LDL-C** | **+0.3891** | 0.01204 | 4.24e-229 | POSITIVE | POSITIVE | 1 | 224.0 | **PASS** |
| **IL6R → CRP** | **+0.3762** | 0.02309 | 1.16e-59 | POSITIVE | **NEGATIVE** | 2 | 134.5 | **FAIL** |

**IL6R → CRP returns the wrong sign under the pre-registration. Reported as BLOCKER. Stopping.**
Nothing was tuned; no threshold, window, tissue or estimator was changed after seeing this.

```bash
cd /d/causal2drug
python src/test_harmonize.py        # Part B2 branch tests
python src/audit_prefix.py          # Part B3 pre/post-fix audit
python src/control_chain_v2.py      # Part C, both controls
```

---

## PART A — estimator rule, frozen before any MR output

Full text in `PREREGISTRATION.md` §R12. Summary:

* **Primary estimator: lead-SNP Wald ratio.** β = β_out/β_exp, SE = |SE_out/β_exp|.
* **Correlated instruments are NOT combined.** LD-aware IVW was considered and rejected as primary
  because Ω would be a 30×30 correlation matrix estimated from 503 EUR samples, the LD reference is
  a proxy for a 31,684-sample meta-analysis rather than the analysis sample, and one cis locus
  supports one instrument's worth of evidence. LD-aware IVW is retained as a **secondary**
  sensitivity analysis with the condition number of Ω reported and flagged unstable above 1,000.
* **By instrument count:** 1 → Wald only; 2 → Wald + LD-aware IVW; ≥3 → Wald + IVW + MR-Egger +
  weighted median; 0 → no estimate, target reported as *not instrumented*, never as null.
* **Explicit unavailability:** of 1,576 instrumented targets, 460 (29.2 %) have one instrument and
  818 (51.9 %) have ≥3. **MR-Egger and weighted median are unavailable for 758 targets — 48.1 %.**
  Any claim of "pleiotropy-checked" must be restricted to the 51.9 % where the check exists.
* Window ±100 kb, p < 5e-8, clumping r² < 0.01, F ≥ 10 — all unchanged.

---

## PART B1 — GRCh37 → GRCh38 liftover

Chain: UCSC `hg19ToHg38.over.chain.gz`, 227,698 bytes, sha256
`5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d`, retrieved
2026-08-27T23:26:55Z. Implementation: `pyliftover` (pure Python).

**Validation against an independent source:** HMGCR GRCh37 start 74,632,154 lifts to GRCh38
75,336,329, which equals the Ensembl GRCh38 gene start used in Deliverable 1 — exact match.

| pair | variants in | lifted | failed | multi-mapping | negative strand | failure rate |
|---|---|---|---|---|---|---|
| HMGCR → LDL-C | 277 | 277 | 0 | 0 | 0 | **0.0 %** |
| IL6R → CRP | 226 | 226 | 0 | 0 | 0 | **0.0 %** |

**Handling of failures (declared, though none occurred here):** a variant that fails to lift is
**dropped and counted**; it is never carried forward on its GRCh37 coordinate. Multi-mapping
variants take the first chain segment and are counted separately. Both counts appear in each pair's
provenance JSON.

A 0 % failure rate is expected for common SNPs in gene-dense autosomal regions and should **not** be
generalised: genome-wide liftover failure is typically 0.1–1 %, concentrated in repetitive and
assembly-corrected regions. **UNVERIFIED at genome-wide scale** — measured here only across 503
variants in two cis-windows.

## PART B2 — the never-executed harmonization branches

In Deliverable 1 all 1,838 merged variants resolved `aligned_same`, so the allele-flip and
strand-flip branches had never run. Harmonization is now isolated in `src/harmonize.py` and
exercised by `src/test_harmonize.py`.

```
case                                       exp    out    beta_in  beta_out  status                             result
(0) control: no flip needed                G/A    A/G    +0.30    +0.30     aligned_same                       PASS
(a) ALLELE FLIP                            G/A    G/A    +0.30    -0.30     aligned_flipped                    PASS
(b) STRAND FLIP                            G/A    T/C    +0.30    +0.30     strand_flipped                     PASS
(c) STRAND FLIP + ALLELE FLIP              G/A    C/T    +0.30    -0.30     strand_flipped_and_allele_flipped  PASS
(d) genuine mismatch                       G/A    C/G    +0.30    dropped   allele_mismatch                    PASS
(e) palindromic, resolvable by frequency   A/T    T/A    +0.30    +0.30     aligned_same                       PASS
(f) palindromic, AF too close to 0.5       C/G    G/C    +0.30    dropped   ambiguous_unresolved               PASS
(g) palindromic, AF discordant             A/T    T/A    +0.30    dropped   ambiguous_unresolved               PASS
(h) palindromic allele-flip, freq-resolved A/T    A/T    +0.30    -0.30     aligned_flipped                    PASS
(i) indel must be rejected                 G/GA   GA/G   +0.30    dropped   not_snp                            PASS
RESULT: 10/10 cases pass
```

Both sign-inverting branches (a) and (c) invert β exactly; the strand-only branch (b) preserves it.
**These branches did execute on real data in Part C**: the IL6R merge resolved 76 variants as
`aligned_flipped`, the first time that branch has run outside a test.

## PART B3 — pre-fix / post-fix audit

The defect: the tabix linear index was used as an **upper** bound as well as a lower bound. It is
valid only as a lower bound, so regions could be silently truncated.

**Unaffected by construction** (never passed through tabix): all `pp.tsv` counts and checksums; the
`.tbi` availability probes (HEAD requests only); the eQTLGen scans, tested-gene universe, mapping
census and therapeutic-area tables (gzip streaming); GTF gene coordinates; all D1b clumping results
(re-run from scratch **after** the fix).

**Computed before the fix and therefore re-run** — `python src/audit_prefix.py`:

| query | before fix | after fix | match |
|---|---|---|---|
| D0 GWAS GCST90239658 chr5 HMGCR window | 3,323 | 3,323 | YES |
| D0 eQTL QTD000266 chr5 HMGCR window (all genes) | 8,990 | 8,990 | YES |
| D0 1000G GRCh38 chr5 HMGCR window | 5,644 | 5,644 | YES |
| D1 HMGCR GWAS | 3,323 | 3,323 | YES |
| D1 HMGCR eQTL | 592 | 592 | YES |
| D1 PCSK9 GWAS | 3,925 | 3,925 | YES |
| D1 PCSK9 eQTL | 846 | 846 | YES |
| D1 IL6R GWAS | 581 | 581 | YES |
| D1 IL6R eQTL | 839 | 839 | YES |

**9/9 reproduce unchanged.** I previously flagged D1's extract row counts as UNVERIFIED because of
this bug; that flag is now **withdrawn** — they are verified correct. The bug did not bite on those
particular windows (each resolved to a single index span; the 1000G window resolved to four spans
and still matched). It bit elsewhere, in D1b, where it was caught and corrected.

`tabix_http.py` now carries all three fixes (no upper bound, per-chunk sliced fetching with early
stop, carry-over across BGZF block boundaries), so the shared module and `d1b_clump.py` no longer
diverge.

---

## PART C — controls re-run

Exposure: eQTLGen phase I whole blood (GRCh37). β and SE derived from Z by the SMR approximation
(Zhu et al. 2016) using eQTLGen's **own** allele frequencies
(`eqtlgen_SNP_AF.txt.gz`, 240,045,342 B, sha256 `edffd43094598d08…`); all 503 required rsIDs were
found. Under this parameterisation β/SE = z exactly, so **F = z²** and instrument strength does not
depend on the allele frequency.

Outcome: GWAS Catalog harmonised (GRCh38), selected by R2+R2a — LDL `GCST90239658` (Graham 2021,
EUR N=1,320,016); CRP `GCST90179146` (Koskeridis 2022, N=575,531).

### HMGCR → LDL-C — PASS

| quantity | value |
|---|---|
| eQTLGen hits p<5e-8 in ±100 kb | 277 |
| independent instruments after clumping (r²<0.01) | **1** (rs6453133) |
| liftover | 277/277, 0 failures |
| **overlapping variants after harmonization + liftover** | **276** |
| harmonization | 276 `aligned_same`, 1 `ambiguous_unresolved` (dropped) |
| lead SNP | rs6453133, p_exposure = 1.21e-50 |
| β_exposure / SE | +0.12776 / 0.008536 |
| β_outcome (harmonized) / SE | +0.049710 / 0.0015383 |
| **MR β / SE / p** | **+0.38908 / 0.012040 / 4.24e-229** |
| **direction** | **POSITIVE** (expected POSITIVE) |
| **F** | **224.0** |
| coloc PP0 / PP1 / PP2 / PP3 / PP4 | 0.0 / 0.0 / 5.77e-44 / **1.000** / 1.00e-11 |

Higher HMGCR expression → higher LDL-C, the direction implied by statin pharmacology. All four
pass criteria met.

### IL6R → CRP — **FAIL (direction)**

| quantity | value |
|---|---|
| eQTLGen hits p<5e-8 in ±100 kb | 226 |
| independent instruments after clumping | **2** (rs10908839, rs58229059) |
| liftover | 226/226, 0 failures |
| **overlapping variants after harmonization + liftover** | **167** |
| harmonization | 91 `aligned_same`, **76 `aligned_flipped`**, 59 `no_position_match` (dropped) |
| lead SNP **used** | rs2229238, p_exposure = 4.27e-31 |
| β_exposure / SE | +0.11175 / 0.009636 |
| β_outcome (harmonized) / SE | +0.042036 / 0.0025804 |
| **MR β / SE / p** | **+0.37616 / 0.023091 / 1.16e-59** |
| **direction** | **POSITIVE — expected NEGATIVE** |
| **F** | **134.5** |
| coloc PP0 / PP1 / PP2 / PP3 / PP4 | 0.0 / 0.0 / 5.90e-25 / **1.000** / 7.96e-17 |

### Is the wrong sign a bug? What was checked

1. **The sign machinery is correct.** 10/10 synthetic branch tests pass, including both
   sign-inverting cases. HMGCR returned the *correct* sign through the identical code path.
2. **It is not a single-SNP artefact.** Among the 102 merged variants with p_exposure < 1e-10,
   **94 give a positive Wald ratio and 8 negative**. The positive direction is consistent across
   the locus.
3. **The instrument set cannot express the mechanism the expectation was based on.** R9 predicted
   NEGATIVE via rs2228145 (Asp358Ala) increasing receptor shedding. **rs2228145 is absent from
   eQTLGen's significant cis-eQTL set for IL6R entirely**, and rs4129267 is present only at
   p = 1.11e-05, above the frozen 5e-8 threshold. The pre-registered mechanism is therefore not
   instrumented.
4. **The outcome GWAS is sparse and lacks the true lead.** `GCST90179146` has only **63 variants in
   a 28 kb window**; the clumping lead **rs10908839 is absent from it altogether** (by position and
   by rsID), so MR fell back to rs2229238. 59 of 226 exposure variants (26 %) found no position
   match at all.

**This does not convert the verdict.** R9a states in advance that a wrong sign is recorded as a
direction failure and will not be retroactively reinterpreted. For completeness and not as a
defence: at *total mRNA* level a positive sign is arguably what tocilizumab pharmacology implies
(less IL6R message → less receptor → less signalling → less CRP), and my pre-registered NEGATIVE
rested on a shedding mechanism that the available instruments cannot capture. That is an argument
that my expectation was miscast — but it was made after seeing the result, so it changes nothing
here. **FAIL stands.**

### A second concern, on the control that PASSED

**Colocalization returns PP3 ≈ 1.000 for both pairs, including HMGCR.** PP3 is "both traits
associated, *different* causal variants"; PP4 (shared causal variant) is ~1e-11 and ~1e-16. For a
relationship as well established as HMGCR→LDL-C, near-zero PP4 is not what a correctly specified
coloc should return. Candidate explanations — none tested, none assumed:

* multiple independent causal variants in the window violating coloc's single-variant assumption;
* the very large GWAS N driving the ABF strongly toward whichever variant is the GWAS lead when it
  differs from the eQTL lead;
* a defect in the coloc implementation or its priors.

The coloc numbers are reported as computed and should be treated as **not yet trustworthy** until
this is resolved. It does not affect the MR verdicts, which do not use coloc.

---

## Intermediates saved

```
deliverables/d1c/
  DELIVERABLE_1C.md
  d1c_summary.json                      every number above
  prefix_audit.json                     Part B3 audit
  HMGCR_LDL_harmonized_merge.tsv        276 rows x 25 cols
  IL6R_CRP_harmonized_merge.tsv         167 rows x 25 cols
  HMGCR_LDL_provenance.json
  IL6R_CRP_provenance.json
src/harmonize.py                        harmonization, isolated for testing
src/test_harmonize.py                   10 synthetic branch cases
src/audit_prefix.py                     pre/post-fix recomputation
src/control_chain_v2.py                 Part C chain
src/tissue_power.py                     Lane B tissue power measurement
```

Merge tables carry, per variant: `rsid, chrom, pos37, pos38, exp_other, exp_assessed,
gwas_effect, gwas_other, harmonization, sign_applied, ambiguous, effect_allele_harmonized,
other_allele_harmonized, beta_exp, se_exp, p_exp, z_exp, n_exp, af_assessed_exp, F,
beta_out_raw, beta_out_harm, se_out, p_out, af_effect_out` — both sides side by side, with
GRCh37 and GRCh38 positions so the liftover is auditable per row.

**Deliverable 2 has not been started.**
