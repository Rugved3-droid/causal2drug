# Part 2 — larger LD reference at HMGCR: **the panel is not the problem**

```bash
cd /d/causal2drug && python src/part2_ukb_ld.py
```

Single-shot test, pre-declared purpose (Amendment 3): answer the reviewer objection that
1000G EUR n = 525 is simply too small, with one number. One reference, one result.

## The reference

| | |
|---|---|
| resource | **UK Biobank LD matrices** released with PolyFun (Weissbrod et al.) |
| access route | `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/` — public S3, no credentials |
| block | `chr5_74000001_77000001` (3 Mb), `.npz` 1,063,183,424 B + `.gz` SNP list 189,955 B (sha256 `75ff611cd9d45c55…`) |
| n | **337,000** unrelated British-ancestry individuals |
| ancestry | British (UK Biobank) — ancestry-matched to the EUR-only GWAS selected under R2a |
| build | **GRCh37** (matched to our stored `pos37`, so no liftover in this test) |
| block contents | 20,921 variants; **all 276** of our HMGCR variants matched, 38,226 off-diagonal LD entries |

The `.npz` holds a 218,854,576-entry sparse COO matrix (~2.6 GB materialised) against ~1.9 GB
free RAM, so members were extracted to disk and memory-mapped, then filtered in chunks.

## Result

| quantity | 1000G EUR n=525 | **UKB n=337,000** |
|---|---|---|
| **DENTIST residual-z RMSE, GWAS side** | 7.80 | **7.119** |
| DENTIST residual-z RMSE, eQTL side | 1.33 | 1.626 |
| coloc.abf PP3 / PP4 | 1 / 1.005e-11 | 1 / 1.005e-11 |
| coloc.susie best PP4 | 0.1482 | 0.9957 — **ARTEFACT, see below** |
| credible sets (exposure / outcome) | 1 / 5 | 3 / 5 |
| **any credible set tops at rs12916 (z = 49.59)** | **NO** | **NO** |

*(coloc.abf does not use LD, so its identical value across references is a correctness check,
not a finding.)*

The residual-z RMSE is the **DENTIST statistic (Chen et al. 2021, Nat Commun 12:7117)**,
reimplemented here and cited as theirs.

## The headline number

**A 640× larger, ancestry-matched LD reference moved the GWAS-side residual-z RMSE from 7.80
to 7.12.** Expected under a well-specified reference is ≈ 1. The reference panel is **not**
the binding constraint.

## The PP4 = 0.9957 is an artefact and is excluded

It looks like recovery above the 0.8 threshold. It is not, and the diagnostic that caught it
is the same one that caught the earlier spurious PP4 = 1.000:

| exposure CS top | outcome CS top | PP4 | same variant? |
|---|---|---|---|
| **rs1544755** | **rs1544755** | **0.9957** | **YES — artefact signature** |
| rs6453133 | rs1544755 | 0.3025 | no |
| rs184692007 | rs1544755 | 0.0176 | no |

**rs1544755 is neither trait's lead** — the eQTL lead is rs6453133 (z = 14.97), the GWAS lead
is rs12916 (z = 49.59). Both traits' SuSiE fits converged a spurious component onto the same
third variant, and colocalizing a component with itself returns PP4 ≈ 1 by construction. No
credible set tops at rs12916 under either reference.

**Excluded from every table and figure**, on the same basis as the withdrawn HMGCR and IL6R
PP4 = 1.000 values. Recorded here only so the record shows it was produced and withdrawn.

## What this means — and it is the stronger finding

Amendment 3 anticipated this: *"If Kanai's mechanism (meta-analysis heterogeneity, not panel
size) is the dominant cause at our locus, a bigger panel may not help — that is itself the
reportable finding, and it is the stronger one."* That is what happened.

The outcome GWAS is **Graham SE et al. 2021 GLGC LDL-C, n = 1,320,016** — a meta-analysis of
many cohorts differing in genotyping array, imputation panel and sample size. Its effective LD
is not the LD of *any single reference population*, so no reference panel, at any size, can
reproduce it. This is precisely the regime **Kanai et al. 2022 (Cell Genomics 2(12):100210,
SLALOM)** characterises when it states that LD reference mismatch is "all but guaranteed in
the case of GWAS meta-analyses", and finds **67 % of GBMI meta-analysis loci suspicious**.

Our locus behaves exactly as their framework predicts, and the 7.80 → 7.12 result is the
direct empirical demonstration that panel *size* is the wrong lever.

**The reviewer objection is answered with a number: no, n = 525 was not the problem.**

## Limits of this test, stated

* **One locus.** HMGCR only. Whether the ~7 residual RMSE is typical is measured at scale in
  Part 3, not here.
* **One reference.** Per Amendment 3, no second reference was tried.
* **Not tested:** in-sample LD from the GWAS cohorts themselves, which is the only thing that
  would fully resolve the mismatch and is not publicly available for GLGC. **UNVERIFIED**
  whether in-sample LD would restore PP4.
* The eQTL-side RMSE rose slightly (1.33 → 1.63). eQTLGen is itself a 31,684-sample
  meta-analysis of dozens of cohorts, so it carries the same heterogeneity problem in milder
  form; a British-only reference is not obviously better matched to it than 1000G EUR.
