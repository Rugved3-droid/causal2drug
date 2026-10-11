# Deliverable 1 — positive controls — **BLOCKER**

Run 2026-08-26 UTC under `deliverables/d1/PREREGISTRATION.md`, written before any result existed.

## VERDICT: all three controls FAIL. Nothing was tuned. Stopping here.

| pair | eGene FDR<0.05 | instruments p<5e-8 | F ≥ 10 | MR | direction | VERDICT |
|---|---|---|---|---|---|---|
| HMGCR → LDL-C | **NO** (p_beta 0.675) | **0** | n/a | not computable | n/a | **FAIL** |
| PCSK9 → LDL-C | **NO** (p_beta 0.0248) | **0** | n/a | not computable | n/a | **FAIL** |
| IL6R → CRP | **NO** (p_beta 0.448) | **0** | n/a | not computable | n/a | **FAIL** |

All three fail at the **exposure** side — the cis-eQTL — in the pre-registered tissue. No MR
estimate, standard error, p-value, direction or F statistic exists for any pair, because zero
variants reached the instrument threshold. Reporting an MR number here would require changing a
pre-registered parameter, which the pass criteria forbid.

Reproduce everything:

```bash
cd /d/causal2drug
python src/control_chain.py          # the whole chain, all three pairs
```

---

## 1. Where exactly it fails

The eGene gate (R4), computed from the dataset's own permutation file with BH FDR across all
20,540 genes (BH threshold on `p_beta` = 0.00632822; 2,603 eGenes in this tissue):

```bash
python -c "
import gzip,csv
rows=list(csv.DictReader(gzip.open(r'D:\causal2drug\data\raw\QTD000266.permuted.tsv.gz','rt'),delimiter='\t'))
v=sorted((float(r['p_beta']),r['molecular_trait_id']) for r in rows); m=len(v); thr=0
for i,(p,g) in enumerate(v,1):
    if p<=0.05*i/m: thr=p
print('BH thr',thr,'n_eGenes',sum(1 for p,g in v if p<=thr))"
```

| gene | cis variants tested | top nominal p | `p_beta` | eGene at FDR<0.05? |
|---|---|---|---|---|
| HMGCR `ENSG00000113161` | 6,674 | 1.09 × 10⁻³ | 0.675 | **NO** |
| PCSK9 `ENSG00000169174` | 7,103 | 1.03 × 10⁻⁵ | 0.0248 | **NO** |
| IL6R `ENSG00000160712` | 4,182 | 5.97 × 10⁻⁴ | 0.448 | **NO** |

And the instrument gate (R5), computed independently from the nominal summary statistics inside the
±100 kb window after harmonization:

| pair | harmonized variants | with p_GWAS < 5e-8 | with p_eQTL < 5e-8 | best eQTL p | best eQTL variant | F of that variant |
|---|---|---|---|---|---|---|
| HMGCR → LDL-C | 529 | **367** | **0** | 1.33 × 10⁻² | `chr5_75268665_G_A` rs143528459 | 6.25 |
| PCSK9 → LDL-C | 763 | **446** | **0** | 4.57 × 10⁻⁶ | `chr1_55054539_G_A` rs1165287 | 22.26 |
| IL6R → CRP | 546 | **416** | **0** | 6.44 × 10⁻³ | `chr1_154408010_G_A` rs6427641 | 7.59 |

PCSK9's best cis-eQTL has F = 22.26, comfortably above the F ≥ 10 criterion — but its p-value is
4.57 × 10⁻⁶, which does not reach the pre-registered 5 × 10⁻⁸. Relaxing that threshold is exactly
the adjustment the pass criteria prohibit, so it was not done.

## 2. The failure is the tissue, not the machinery — evidence

R11 defined this diagnostic in advance. The evidence that the retrieval, build-alignment and
harmonization code is working correctly is that the **outcome** side recovers the known signals
massively, through the identical code path:

| pair | GWAS variants p < 5e-8 in window | strongest GWAS variant | β (harmonized) | SE | p |
|---|---|---|---|---|---|
| HMGCR → LDL-C | 367 / 529 | `chr5_75266739_T_C` | +0.0601913 | 0.00144448 | underflow (< 5e-324) |
| PCSK9 → LDL-C | 446 / 763 | `chr1_55021673_C_G` | −0.0659791 | 0.00164207 | underflow |
| IL6R → CRP | 416 / 546 | `chr1_154423363_C_T` | −0.0814838 | 0.00204503 | underflow |

So: the LDL-C GWAS has hundreds of genome-wide significant variants at both HMGCR and PCSK9, and
the CRP GWAS likewise at IL6R, all retrieved by range query and correctly allele-aligned. The
exposure side — GTEx v10 liver, `sample_size` 269, effective 2·208 = 416 haplotypes as read from
the `an` field — has **zero** genome-wide significant cis-eQTLs for any of the three genes.

This is a power failure of the tissue fixed by pre-registration rule R3, not a broken chain.
**It does not convert any control into a pass, and under the pass criteria no tissue swap is
permitted.**

## 3. Colocalization (computed; does not constitute a pass)

Coloc does not require instruments, so it was computed to exercise that stage. Priors
p1 = p2 = 1e-4, p12 = 1e-5, W = 0.15², trait 1 = eQTL, trait 2 = GWAS.

| pair | n variants | PP0 | PP1 | PP2 | PP3 | PP4 |
|---|---|---|---|---|---|---|
| HMGCR → LDL-C | 529 | 0.000 | 0.000 | **0.9002** | 0.0385 | 0.0613 |
| PCSK9 → LDL-C | 763 | 0.000 | 0.000 | 0.4546 | **0.5033** | 0.0421 |
| IL6R → CRP | 546 | 0.000 | 0.000 | **0.9018** | 0.0361 | 0.0621 |

PP2 is "association with trait 2 (GWAS) only". Its dominance for HMGCR and IL6R is the same
finding restated: strong outcome signal, no exposure signal. PCSK9's PP3 = 0.50 reflects its
partial eQTL signal being distinct from the GWAS peak. None of this is evidence of a positive
control recovering.

## 4. Harmonization — counts, and what the data did *not* exercise

| | HMGCR→LDL | PCSK9→LDL | IL6R→CRP |
|---|---|---|---|
| eQTL variants in window | 592 | 846 | 839 |
| GWAS variants in window | 3,323 | 3,925 | 581 |
| position-matched | 578 | 828 | 546 |
| indel / multi-allelic dropped | 44 | 59 | 0 |
| `aligned_same` | 529 | 763 | 546 |
| `aligned_flipped` | 0 | 0 | 0 |
| `strand_flipped` | 0 | 0 | 0 |
| `allele_mismatch` | 0 | 0 | 0 |
| strand-ambiguous total | 85 | 106 | 0 |
| ↳ resolved by frequency | 80 | 100 | 0 |
| ↳ unresolved → dropped | 5 | 6 | 0 |
| **overlapping variants after harmonization** | **529** | **763** | **546** |

**Caveat you should know about.** Every single merge resolved as `aligned_same`. That is expected —
GWAS Catalog harmonised files and eQTL Catalogue both express alleles relative to the GRCh38
reference, so `effect_allele` equals `alt` throughout. The consequence is that the
**allele-flip and strand-flip branches of the harmonization code were never executed by these
data and remain UNVERIFIED.** They will need a deliberate test (e.g. a synthetic flipped record)
before the pipeline is trusted at scale.

The strand-ambiguous handling *was* exercised: 191 palindromic variants across the three pairs,
180 resolved by allele frequency and 11 dropped as unresolvable.

## 5. Tissue selection rule — stated in the pre-registration, before results

**Rule R3:** GTEx_v10 liver, eQTL Catalogue R8 dataset `QTD000266`, `quant_method = ge`, fixed for
all three pairs; no per-pair selection, no maximum over tissues. Chosen because all three traits
are hepatically determined (LDL-C by hepatic synthesis and LDL-receptor recycling; CRP by
hepatocyte acute-phase synthesis downstream of IL6R). It is also the unique liver `ge` dataset in
R8, so the rule has no discretion:

```bash
awk -F'\t' 'NR>1 && $9=="ge" && tolower($6) ~ /liver/ {print $1,$2,$3,$6,$8}' \
  data/raw/dataset_metadata_r8.tsv     # QTS000015 QTD000266 GTEx_v10 liver 269
```

## 6. Sources, accessions, versions, timestamps

| role | source | accession / version | retrieved (UTC) |
|---|---|---|---|
| outcome GWAS (LDL-C) | GWAS Catalog harmonised | `GCST90239658`, Graham SE 2021, EUR N=1,320,016 | 22:57:27Z / 22:58:35Z |
| outcome GWAS (CRP) | GWAS Catalog harmonised | `GCST90179146`, Koskeridis F 2022, N=575,531 | 2026-08-26T22:59:25Z |
| GWAS study table | GWAS Catalog release | `v1.0.3.1` | 2026-08-26T02:49:13Z |
| cis-eQTL | eQTL Catalogue **R8** | `QTD000266` GTEx_v10 liver `ge`, n=269, PMID 32913098 | 22:57:31Z / 22:58:39Z / 22:59:30Z |
| eQTL permutations | eQTL Catalogue R8 | `QTD000266.permuted.tsv.gz` sha256 `9ac3869a…` | 2026-08-26T22:47:43Z |
| LD panel | 1000 Genomes high-coverage phased | `20220422_3202_phased_SNV_INDEL_SV`, GRCh38 | 22:58:24Z / 22:58:50Z / 23:00:05Z |
| LD sample panel | 1000 Genomes | `20130606_g1k_3202_samples_ped_population.txt` sha256 `4e164b71…` | 2026-08-26T22:40:04Z |
| gene coordinates | Ensembl REST | GRCh38, per-gene, recorded in provenance JSON | 22:57:22Z / 22:58:29Z / 22:59:20Z |

Exact per-stage timestamps are in `deliverables/d1/provenance/*_provenance.json`.

## 7. Retrieval cost actually measured

525 EUR unrelated samples used for LD (`Superpopulation == EUR`, `FatherID == MotherID == 0`).

| pair | window | GWAS bytes (% of file) | eQTL bytes (%) | LD bytes (%) | total | LD wall-clock |
|---|---|---|---|---|---|---|
| HMGCR→LDL | chr5:75,236,329-75,464,241 | 2,028,277 (0.1385 %) | 2,356,944 (0.0811 %) | 2,685,130 (0.1276 %) | 7,070,351 | 52.85 s |
| PCSK9→LDL | chr1:54,939,445-55,164,852 | 2,075,062 (0.1418 %) | 2,519,378 (0.0867 %) | 3,052,260 (0.1167 %) | 7,646,700 | 10.43 s |
| IL6R→CRP | chr1:154,305,141-154,569,450 | 1,673,467 (0.4834 %) | 3,671,276 (0.1263 %) | 3,486,224 (0.1349 %) | 8,830,967 | 35.35 s |

Build alignment verified per pair, not assumed: eQTL `variant` ID decomposed to
`chr_pos_ref_alt` matched the row's own position/ref/alt in **200/200** sampled rows for every
pair; GWAS files are the `.h.` harmonised GRCh38 product; Ensembl reported `GRCh38`; LD panel is
the GRCh38 high-coverage release.

## 8. Files written — all unmodified, open them directly

```
deliverables/d1/
  PREREGISTRATION.md                 rules + amendment R2a, written before results
  egene_status.json                  eGene test, all three genes
  gwas_selection_EFO_0004611.csv     full ranked shortlist, LDL  (292 studies)
  gwas_selection_EFO_0004458.csv     full ranked shortlist, CRP
  d1_summary.json                    every number in this report
  extracts/
    HMGCR_LDL_GWAS_GCST90239658_5_75236329_75464241.tsv
    HMGCR_LDL_EQTL_QTD000266_ENSG00000113161_5_75236329_75464241.tsv
    PCSK9_LDL_GWAS_GCST90239658_1_54939445_55164852.tsv
    PCSK9_LDL_EQTL_QTD000266_ENSG00000169174_1_54939445_55164852.tsv
    IL6R_CRP_GWAS_GCST90179146_1_154305141_154569450.tsv
    IL6R_CRP_EQTL_QTD000266_ENSG00000160712_1_154305141_154569450.tsv
  harmonized/
    HMGCR_LDL_harmonized_merge.tsv   529 rows, 28 cols
    PCSK9_LDL_harmonized_merge.tsv   763 rows
    IL6R_CRP_harmonized_merge.tsv    546 rows
  provenance/
    {HMGCR_LDL,PCSK9_LDL,IL6R_CRP}_provenance.json
```

The harmonized merge tables carry, per variant: `eqtl_ref`, `eqtl_alt`, `gwas_effect_allele`,
`gwas_other_allele`, `harmonization`, `sign_applied`, `ambiguous`, `ambiguous_note`,
`beta_eqtl`, `se_eqtl`, `p_eqtl`, `beta_gwas_raw`, `beta_gwas_harmonized`, `se_gwas`, `p_gwas`,
`af_alt_eqtl`, `af_alt_gwas`, `F_eqtl` — both sides side by side, as requested.

## 9. What this blocks, and the decision that is yours

The chain is verified end-to-end for: region retrieval by range query (three independent sources),
build alignment, allele harmonization with explicit palindromic handling, LD panel retrieval and
r² computation, clumping, MR estimators and coloc. What is **not** verified is that any of it
produces a correct causal estimate, because no control ever reached the MR stage.

The binding constraint is that **GTEx v10 liver is too small (n=269) to instrument these genes**.
That is a study-design question, not a bug, and resolving it means changing a pre-registered rule —
which is your call, not mine. I have not made it.

Non-exhaustive options, offered without a recommendation being acted on:

1. Change R3 to a larger single tissue fixed in advance (e.g. a whole-blood dataset with n in the
   thousands), re-run the same controls. Still one tissue, still no maximum over tissues.
2. Change R3 to a non-GTEx resource with far larger n for these genes.
3. Relax R5's instrument threshold to a pre-registered cis-appropriate value.
4. Accept that cis-eQTL is the wrong exposure for these particular controls and choose controls
   whose canonical instrument is a coding/protein-level variant.

Each of these is a change to the pre-registration and must be recorded in `DECISIONS.md` as such,
with the old rule and the reason, so the audit trail shows the rule changed **before** the
re-run and not because a particular answer was wanted.

**Deliverable 2 has not been started.**
