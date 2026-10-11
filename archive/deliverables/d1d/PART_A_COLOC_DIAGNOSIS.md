# Part A — colocalization defect: CAUSE FOUND

Run 2026-08-30. Locus: HMGCR – LDL-C, the same locus that produced PP3 ≈ 1.000 / PP4 ≈ 1e-11.

```bash
cd /d/causal2drug
python src/test_coloc.py        # implementation validation, 6 synthetic cases
python src/coloc_diagnose.py    # eQTLGen-converted vs GTEx-blood-native vs GTEx-blood-SMR
python src/coloc_signals.py     # independent signals + prior sensitivity
```

## Verdict

**The implementation is correct. The Z-score conversion is not the cause. The cause is
coloc-ABF's single-variant assumption applied to two very-well-powered traits whose observed
lead variants differ inside one LD block.**

Gate D remains **NOT TRUSTWORTHY** until an LD-aware colocalization is implemented. The cause is
identified, not yet fixed.

---

## Step 1 — how beta and SE are currently derived (exact code)

```python
def z_to_beta_se(z, af, n):          # src/control_chain_v2.py
    v = 2.0 * af * (1.0 - af) * (n + z * z)
    se = 1.0 / math.sqrt(v)
    return z * se, se
```

Formula (Zhu et al. 2016, SMR):  β = z / √(2p(1−p)(n+z²)),  SE = 1 / √(2p(1−p)(n+z²)).

* **N**: **per-variant**, eQTLGen's own `NrSamples` column — not a global constant.
* **MAF**: eQTLGen's **own** allele-frequency file
  (`2018-07-18_SNP_AF_...txt.gz`, 240,045,342 B, sha256 `edffd43094598d08…`), column
  `AlleleB_all`, explicitly re-oriented to the `AssessedAllele` (if `AssessedAllele == AlleleA`
  the frequency is taken as `1 − AlleleB_all`). All 503 required rsIDs were found; none was
  imputed or defaulted.
* Consequence: β/SE = z **exactly**, so F = z² and instrument strength is independent of MAF.

## Step 2 — decisive test: native beta/SE vs reconstructed

Same locus, same outcome GWAS, three colocalizations. INTERVAL (n=4,729) was the intended
native-beta/SE comparator but is **not yet released on the FTP** (`QTS000053` → HTTP 404, it is in
the "upcoming" table only), so the largest *available* whole-blood `ge` dataset was used: GTEx v10
whole blood, `QTD000356`, n = 853.

| colocalization | n | PP2 | PP3 | PP4 | s1+s2−s12 |
|---|---|---|---|---|---|
| (1) eQTLGen, beta/SE **reconstructed** from Z | 276 | 5.77e-44 | 1.0000 | 1.0e-11 | 32.23 |
| (2) GTEx blood, **NATIVE** beta/SE | 644 | 0.918 | 0.0550 | **0.0271** | 7.62 |
| (3) GTEx blood, beta/SE **re-derived by SMR** | 644 | 0.899 | 0.0605 | **0.0403** | 7.31 |

**(2) vs (3) is the stronger test**: same dataset, same locus, same GWAS, differing *only* by
whether beta/SE are native or reconstructed. PP4 = 0.0271 vs 0.0403 — materially identical.

Scale check: `beta_smr/beta_native` median **1.6012**, `se_smr/se_native` median **1.6012** — the
*same* constant. A common factor on β and SE leaves z unchanged, which is why the ABFs barely move.
The 1.60 factor is a units difference (GTEx expression scaling vs the SMR standardisation), not an
error.

**→ The Z-score conversion is NOT the cause. Primary hypothesis refuted.**

## Step 3 — implementation validation and the single-causal-variant assumption

**Implementation validated against known truth, 6/6** (`src/test_coloc.py`):

| case | want | result |
|---|---|---|
| A. shared causal variant, both strong | PP4 | PASS (PP4 = 1) |
| B. distinct causal variants, both strong | PP3 | PASS (PP3 = 1) |
| C. outcome only | PP2 | PASS (0.992) |
| D. exposure only | PP1 | PASS (0.991) |
| E. neither | PP0 | PASS (0.997) |
| F. **shared causal, GWAS 1.3M vs eQTL 32k** | PP4 | **PASS (PP4 = 1)** |

Case F is the test of the N-imbalance hypothesis (step 4): with the *same* precision imbalance as
the real data and a shared causal variant, PP4 = 1. **N imbalance does not break H4. Hypothesis
refuted.**

**Independent signals — the multiple-causal-variant hypothesis is also refuted.** LD-clumping each
trait separately inside the window (1000G GRCh38, EUR unrelated n = 525, r² < 0.01):

* variants with p < 5e-8: eQTL **276**, GWAS **203**
* **independent signals: eQTL 1 (rs6453133), GWAS 1 (rs12916)**

There is exactly **one signal per trait**. coloc-SuSiE would not help *for the reason of multiple
signals*, because there are none.

## Step 4 — the actual cause, closed arithmetically

| | |
|---|---|
| eQTL lead | **rs6453133**, GRCh38 75,396,951, z = 14.97, p = 1.21e-50 |
| GWAS lead | **rs12916**, GRCh38 75,360,714, z = 49.59, p underflowed to 0 |
| distance between leads | **36,237 bp** |
| **LD between leads (1000G EUR)** | **r = 0.681, r² = 0.464 — same LD block** |
| rank of rs12916 on the eQTL side | **3 of 276** |
| rank of rs6453133 on the GWAS side | 87 of 276 |

Because lABF ≈ ½z² when the effect is well determined, the eQTL log-ABF difference between the two
leads is

```
lABF_eQTL(rs6453133) = 108.77   (z = 14.97)
lABF_eQTL(rs12916)   =  76.54   (z = 12.63)
difference           =  32.23   ->  e^32.23 = 1.0e14
observed s1+s2-s12   =  32.23   -- identical
```

and PP3/PP4 = e^32.23 × (p1·p2/p12) = 1e14 × 1e-3 = 1e11, giving PP4 = 1.0e-11. The observation is
reproduced to the decimal by this single quantity.

**Controlled demonstration.** Holding a single shared causal signal fixed by construction and
displacing only the *observed peak* of trait 1:

| peak displacement | PP3 | PP4 | s1+s2−s12 |
|---|---|---|---|
| 0 | 6.2e-09 | **1.000** | 0.00 |
| 1 | 0.9969 | 3.1e-03 | 12.68 |
| 2 | 1.000 | 4.1e-08 | 23.93 |
| **3** | **1.000** | **1.9e-12** | **33.90** |
| 5 | 1.000 | 1.1e-19 | 50.60 |

A displacement of three variants — with one shared causal signal throughout — reproduces the HMGCR
observation almost exactly. This is a property of coloc-ABF, not of this code.

**Prior sensitivity (diagnostic only — the pre-registered p12 = 1e-5 stands, nothing was tuned):**

| p12 | 1e-6 | 1e-5 (pre-registered) | 1e-4 | 1e-3 | 1e-2 |
|---|---|---|---|---|---|
| PP4 | 1.0e-12 | 1.0e-11 | 1.0e-10 | 1.0e-09 | 1.0e-08 |

PP4 stays negligible across four orders of magnitude. **No choice of prior rescues it** — which is
exactly why tuning the prior would have been both futile and dishonest.

## What this means, and why it is not only a coding matter

coloc-ABF compares traits **one variant at a time**. Its discriminating power scales as ½z², so the
more powerful the datasets, the more brittle it becomes to which variant happens to peak. With
eQTLGen at n = 31,684, eQTL z-scores of ~15 are routine; a lead mismatch of a few variants inside a
single LD block — caused by differential power and sampling noise, **not** by different causal
variants — is then enough to force PP3 → 1.

This will affect **every** locus in this project, not only HMGCR, and it gets worse as the exposure
resource gets larger. It is a systematic property of applying coloc-ABF to very large eQTL
resources.

**Indicated remedy:** LD-aware / credible-set colocalization (coloc-SuSiE or equivalent), which
compares *credible sets* using the LD matrix instead of comparing single variants, and is therefore
insensitive to which variant inside a block happens to peak. Note this is needed even though there
is only one signal per trait — the reason is lead-variant mismatch, not multiplicity.

**Not implemented in this pass.** Gate D therefore stays **NOT TRUSTWORTHY** in `PROTOCOL.md`. The
cause is found and closed arithmetically; the fix is new work and is not being smuggled in
untested.

## Incidental finding

**53 of 276 GWAS p-values in this window underflowed to exactly 0.0** in the harmonised file.
Colocalization and MR both use β and SE rather than p, so nothing here is affected — but any future
step that reads `p_value` numerically must not treat 0.0 as a valid probability. Recorded so it is
not rediscovered as a bug later.
