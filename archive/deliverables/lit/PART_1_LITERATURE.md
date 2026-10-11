# Part 1 — literature check: what is already published

Searched 2026-08-30. **Bottom line: both mechanisms are already described in the literature.
Our contribution is not the mechanism. What is NOT published is the footprint — nobody has
quantified how much of the eQTL–GWAS colocalization literature is affected.**

---

## 1. Your recollection: CONFIRMED on both counts, with one correction

### Kanai / SLALOM — confirmed, and it is more directly on point than you may realise

**Kanai M, et al. "Meta-analysis fine-mapping is often miscalibrated at single-variant
resolution." *Cell Genomics* 2022; 2(12):100210.** (medRxiv 2022.03.16.22272457; PMC9839193)

* SLALOM = **s**uspicious **l**oci **a**nalysis of **m**eta-analysis summary statistics.
* Detects association-statistic outliers against local LD structure — i.e. the *same*
  observed-vs-LD-predicted residual logic I implemented independently as the residual-z RMSE
  check.
* Explicitly frames LD reference mismatch as **"all but guaranteed in the case of GWAS
  meta-analyses."**
* Headline number: applied to 14 Global Biobank Meta-analysis Initiative meta-analyses,
  **67 % of loci show suspicious patterns that call fine-mapping accuracy into question.**
* Predicted-suspicious loci are depleted 2.7× for nonsynonymous lead variants
  (Fisher p = 7.3 × 10⁻⁴) — an external validity check that the flag is real.

**Correction to your framing:** SLALOM is about *meta-analysis* fine-mapping miscalibration
driven by cross-cohort heterogeneity (arrays, imputation panels, sample size), not solely
about the reference-panel-too-small problem. That distinction matters for us, because our
GWAS (Graham 2021 GLGC, n = 1,320,016) **is** a meta-analysis of heterogeneous cohorts, so
SLALOM's mechanism applies to our locus on its own terms — independently of the n = 525
panel size.

### DENTIST — confirmed, and it is the published version of my diagnostic

**Chen W, Wu Y, Zheng Z, Qi T, Visscher PM, Zhu Z, Yang J. "Improved analyses of GWAS summary
statistics by reducing data heterogeneity and errors." *Nature Communications* 2021;
12:7117.** (PMC8654883)

* DENTIST = **D**etecting **E**rrors i**N** analyses of summary s**t**at**ist**ics.
* Method: compares each variant's observed GWAS test statistic against the value **predicted
  from neighbouring variants using LD from a reference panel**, and removes discordant
  variants. That is precisely the `z_j ≈ z_lead · r(lead,j)` residual test I built.
* Detects genotyping/imputation error, allele mislabelling, **and heterogeneity between the
  GWAS and the LD reference sample** — the third being our case.
* Reports false-positive rate for secondary signals in conditional analysis falling from
  **> 28 % to < 2 %** in the presence of GWAS-vs-LD-reference heterogeneity.

**So my "LD-consistency diagnostic" is an unwitting reimplementation of DENTIST.** I should
have cited it rather than presented the residual-z RMSE as a bespoke check. Recording that
plainly.

---

## 2. Is the ABF peakiness / lead-displacement mechanism already described? **YES**

**Hukku A, Pividori M, Luca F, Pique-Regi R, Im HK, Wen X. "Probabilistic colocalization of
genetic variants from complex and molecular traits: promise and limitations." *Am J Hum
Genet* 2021; 108(1):25–35.** (PMC7820626)

This paper describes the mechanism directly:

> "many mismatches between causal and lead SNPs lead to **false negatives** in identifying the
> colocalized signals"

and classifies it as **"class II false negatives"** — both association signals are genuinely
present, but SNP-level quantification fails to detect the colocalization. That is exactly the
HMGCR pattern we measured (PP3 ≈ 1, PP4 = 1e-11, one signal per trait, leads 36 kb apart at
r² = 0.464).

**What is NOT in that paper:** the explicit algebra I derived — that lABF ≈ ½z² so the lBF gap
between two variants scales as ½(z_j² − z_k²), making the failure *worse as z grows*, i.e. as
the exposure resource becomes more powerful. Hukku et al. describe the failure; they do not
derive its **power dependence**. I could not find that derivation stated anywhere.

Related and relevant:

* **Wallace C. "Eliciting priors and relaxing the single causal variant assumption in
  colocalisation analyses." *PLOS Genetics* 2020; 16(4):e1008720.** (PMC7192519)
* **Wallace C. "A more accurate method for colocalisation analysis allowing for multiple
  causal variants." *PLOS Genetics* 2021; 17(9):e1009440.** (PMC8504726) — the coloc+SuSiE
  paper. Notably it frames SuSiE as the fix for **multiple causal variants**, consistent with
  my own correction that it is *not* the fix for lead displacement with a single signal.
* The coloc documentation itself notes that the single-variant route "avoids the need to
  carefully match reference LD matrices," i.e. the LD-matching burden is a known cost of the
  SuSiE route — which is exactly the wall we hit.

---

## 3. Has anyone QUANTIFIED the affected fraction of the literature? **NO — this is the gap**

* **Hukku et al. 2021** give the closest thing: for well-powered trait–tissue pairs the
  average rejection-to-expected ratio is **10.7 %**, implying ~89 % of expected
  colocalizations go undetected. But that is aggregate power loss, **not** attributed to our
  two mechanisms specifically, and not a survey of published results.
* **Kanai et al. 2022** quantify **67 % of GBMI loci as suspicious** — but for *fine-mapping*
  in meta-analyses, not for *colocalization*, and not for eQTL–GWAS pairs.
* Searches for a quantification of how many published eQTL–GWAS colocalizations are affected
  by LD-reference mismatch or by lead-displacement returned **nothing**. The recurring
  statement in the literature is that there is "a lack of empirical evaluation of how
  different levels of LD mismatch affect colocalization outcomes."

**One paper runs directly against our Part 3 claim and must be confronted, not ignored:**

**Rosen JD, Broadaway KA, Brotman SM, Mohlke KL, Love MI. "Higher eQTL power reveals signals
that boost GWAS colocalization." bioRxiv 2025.08.05.668745.**

* Central claim: **higher eQTL power *increases* colocalization detection.** Detection of
  eQTLs rises from 15.6 % at n = 500 to 36.8 % at n = 2,000; colocalizing GWAS signals
  roughly **double** going from METSIM-S (n = 420) to AdipoExpress (n = 2,256).
* They report **no** instance of PP4 decreasing with higher power, and do **not** discuss
  lead-variant mismatch. They used single-variant coloc with PP4 ≥ 0.5.

**This is not necessarily a contradiction, and I will not pretend it is settled.** Their
comparison spans n = 420 → 2,256; ours would span n = 853 → 31,684. At low power the binding
constraint is detecting an eQTL *at all*, so more power plainly helps. Our proposed mechanism
only bites once z is large enough for the ½z² gap to dominate (z ≈ 15–50, reached around
n ≈ 30,000). The two claims could both be true with a **non-monotone** relationship — PP4
rising with power, then falling. That is a sharper and more interesting hypothesis than the
one we started with, **and it is also a live route to a null**: the head-to-head in Part 3 may
simply show PP4 higher in eQTLGen, which would refute our claim outright.

---

## 4. What this means for the project — plainly

**Established, not ours:**
1. LD reference mismatch miscalibrates summary-statistic fine-mapping (Kanai 2022; Chen 2021).
2. The observed-vs-LD-predicted residual is the standard diagnostic for it (DENTIST).
3. Lead/causal-variant mismatch causes colocalization false negatives (Hukku 2021, "class II").
4. coloc+SuSiE addresses multiple causal variants, at the cost of requiring well-matched LD
   (Wallace 2021).

**Possibly novel, stated cautiously:**
5. The explicit **power dependence** — that the lBF gap scales as ½(z_j² − z_k²), so the
   failure *worsens monotonically as the exposure resource improves*. Described nowhere I
   could find. **UNVERIFIED as novel** — absence of evidence from four searches is weak
   evidence of absence, and this is exactly the kind of result that hides in a methods
   appendix.
6. The **footprint**: nobody has quantified what fraction of eQTL–GWAS colocalization results
   are affected by either mechanism. This is a genuine, clearly-stated gap.

**The honest recommendation:** our contribution should be framed as the *footprint*, exactly
as you anticipated — and the mechanism sections should cite Hukku, Kanai and Chen rather than
present the mechanisms as discoveries. The one piece worth claiming as new is the power
dependence, and Part 3's head-to-head is what would earn it. Rosen et al. means Part 3 has a
real chance of returning a null, and that null would be worth reporting.

---

## Sources

- [Kanai et al. 2022, Cell Genomics — Meta-analysis fine-mapping is often miscalibrated at single-variant resolution](https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00163-X) · [PMC9839193](https://pmc.ncbi.nlm.nih.gov/articles/PMC9839193/) · [SLALOM code](https://github.com/mkanai/slalom)
- [Chen et al. 2021, Nature Communications — Improved analyses of GWAS summary statistics by reducing data heterogeneity and errors (DENTIST)](https://www.nature.com/articles/s41467-021-27438-7) · [DENTIST code](https://github.com/Yves-CHEN/DENTIST)
- [Hukku et al. 2021, AJHG — Probabilistic colocalization of genetic variants from complex and molecular traits: promise and limitations](https://www.cell.com/ajhg/fulltext/S0002-9297(20)30409-2) · [PMC7820626](https://pmc.ncbi.nlm.nih.gov/articles/PMC7820626/)
- [Wallace 2021, PLOS Genetics — A more accurate method for colocalisation analysis allowing for multiple causal variants](https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1009440)
- [Wallace 2020, PLOS Genetics — Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7192519/)
- [Reales et al., PLOS Genetics — Design and interpretation of eQTL-GWAS colocalisation studies: Lessons from a large-scale evaluation](https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1012141)
- [Rosen et al. 2025, bioRxiv — Higher eQTL power reveals signals that boost GWAS colocalization](https://pmc.ncbi.nlm.nih.gov/articles/PMC12340864/)
