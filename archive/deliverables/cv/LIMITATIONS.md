# Limitations

Every number below ships with the command that reproduces it. Anything not established by
a command in this repository is marked **UNVERIFIED** and is not used to support a claim.

Run from the repository root (`D:\causal2drug`).

---

## 1. Class B is under-counted at N = 1, and is bounded rather than assigned

```bash
python src/cv_atlas.py           # the atlas, including untestability_class
python src/a1_classB_bound.py    # the bound
```

`untestability_class = CLASS_B` is assigned to **exactly one target: IL6R**. That is not an
estimate of how many Class B targets exist. It is the number we could assign, because no
resource labels "this drug's mechanism operates at a level total mRNA cannot represent".

**The true count is UNVERIFIED.** It is bounded, not measured:

| bound | flagged | of instrumented CV targets | reading |
|---|---|---|---|
| **primary — protein-altering association for a cardiovascular trait** | **62** | **386 (16.1 %)** | the figure to cite |
| loose — protein-altering association for any trait | 196 | 386 (50.8 %) | secondary, discounted |
| primary, restricted to launched-drug targets | 6 | 27 (22.2 %) | |
| loose, restricted to launched-drug targets | 16 | 27 (59.3 %) | |
| loose, Phase III or launched | 36 | 60 (60.0 %) | |

**Why a coding association bounds Class B rather than assigning it.** The rule flags a
target when the GWAS Catalog reports an association whose variant has a protein-altering
consequence in that gene. That shows the gene's *mechanism can* act at the protein level.
It does not show that *the drug's* mechanism does. IL6R is Class B because tocilizumab's
effect derives specifically from rs2228145 shedding — a fact established by pharmacology,
not by the presence of a coding variant. Every flagged target is a candidate; none is a
case. The atlas therefore keeps the assignment at N = 1 and carries the bound separately
in `classB_upper_bound.csv` (196 rows).

The loose bound is reported for completeness and **should be discounted**: a
protein-altering association for an unrelated trait says little about a cardiovascular
drug mechanism. Lead with 16.1 %.

**Consequence for the reader.** Class A is a floor on untestability, not a ceiling. Every
one of the 386 instrumented targets looks normal in every table in this paper, and an
unknown fraction of them are Class B — analyses that run, return confident numbers, and
do not mean what they appear to mean.

---

## 2. The FIX 1 assertion rejects genuine colocalizations

```bash
python src/test_coloc_assertion.py   # the assertion and its cost
python src/test_susie_coloc.py       # synthetic suite, 4/4 with strict=False
```

`ColocArtefact` in `src/susie_coloc.py` is active by default (`strict=True`) and **raises**
rather than returning a number. It rejects a credible-set pair when either:

* **(a)** both credible sets top at the **same variant** — a component colocalizing with
  itself; or
* **(b)** the winning pair's tops match **neither trait's lead**.

**The cost, measured rather than assumed.** Condition (a) **also rejects a true
colocalization**: synthetic case **B1**, whose ground truth is a single shared causal
variant, is rejected at **PP4 = 1**. Condition (a) cannot distinguish a real shared causal
variant from a component colocalizing with itself, because the two produce the same
pattern. This is not a hypothetical cost; it is a demonstrated one, and it is the exact
case colocalization exists to detect.

**Why conservative is the safe direction here — and only here.** Gate D is closed and
colocalization is excluded from the atlas as an evidence gate (§3). Nothing in this paper
depends on a colocalization being *accepted*, so an assertion that refuses true positives
costs this analysis nothing and protects it from the four artefacts in §9.

**Cost to anyone reusing this code.** If you use `susie_coloc.py` to generate evidence
rather than to guard against artefacts, this assertion will silently cost you true
positives at an unmeasured rate. You would need a better discriminator — one that
separates a shared causal variant from an over-fitted component, which condition (a) does
not. The synthetic validation suite calls `strict=False` with a comment saying why, and
passes 4/4 under that setting. **The rate at which the assertion rejects true positives on
real data is UNVERIFIED.**

On live data the assertion is doing work, not sitting idle: under the corrected R15 locus
set it rejected one real locus (**TNF**), and both surviving Part 3 results pass
(`ADORA1` PP4 = 0.9587, `F10` PP4 = 0.0904).

---

## 3. Gate D is closed as a documented negative result

```bash
python src/coloc_diagnose.py    # deliverables/d1d/PART_A_COLOC_DIAGNOSIS.md
```

Colocalization returned **PP3 ≈ 1.000, PP4 ≈ 1 × 10⁻¹¹ at HMGCR–LDL-C**, one of the most
robust colocalizations in human genetics. Near-zero PP4 there indicates a specification
defect, not biology. Two independent mechanisms were identified, **both already published**
— this project's contribution is the footprint, not the mechanism:

**Mechanism 1 — single-variant discrimination (lead/causal-variant mismatch).**
coloc-ABF compares single variants and its discrimination scales as ½z². The eQTL lead
(rs6453133, z = 14.97) and the GWAS lead (rs12916, z = 49.59) sit **36,237 bp apart at
r² = 0.464** — same LD block, one signal each. The eQTL log-ABF difference between them is
**32.23**, which equals the observed s1 + s2 − s12 exactly and yields PP4 = 1 × 10⁻¹¹.
This is the class II false-negative mechanism of **Hukku et al. 2021, AJHG 108(1):25–35**.
It affects every locus in this project and **worsens as the exposure resource grows**.

**Mechanism 2 — LD reference misspecification.** SuSiE-RSS + coloc.susie was implemented
from scratch and validated 4/4 against synthetic ground truth, including the multi-signal
case it is designed for. On the real HMGCR–LDL-C locus it gives **PP4 = 0.148** against a
pre-fixed pass criterion of > 0.8 — a **BLOCKER**. Cause: a **525-sample** 1000G EUR panel
cannot reproduce the LD of a **1,320,016-sample** meta-analysis at z ≈ 50. The residual-z
RMSE of the single-causal-variant relation is **7.80 on the GWAS side** (expected ≈ 1)
versus **1.33 on the eQTL side**, so SuSiE-RSS attributes the mismatch to spurious
components and no credible set tops the true GWAS lead. This is the miscalibration
described by **Kanai et al. 2022, Cell Genomics 2(12):100210** (SLALOM), which reports LD
reference mismatch as "all but guaranteed" in GWAS meta-analyses. The residual-z diagnostic
is **DENTIST (Chen et al. 2021, Nat Commun 12:7117)**, independently reimplemented here and
cited as theirs wherever it appears. The credible-set method is **Wallace 2021, PLOS Genet
17(9):e1009440**.

**Gate D therefore has no trustworthy implementation.** Single-variant coloc is unusable
(peakiness); credible-set coloc is unusable with this LD panel (misspecification). What
would be needed is an LD reference matched to the outcome GWAS. Whether one is obtainable
for GWAS Catalog studies is **UNVERIFIED and was not investigated** — that is the "switch
resources" move the pre-registered criterion forbids without an explicit decision.

**Consequence.** Colocalization appears in this paper only as *method availability*
(Table 2), never as evidence. No colocalization output supports any claim in the atlas.

---

## 4. eQTLGen is a single whole-blood exposure resource

```bash
python src/d1b_census.py    # yield census
python src/cv_atlas.py      # CV-restricted yield
```

Every instrument in the atlas comes from **one resource**: eQTLGen phase I cis-eQTL, whole
blood, n = 31,684. Yield is **61.3 % of CV targets** (386/630), against 62.6 % genome-wide.

**What a single whole-blood total-mRNA resource structurally cannot see**, regardless of
sample size:

* **Genes not expressed in blood.** 28.1 % of CV targets (177/630) are absent from the
  file entirely. No increase in n fixes this; the transcript is not there to measure.
* **Splicing.** A variant that shifts isoform usage without changing total transcript
  abundance is invisible. This is one of the two mechanisms that make IL6R Class B.
* **Protein-level regulation.** Shedding, cleavage, secretion and post-translational
  control leave total mRNA unchanged. The other IL6R mechanism.
* **Cell-type-specific effects.** Whole blood is a mixture; an effect confined to one cell
  type is diluted and can be lost.
* **Context-dependence.** Effects present only under stimulation, disease state, or drug
  exposure are absent from a resting-state resource.
* **The tissue where cardiovascular drugs act.** Blood is used because samples are
  abundant, not because it is mechanistically right (Figure 2).

The cardiac, vascular and hepatic columns partially address the last point but are **upper
bounds** (§12), so they cannot be substituted for blood as a primary exposure.

---

## 5. Class C was proposed, pre-registered, tested and refuted

```bash
python src/part3_headtohead.py
```

Class C was our own proposed third untestability class: that colocalization power-dependence
would make some targets untestable. A refutation criterion was **pre-registered before the
measurement was computed** (`PROTOCOL.md`): *if PP4 increases with exposure power in the
paired head-to-head, that refutes the claim, and Class C is dropped — not reframed, not
narrowed, not explained away.*

Paired head-to-head, 50 targets, same outcome GWAS / window / harmonization, only the
exposure resource differing:

* PP4 **higher** in the more powerful resource: **32 / 50**; lower: 18 / 50
* exact binomial **p = 0.0649**; **Wilcoxon signed-rank on log₁₀(PP4) p = 0.0262**
* median PP4: GTEx blood (n = 853) **0.00821** → eQTLGen (n = 31,684) **0.0139**

**PP4 increases with exposure power.** The criterion fired and **Class C is dropped.**
This runs in the same direction as **Rosen et al. 2025 (bioRxiv 2025.08.05.668745)**,
"Higher eQTL power reveals signals that boost GWAS colocalization", which reports the same
effect over n = 420 → 2,256. Classes A and B are unaffected.

The refutation is reported here rather than omitted because a dropped hypothesis that was
pre-registered and tested is a result, and because omitting it would misrepresent how the
two surviving classes were arrived at.

---

## 6. No intermediate power point is obtainable, so non-monotonicity is unverifiable

```bash
# census the R8 'ge' datasets and probe the FTP for the intermediate ones
python - <<'EOF'
import csv, urllib.request
rows = list(csv.DictReader(open("data/raw/dataset_metadata_r8.tsv", encoding="utf-8"),
                           delimiter="	"))
ge = [(int(r["sample_size"]), r["dataset_id"], r["study_id"], r["study_label"])
      for r in rows if r.get("quant_method") == "ge" and r["sample_size"].isdigit()]
for n, d, st, lab in sorted(ge, reverse=True)[:5]:
    u = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
         "%s/%s/%s.all.tsv.gz" % (st, d, d))
    try:
        urllib.request.urlopen(urllib.request.Request(u, method="HEAD"), timeout=45)
        code = "200 available"
    except Exception as e:
        code = str(e)
    print(n, d, st, lab, code)
EOF
```

The Class C measurement rests on **two power points**: GTEx v10 whole blood
(`QTD000356`, **n = 853**) and eQTLGen phase I (**n = 31,684**). A third, intermediate
point was sought.

**Exactly two R8 `ge` datasets fall in n ∈ [2,000, 10,000]**, and **neither is
retrievable**:

| dataset | study | label | n | FTP status |
|---|---|---|---|---|
| QTD001000 | QTS000053 | INTERVAL_RNA | 4,729 | **HTTP 404** — in the "upcoming" table only |
| QTD001003 | QTS000054 | INTERVAL_RNA_WGS | 2,892 | **HTTP 404** |

The next two largest (`QTD000891`, `QTD000881`, OneK1K_reannotated, n = 997 and 991,
QTS000048) are **also 404**. The largest **retrievable** `ge` dataset is `QTD000613`
(OneK1K, QTS000038, **n = 954**), verified as HTTP 200.

The precise claim is therefore *not* that no intermediate dataset exists — two are listed
in the R8 metadata table — but that **none is obtainable**, which was confirmed by probing
the FTP directly rather than inferred from the metadata.

**Two power points give a slope, not a shape.** A non-monotonic relationship — PP4 falling
then rising, or rising then falling — would be invisible to this design. **No curve is
claimed, and non-monotonicity is UNVERIFIABLE with obtainable data.** The refutation in §5
is a directional result between two points and is not weakened by this; a claim about the
functional form would be, and none is made.

The same unavailability shaped the Gate D diagnosis: INTERVAL (n = 4,729) was the intended
native-beta/SE comparator there too, and the largest *available* whole-blood `ge` dataset
(GTEx v10, n = 853) was substituted. That substitution is recorded in
`deliverables/d1d/PART_A_COLOC_DIAGNOSIS.md`.

## 7. The cardiovascular set definition is a decision, and it is shown to be immaterial

```bash
python src/a2_robustness.py
```

Minikel classifies the atherogenic lipid disorders as **metabolic (C18), not cardiovascular
(C14)**. A strict C14 definition would therefore exclude hypercholesterolemia and
hyperlipoproteinemia type II — the PCSK9 and statin indications — and so exclude the atlas's
own anchor example and the entire lipid-lowering class. **CV-EXTENDED** (C14 ∪ a fixed
eight-token match on the MeSH term) is the atlas definition. The rule is mechanical, not
case-by-case.

**Every headline result was recomputed under both definitions and every conclusion
survives:**

| | **CV-PRIMARY** (C14 only) | **CV-EXTENDED** (atlas) |
|---|---|---|
| distinct targets | 594 | 630 |
| CLASS A total | 215 (36.2 %) | 232 (36.8 %) |
| yield: whole blood | 62.0 % | 61.3 % |
| yield: heart LV | 7.4 % | 7.1 % |
| yield: cardiac union | 18.0 % | 17.5 % |
| Class A enrichment OR | 10.60 [4.81, 23.32] | 10.02 [4.77, 21.03] |
| Fisher p | 1.83 × 10⁻¹⁴ | 4.05 × 10⁻¹⁵ |
| CLASS A with launched drug | 17 (corrected 2026-10-10; was 20) | 21 |
| infrastructure attrition | 218/401 = 54.4 % | 265/448 = 59.2 % |

The extension adds 36 targets and shifts no headline figure by more than about one
percentage point. **PCSK9 is in the set under both definitions**, so the anchor does not
depend on the extension. The enrichment is if anything *stronger* under the stricter
definition.

---

## 8. The R14.1 selection rule was flawed; what it produced is withdrawn

```bash
python src/fix2_r15.py
```

**The flaw.** R14.1 selected loci on `n_overlap` and **did not require the outcome GWAS to
be genome-wide significant in the window**. It therefore ran the LD-dependent arms mostly on
loci with no disease signal at all: median GWAS lead |z| = **2.6**, with only **1 of 28**
loci genome-wide significant. Colocalization is only interpretable when *both* traits are
associated; the exposure side is significant by construction, so omitting the outcome-side
criterion admitted loci where there was nothing to colocalize and nothing to misspecify.

**What was withdrawn.** The claim that credible-set colocalization "returns an answer at
only 7 % of loci" (**2 / 28**) — an artefact of selecting loci with no GWAS signal. Also
withdrawn: the disease-side residual-z figures under that rule (median **0.71**, **0 %** of
loci above 2, n = 27), which looked clean only because those loci had no signal to
misspecify.

**What replaced it.** **R15**, declared in `PREREGISTRATION.md` before running and
descriptive rather than confirmatory (Class C was already dead, so no live hypothesis could
be contaminated): `n_overlap ≥ 100`; **the outcome GWAS has ≥ 1 variant with p < 5 × 10⁻⁸
inside the merged window**; one locus per target by largest overlap; **no cap**.

| quantity | R14.1 (withdrawn) | **R15 (corrected)** |
|---|---|---|
| N loci | 28 | **6** |
| genome-wide significant on disease side | 1 / 28 | **6 / 6 (100 %)** |
| DENTIST RMSE, expression side — median | 4.46 | **3.15** |
| DENTIST RMSE, expression side — % > 2 | 93 % | **100 %** |
| DENTIST RMSE, disease side — median | 0.71 (n = 27) | **1.70** |
| DENTIST RMSE, disease side — % > 2 | 0 % | **33 %** |
| credible-set pair returned | 2 / 28 (7 %) | **4 / 6 (67 %)** |

**What survives the correction:** LD misspecification on the high-z expression side is if
anything *more* pervasive — 100 % of R15 loci exceed RMSE 2. **What does not:** the 7 %
figure, withdrawn above.

Note that R15 screened **106** candidate loci and only **6** qualified. That 94 % of
instrumented target loci have no disease signal in the window is itself a finding about
method availability, and is reported as such in Table 2.

One R14.1 locus (**CA2**) carries a non-finite disease-side residual-z and is excluded from
that row only, giving n = 27 there against 28 elsewhere. `float("nan")` parses without
error, so this value silently corrupts a median if it is not screened out; `stat()` in
`src/make_table2.py` excludes non-finite values explicitly for that reason.

---

## 9. Four artefactual PP4 values were produced and withdrawn

```bash
python src/test_coloc_assertion.py
```

| artefact | configuration | PP4 | caught by |
|---|---|---|---|
| HMGCR–LDL-C | 1000G n = 525, s² = 1 | 0.998 | (a) — both top at index 108 |
| **IL6R–CRP** | 1000G n = 525 | **1.000** | (a) — both top at index 93 |
| HMGCR–LDL-C | UK Biobank LD n = 337,000 | 0.9957 | (a) — both top at index 62 (rs1544755, neither lead) |
| HMGCR–LDL-C | 1000G, first draft | 1.000 | (a) — the originally withdrawn value |

**The signature that identified them:** both credible sets top at the **same variant index**
— a SuSiE component colocalizing with itself, which manufactures PP4 ≈ 1 out of an
over-fitted decomposition. In the UK Biobank case the shared top was *neither trait's lead*,
so condition (b) fires too; in the HMGCR/1000G case the shared top *was* the exposure's
lead, so condition (b) alone would have missed it. **Both conditions are required**, which
is why the assertion carries the true-positive cost described in §2.

**IL6R–CRP PP4 = 1.000 is excluded from every table and figure in this paper.** It is the
same over-fitting artefact as the withdrawn HMGCR value, on a locus with no independent
sanity check to expose it. It is recorded here so that its absence is documented rather than
silent.

---

## 10. The liftover failure rate is verified on two loci only

```bash
python src/control_chain_v2.py    # per-pair provenance JSON in deliverables/d1c/
```

GRCh37 → GRCh38 liftover (pyliftover, UCSC `hg19ToHg38` chain) recorded **0 failures across
503 variants in two cis-windows** (277/277 and 226/226). Declared handling, though no case
arose: a variant that fails to lift is **dropped and counted**, never carried forward on its
GRCh37 coordinate; multi-mapping variants take the first chain segment and are counted
separately.

**A 0 % failure rate must not be generalised.** It is expected for common SNPs in gene-dense
autosomal regions. Genome-wide liftover failure is typically 0.1–1 %, concentrated in
repetitive and assembly-corrected regions. **The genome-wide failure rate for this pipeline
is UNVERIFIED** — it was measured only across those 503 variants.

---

## 11. Thirty-two LD loci are excluded as superseded

```bash
python src/fix2_r15.py
```

The R14.1 set selected 60 loci; **28 completed and 32 did not** — the 1000 Genomes LD
endpoint stalled repeatedly, blocking 15–20 minutes on individual loci. The 32 are
**excluded with reason: superseded**, not left pending. They were selected by a rule now
known to be inadequate (§8), so completing them would add loci that mostly fail R15's
criterion anyway.

This is a disclosed exclusion, not a silent one: the final N for the flawed-rule comparison
is **28**, labelled as such, and the final N for the corrected measurement is **6**, complete
with no cap and nothing pending. Checkpoints for the completed loci are in
`deliverables/d1h/ckpt/`.

---

## 12. Atlas-level caveats

```bash
python src/cv_atlas.py
```

* **Cardiac, vascular and hepatic instrumentation figures are UPPER BOUNDS.** They use the
  eQTL Catalogue permuted files' top cis p-value over a **±1 Mb** window — wider than the
  pre-registered ±100 kb — and are **unclumped**. The true yields are lower than the 7.1 %
  (heart LV) and 17.5 % (cardiac union) reported. This makes the blood-versus-tissue gap in
  Figure 2 a **conservative** statement of the problem.
* **12 targets (1.9 %) are UNMAPPED** — no unique GRCh37 gene record — and are excluded from
  the class statistics rather than assigned to a class.
* **`required_modality` is a declared heuristic** over expression and secretion, not a
  curation of drug mechanism. It says where a signal *could* exist, not that one does.
* **The atlas records whether a question can be asked, not its answer.** No column states
  that any target does or does not affect disease.
* **Counts are genes and target–disease pairs, never patients.** Nothing here is a clinical
  decision aid.

---

## Summary of what is UNVERIFIED

| claim | status |
|---|---|
| true Class B count among instrumented CV targets | **UNVERIFIED** — bounded at 16.1 % (primary) / 50.8 % (loose), not measured |
| rate at which the FIX 1 assertion rejects true positives on real data | **UNVERIFIED** — demonstrated on synthetic B1, unmeasured at scale |
| whether an LD reference matched to GWAS Catalog outcome studies is obtainable | **UNVERIFIED** — not investigated, by pre-registered rule |
| non-monotonicity of PP4 with exposure power | **UNVERIFIABLE** — no intermediate power point exists |
| genome-wide liftover failure rate | **UNVERIFIED** — measured on 503 variants in two loci |
| retrieval URLs for GTEx v10, Human Protein Atlas, GWAS Catalog associations, Open Targets 26.06 | **UNVERIFIED** — not captured in `logs/retrieval_log.csv` |
| licence terms of incorporated sources | **UNVERIFIED** — stated as commonly documented, not re-verified against source terms pages |
