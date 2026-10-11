# Figures

One script makes all four. Each is written as PDF and SVG (vector) and PNG (300 dpi).

```powershell
cd D:\causal2drug
python manuscript_prep\analysis\classify.py          # inputs
python manuscript_prep\figures\make_figures.py       # all; or add fig1 | fig2 | fig3 | figS1
```

Colours are three slots checked for colour-vision deficiency (blue, orange, aqua; worst pair
ΔE 9.2 under simulated deuteranopia). Aqua has low contrast against white, so every aqua mark
carries a printed value or a distinct marker shape. Status is never shown by colour alone.

## Figure 1 — `figure1_flow`

**Derivation of the cardiovascular target set and instrument status.** Rows of the Minikel et
al. target–indication table were reduced to unique pairs; four programmes with conflicting
phase labels were removed; pairs were kept if the indication is tagged cardiovascular or
matches a lipid term. Each of the 726 targets was classified by whether eQTLGen whole blood
provides a cis-eQTL instrument (p < 5 × 10⁻⁸, gene body ± 100 kb). "No eQTL data": the gene is
not among the 19,250 genes eQTLGen tested. "No qualifying instrument": tested, no variant meets
the criterion. The no-eQTL-data box gives its sub-types: autosomal genes not tested, and
X-linked (27) and mitochondrial (2) genes, which eQTLGen did not analyse. The bracket gives the lower and upper bound of
uninstrumentable targets; the line beneath gives the autosomal-only lower bound. Numbers:
`classify.py` [K01], [K04], [K10].

## Figure 2 — `figure2_status_by_phase`

**Instrument status by highest clinical phase.** Share of targets in each status, by the
highest phase reached in a cardiovascular indication. Solid orange: no eQTL data, autosomal
gene not tested. Hatched: no eQTL data, X-linked or mitochondrial gene not analysed by eQTLGen. The value
printed on the orange segment is the total of both sub-types. The black tick marks the upper
bound (no eQTL data plus no qualifying instrument). Bottom two rows: launched targets with an
approved cardiovascular drug in ChEMBL 37 under the strict and the permissive rule; **both are
provisional** pending manual review of the drug links. Numbers: [K04], [K05], [K06]; Table 2.

## Supplementary Figure S2 — `figureS2_expression_launched_phase3` (was Figure 3 until 2026-10-11)

Figure 3 is now the recovery map: `manuscript/display/figure3_recovery_map`, legend in `manuscript/manuscript.md`.

**Launched (a) and Phase III (b) targets without a blood instrument.** One row per target.
Filled marker: median expression in whole blood (circle, no eQTL data; diamond, no qualifying
instrument). Open circle: the GTEx tissue with the highest median expression, named beside it.
Dashed line: 1 TPM. Values at or below 0.01 TPM are drawn at the left edge. (X) X-linked gene,
not analysed by eQTLGen. Rows are ordered by status, then by highest tissue expression. Tissue
names are GTEx v10 labels. Data: Supplementary Tables S1 and S7.

**Only example-eligible targets are drawn:** an approved drug linked directly to the gene has
an approved, strictly cardiovascular indication in ChEMBL 37 (strict rule), and the target holds
the same phase in sensitivity set 2. Panel a: 13 launched targets (10 with no eQTL data, 3 with
no qualifying instrument). Panel b: 4 targets labelled Phase III in the source data for which
ChEMBL lists such an approved drug. GHR, GLA, PDE4C, the adrenoceptor and endothelin subtypes,
and Phase III targets without an approved drug (LPA among them) are not shown. The selection is
provisional: it follows the strict rule and will change if the manual review changes it.

## Supplementary Figure S1 — `figureS1_classC_negative_result`

**Pre-registered test of the Class C hypothesis: refuted.** Posterior probability of a shared
causal variant (coloc.abf PP4) at 50 loci, computed with the same outcome GWAS, window and
harmonisation and two exposure resources. Points above the diagonal have higher PP4 with the
larger resource (32 loci); below, lower (18). Wilcoxon signed-rank on log10 PP4, p = 0.0262.
The hypothesis predicted lower PP4 with the larger resource. No locus exceeds PP4 = 0.15 in
either arm; the shaded band marks PP4 > 0.8. *Colocalization is not used as evidence for any
target in this study; this figure reports only the pre-registered test of the Class C
hypothesis, which was refuted.* Loci are not labelled, by the terms of
`DEVIATION_CLASS_C_FIGURE.md`. Numbers: `verify_numbers.py` [V09].

## Not carried forward

The four figures in `deliverables/figures/` describe the original 630-target set with the old
class names. Their "required modality" and colocalization content is not in the paper.
