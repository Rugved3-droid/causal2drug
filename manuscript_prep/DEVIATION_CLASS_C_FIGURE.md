# Deviation note — reporting the Class C negative result

**Dated 2026-10-10. Recorded by decision of the project owner (Rugved Parmar).**
This note lives in `manuscript_prep/`. `PROTOCOL.md` and `DECISIONS.md` are not edited.

## The rule being deviated from

`PROTOCOL.md:135`, standing rule under "Gate D closed", 2026-08-31:

> No colocalization output from this project may be used, cited, tabulated or plotted. MR
> verdicts do not depend on coloc and are unaffected.

`DECISIONS.md:259-262` states the same prohibition in shorter form.

## The deviation

The pre-registered Class C refutation test is a paired comparison of coloc.abf PP4 values
(`deliverables/d1h/part3_headtohead.tsv`). Reporting it means tabulating and plotting
colocalization output, which the rule's wording forbids.

**Reporting the pre-registered Class C negative result is an exception to the rule.**

## Reason

The prohibition was written to stop colocalization being used **as evidence about a target**:
no PP4 value may support, weaken or rank any target–indication pair, and none enters the atlas.
That purpose is unchanged.

The Class C result is a different use. It is a pre-registered test about the *method*, with a
refutation criterion written before the data were computed (`PROTOCOL.md:62-72`). It makes no
claim about any target. Leaving it out would hide a hypothesis the project proposed, tested and
lost.

## Scope of the exception

| permitted | not permitted |
|---|---|
| One **supplementary** figure of the paired PP4 comparison (50 loci, GTEx whole blood vs eQTLGen) | Any main-text figure of colocalization output |
| The summary statistics of that test: 32 / 18, medians, Wilcoxon and binomial p | Any per-target PP4 value presented as evidence for or against a target |
| A plain statement that no locus exceeds PP4 = 0.15 in either arm | Gene labels on the figure that invite reading a locus as "colocalized" or "not" |
| | Any SuSiE / credible-set PP4, including the four withdrawn artefact values |

The figure will carry this sentence in its legend: *"Colocalization is not used as evidence for
any target in this study; this figure reports only the pre-registered test of the Class C
hypothesis, which was refuted."*

## What this note does not change

* Gate D stays closed. Colocalization is not an evidence gate.
* The existing Table 2 ("method availability") also tabulates colocalization counts. It is
  **not** covered by this exception. Its headline figure is being corrected or removed under
  `COLOC_CLAIM_CORRECTION.md`, and whether any of it remains is a separate decision recorded
  there.

## Numbers the exception covers

Reproduce with `python manuscript_prep\verify\verify_numbers.py`, block [V09]: 60 loci
attempted, 50 paired; PP4 higher in eQTLGen at 32, lower at 18; median PP4 0.00821 → 0.01390;
Wilcoxon signed-rank on log10 PP4, W = 408.0, p = 0.0262; exact binomial p = 0.0649; maximum
PP4 0.091 (eQTLGen) and 0.142 (GTEx).
