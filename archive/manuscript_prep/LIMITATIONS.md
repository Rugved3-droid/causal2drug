# LIMITATIONS

> **Round 3, 2026-10-10.** Numbers about approved drugs are **provisional** until `drug_links/MANUAL_REVIEW.csv` is reviewed. A strict verification rule now sits beside the permissive one; "no eQTL data" has two named sub-types (autosomal, not tested; X-linked, not analysed) and the autosomal-only lower bound is reported with the overall one. Current figures: `KEY_NUMBERS.md`.

Each limitation states what the data show about its size. Numbers are from `KEY_NUMBERS.md`;
reproduce with `python manuscript_prep\analysis\classify.py` (tags [K..]). The project's earlier
limitations file (`deliverables/cv/LIMITATIONS.md`) covers the colocalization work, which is no
longer part of the paper's claims.

## 1. One resource, one tissue, total mRNA

Every status comes from eQTLGen phase I whole blood. "No eQTL data" means "not in this
resource", not "cannot be instrumented anywhere".

- 180 of 182 autosomal no-eQTL-data targets have GTEx whole-blood TPM below 1 [K10]. A larger
  blood study would not recover them.
- The other tissues examined do not fill the gap: 28 of the 217 have a genome-wide significant
  top cis-eQTL in heart, aorta or liver, and those figures are upper bounds (± 1 Mb, unclumped,
  n = 269–492) [K08].
- Splicing, protein-level regulation, cell-type-specific and stimulus-dependent effects are
  invisible to total mRNA in resting blood. IL6R is the one worked example; how many
  instrumented targets share the problem was not measured.
- pQTL resources were not used. UK Biobank PPP needs credentials; no pQTL coverage figure is
  available from this project.

## 2. eQTLGen did not analyse the X chromosome

None of the 19,250 tested genes is on X, Y or the mitochondrial genome [K10]. All 27 X-linked
and 2 mitochondrial targets in the primary set are therefore "no eQTL data"; 12 of the 27 and
both mitochondrial targets are expressed in blood (relabelled 2026-10-11; previously "29 X-linked").
That is 13 % of the lower bound. They are reported as a named sub-type, and the
autosomal-only lower bound, 188 of 697 (27.0 %), is given beside the overall one everywhere.
Two launched targets are affected (AVPR2, GLA).

## 3. European ancestry

eQTLGen is predominantly European-ancestry; the exact composition is **not in the repository**.
LD for clumping is 1000 Genomes EUR (503 samples). Whether a target has an instrument in other
ancestries was not examined. Status does not depend on LD, so the classification itself is not
sensitive to the reference panel.

## 4. Cis-only, common variants, one threshold

- Trans-eQTLs and rare variants are out of scope.
- The threshold p < 5 × 10⁻⁸ is frozen. 26 of the 63 "no qualifying instrument" targets have a
  cis-eQTL at FDR < 0.05 that misses it [K04]. A laxer threshold would move them, which is why
  the paper reports a lower and an upper bound.
- The window matters little: widening from ± 100 kb to ± 1 Mb changes one target (IL1B) [K07].
  The ± 1 Mb window is not fully covered, because eQTLGen tested variants within 1 Mb of the
  gene centre.

## 5. Dependence on Minikel phase labels

- "Launched" is Minikel's label for a target–indication pair, taken as published.
- It disagrees with ChEMBL for 23 primary-set targets: 8 launched targets have no approved
  cardiovascular drug in ChEMBL (AGT, EGFR, ITGB2, LPL, PPARG, RXFP1, SLC12A2, TNF), and 15
  Phase III targets do have one (HCN4 among them) [K06]. Labels were not edited.
- The label attaches to the target, so a non-selective drug marks several genes. Five of the
  20 launched no-eQTL-data targets are α-adrenoceptor subtypes and two are endothelin
  receptors. **Counts are genes, not drugs.**
- Four programmes with conflicting labels were excluded.

## 6. Drug verification depends on the rule, and is provisional

- Two rules are reported. Permissive: any ChEMBL link, any approved indication under MeSH C14
  or a lipid term. Strict: a direct link to the gene and an indication that sits only in
  cardiovascular MeSH branches. They give 58 and 35 verified launched targets [K06].
- The strict rule rejects 16 launched targets only because ChEMBL records their drugs against a
  receptor family or complex group (α-adrenoceptor and endothelin subtypes, sGC subunits). It
  is therefore conservative for non-selective drugs.
- Neither rule requires the approved indication to be the one Minikel labels. IL1B is verified
  through pericarditis; its Launched label is for atherosclerosis.
- The numbers are provisional until `drug_links/MANUAL_REVIEW.csv` has been reviewed.
- ChEMBL's indication rows include weak entries. GRIN1 is "verified" by amantadine for
  intracranial arteriosclerosis. Rows were not filtered.
- Drugs annotated to a receptor family or complex are linked to every member gene.
- Antisense and RNAi drugs were linked through a hand-made map of 62 mechanisms.
- Two launched targets (LPL, SLC12A2) have no mechanism row in ChEMBL, so they count as
  unverified for lack of a link, not for lack of a drug.
- Only targets at Phase III or launched were checked.

## 7. The cardiovascular set

- The rule follows Minikel's tag. It brings in 15 terms that are not primarily cardiac
  (170 pairs; Diabetic Retinopathy and Diabetic Foot are 136), kept by decision. 46 targets are
  in the set only through one of them.
- Three launched targets come from such additions (GHR, GLA, PDE4C). They fail the strict
  rule and are never named as examples.
- The seven lipid terms are in through a token rule, because Minikel tags them metabolic.
- Brain Ischemia and Stroke have identical target lists.
- In Table 3, "Other" holds 46 of 93 terms, and a target is counted in every area it appears
  in, so area rows overlap and cannot be compared by a simple test.

## 8. What was pre-registered and what was not

- Pre-registered: the instrument threshold, the ± 100 kb window, the controls, the Class C
  refutation criterion.
- **Post hoc, decided 2026-10-10:** the cardiovascular set rule, the two sub-type names and the
  bound framing, the window sensitivity analysis, the gene-record rule, ChEMBL verification, the
  indication areas.
- The pre-registration is a dated document in a folder. The project is not under version
  control and was not lodged with a registry, so its timestamps cannot be shown independently.
- The exposure resource was changed from GTEx liver to eQTLGen after all three controls failed
  in liver. The change was recorded before the re-run, but it was made with those results
  known.

## 9. Validation of the MR pipeline is thin

- One positive control passes (HMGCR → LDL-C). IL6R → CRP fails on direction. PCSK9 has no
  instrument.
- The rule-based replacement (concordance of MR direction with drug mechanism) was null:
  287 of 571, 50.3 %, p = 0.93 (`PART_C_CONCORDANCE.md:11-12`; not recomputed here).
- This does not touch the classification, which needs only a lookup. It does mean the paper
  should make no claim that rests on an MR estimate.

## 10. The Class C result

- The refutation compares coloc.abf PP4 values that are all below 0.15. It shows the predicted
  decline does not occur; it does not show colocalization works at these loci.
- Two sample sizes only, so no statement about the shape of the relationship.
- 31 of the 50 loci are cardiovascular targets; the test is not specific to them.
- The project found its own colocalization implementation unreliable at HMGCR. The test is
  paired, which removes outcome-side weakness, but a reviewer may still discount it.

## 11. Unresolved and unmeasured

- Three microRNA targets have no usable gene record (MIR29B1, MIR92A1, MIR195). At least four
  other microRNA targets (MIR122, MIR143, MIR145, MIR21) are classed "no eQTL data"; whether
  eQTLGen's platforms could measure microRNAs at all is **unverified**.
- ADRA2B has no GTEx or Human Protein Atlas record under its GRCh37 ID.
- Outcome GWAS availability is a separate barrier not counted in the bounds: 404 of 771
  primary pairs with a matching European GWAS have a file that cannot be range-queried [K09].
- Retrieval URLs for GTEx, Human Protein Atlas and Open Targets were not logged; licences were
  not re-verified; installed package versions may differ from those used in August.

## 12. What the atlas is not

It records whether a question can be asked in one resource. It says nothing about whether any
target affects any disease, and it is not a guide to clinical decisions.
