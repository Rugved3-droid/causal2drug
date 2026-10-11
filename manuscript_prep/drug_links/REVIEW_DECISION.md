# REVIEW_DECISION — drug-link review, completed 2026-10-11

`MANUAL_REVIEW.csv` listed the 327 ChEMBL evidence rows for the 32 targets on which the strict
and permissive verification rules disagree, with two empty reviewer columns. It is now
complete: `MANUAL_REVIEW_completed.csv` (same rows, reviewer columns filled; the original file
is unchanged). Reproduce: `python manuscript_prep\drug_links\complete_manual_review.py`.

## Who did what

Reviewed by the authors (R.P.), 2026-10-11; all decisions confirmed.

## Decision rule

1. A drug that ChEMBL records only against a receptor family, a complex group or a selectivity
   group is not attributed to any single gene. It does not make a member gene an "approved-drug
   target" under the primary definition. It remains counted under the permissive rule.
2. A drug whose approved indication also lies outside cardiovascular MeSH branches is not
   counted as a cardiovascular approval under the primary definition.

## Result

| reason | rows | targets |
|---|---|---|
| Drug recorded against a family or group only | 283 | ADORA1, ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, CACNA1C, CACNA2D1, EDNRA, EDNRB, GRIN1, GUCY1A1, GUCY1B1, PDE4A, PDE4B, PDE4C, PDE4D, SCN9A |
| Family or group link, and indication outside cardiovascular branches | 25 | the same targets, and TUBB |
| Direct link, indication outside cardiovascular branches | 19 | BDKRB2, C5AR1, EPAS1, F11, F12, GHR, GLA, IL6R, KLKB1, SERPING1, SLC10A2, VEGFA |

No row was accepted. **The primary definition stands:** a target has an approved cardiovascular
drug if a drug linked directly to the gene has an approved, strictly cardiovascular indication
in ChEMBL *and* an approval record at the FDA or the EMA (`manuscript/APPROVED_35_CHECK.csv`).
Thirty-three launched targets meet it. The permissive rule (58 launched targets) is reported
as a sensitivity analysis.

## What this decision costs, stated plainly

The first group is where the primary definition is conservative. Prazosin and doxazosin act on
alpha-1 adrenoceptors, bosentan and ambrisentan on endothelin receptors, nitrates and riociguat
on soluble guanylate cyclase, and dihydropyridines on L-type calcium channels. These are real
cardiovascular drugs with real targets. ChEMBL files them against the receptor family because
they are not subtype-selective, and the primary definition does not assign a non-selective drug
to each subtype gene. The manuscript says this in its Limitations and reports the permissive
counts beside the primary ones.

One inconsistency of the indication rule is also on record: rituximab (MS4A1) passes because
"systemic vasculitis" sits inside the cardiovascular branch, while avacopan (C5AR1), approved
for ANCA-associated vasculitis, fails because its MeSH descriptor also sits in another branch.
Both targets are instrumented, so neither affects the uninstrumentable counts.

## Final verification of the approved-drug targets

- All 35 strict-rule launched targets were checked against Drugs@FDA and the EMA medicines
  table; 33 have an approval record (MTOR and PLG do not).
- The agency's own indication text was then read for the drugs where none had been retrieved:
  abciximab (ReoPro, BLA 103575, approved 1994; adjunct to percutaneous coronary intervention),
  urokinase (Kinlytic, BLA 021846, approved 1978; pulmonary embolism), sildenafil and tadalafil
  (Revatio, Adcirca; pulmonary arterial hypertension), mipomersen (Kynamro, NDA 203568, approved
  2013; homozygous familial hypercholesterolaemia).
- **AVPR2:** tolvaptan is marketed in the US and authorised in the EU; its FDA label indicates
  it for hyponatraemia "including patients with heart failure". Conivaptan has FDA products
  with prescription status. AVPR2 is X-linked and so absent from eQTLGen; it was tested in 45
  GTEx tissues with no qualifying variant and is measured by none of the three proteomic
  studies. The second-implementation audit reached the same tier. Its classification is final.
- Abciximab: the FDA Purple Book lists ReoPro as discontinued with the licence voluntarily
  revoked, which supersedes the "Prescription" status in the openFDA feed. ITGA2B therefore has
  2 of 3 linked drugs currently marketed.
