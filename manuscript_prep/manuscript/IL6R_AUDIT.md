# IL6R_AUDIT — the IL6R → C-reactive protein control

Audited 2026-10-11. **The control is not reclassified. It stays a FAIL on the pre-specified
sign.** This file reports what the audit found.

Reproduce: `python manuscript_prep\manuscript\il6r_audit\il6r_audit.py`
Output: `il6r_audit\il6r_audit_output.txt` (tags [A01]–[A10]). Inputs read, not changed:
`deliverables\d1c\IL6R_CRP_provenance.json`, `IL6R_CRP_harmonized_merge.tsv`, the eQTLGen full
file, 1000 Genomes phase 3 (503 EUR, the pipeline's panel) and the CRP file GCST90179146.

## Short answer

1. **Harmonization is correct.** The sign was not produced by an allele error.
2. **The instrument does not tag rs2228145 (Asp358Ala).** r² = 0.16 for the variant the MR
   used, 0.23 for the eQTLGen lead.
3. **The pre-specified sign rested on a premise the data do not support.** The protocol
   expected the lead signal to tag rs2228145 and expected the 358Ala allele to *raise* IL6R
   mRNA. In eQTLGen the 358Ala allele is associated with slightly *lower* IL6R expression and
   is not genome-wide significant.
4. **Every variant examined gives the same sign:** the allele with higher IL6R expression in
   blood has higher CRP.

## 1. What the eQTLGen IL6R instrument is  [A01]–[A03], [A09]

| | variant | GRCh37 | allele assessed | expression z | p |
|---|---|---|---|---|---|
| eQTLGen lead in the window | rs10908839 | chr1:154,430,798 | C | +12.18 | 4.1 × 10⁻³⁴ |
| Variant the MR used | **rs2229238** | chr1:154,437,896 | T (frequency 0.216) | +11.60 | 4.3 × 10⁻³¹ |

IL6R has 4,457 rows in the eQTLGen full file, 841 within the gene body ± 100 kb, and 226 of
those at p < 5 × 10⁻⁸. The lead, rs10908839, has no row at its lifted position in the CRP
file, so the pipeline's rule (the lowest-p variant that survives harmonization) took the
second-ranked variant, rs2229238. That is the rule working as written, not a substitution.

## 2. Allele harmonization of the MR variant  [A08]

| | value |
|---|---|
| eQTLGen | assessed allele **T**, other C, z +11.60; frequency of T 0.216 |
| CRP file | effect allele **C**, other T, β −0.04204, SE 0.00258; frequency of C 0.814, so T 0.186 |
| By hand | the outcome effect allele is the eQTLGen *other* allele, so β for T is +0.04204 |
| Pipeline | status `aligned_flipped`, β for T +0.04204 |

The alleles are C/T, not palindromic, and the two allele frequencies agree (0.216 and 0.186).
The Wald ratio is +0.0420 / +0.1118 = +0.376, as reported. Across all 167 harmonized variants
in the control's merge file, 101 give a positive single-variant ratio and 66 a negative one
[A10]; the lead-variant estimate is not an outlier produced by one row.

## 3. Linkage disequilibrium with rs2228145  [A04], [A05]

rs2228145 is at chr1:154,426,970 (A>C; C encodes 358Ala; C frequency 0.360 in the panel).

| variant | r² with rs2228145 | which alleles travel together |
|---|---|---|
| rs10908839 (eQTLGen lead) | 0.225 | its G allele with 358Ala; so its expression-raising C allele with 358Asp |
| **rs2229238 (MR variant)** | **0.161** | its C allele with 358Ala; so its expression-raising T allele with 358Asp |
| rs4129267 | 1.000 | T with 358Ala |
| rs12126142 | 1.000 | A with 358Ala |
| rs4537545 | 0.935 | T with 358Ala |

The instrument is in weak LD with Asp358Ala. It is neither the missense variant nor a proxy
for it.

## 4. Direction of each variant with IL6R expression and with CRP  [A03], [A06], [A07]

Both associations are expressed for the same allele.

| variant | allele | IL6R expression in blood (eQTLGen z, p) | CRP β (SE) | sign of the ratio |
|---|---|---|---|---|
| rs2229238 (MR variant) | T | +11.60, 4.3 × 10⁻³¹ | +0.0420 (0.0026) | positive |
| **rs2228145 (358Ala)** | C | **−3.17, p = 0.0015** (not significant at 5 × 10⁻⁸) | −0.0835 (0.0020), p = 2.5 × 10⁻³⁶⁵ | positive |
| rs4129267 | T | −4.39, 1.1 × 10⁻⁵ | −0.0836 (0.0020) | positive |
| rs12126142 | A | −4.35, 1.4 × 10⁻⁵ | −0.0863 (0.0021) | positive |
| rs4537545 | T | −6.80, 1.0 × 10⁻¹¹ | −0.0810 (0.0020) | positive |

The 358Ala allele lowers CRP strongly, as expected from the literature the protocol relied
on. But in eQTLGen whole blood it goes with *lower*, not higher, IL6R expression, and only
weakly.

## 5. Why the control was pre-specified as negative

`deliverables\d1\PREREGISTRATION.md`, lines 112–120 and 382, set the expected sign to negative
and called it "declared ambiguous in advance". The reasoning written there: tocilizumab lowers
CRP, which naively implies a positive sign; but the dominant cis signal was expected to tag
rs2228145, whose allele "increases IL6R mRNA and increases receptor shedding", reducing
classical signalling and lowering CRP.

Against that reasoning, the audit shows: the lead signal does not tag rs2228145 (r² 0.16–0.23),
and in this dataset the 358Ala allele is not associated with higher IL6R mRNA.

## 6. What the audit does not settle

- **Whether the rs2229238 CRP association is its own signal.** A rough check: if rs2229238 had
  no effect of its own, LD with rs2228145 alone would give it a CRP β of about
  r × β × √(p₁q₁ / p₂q₂) = 0.401 × 0.0835 × √(0.360 × 0.640 / (0.780 × 0.220)) ≈ 0.039, against
  0.042 observed. So most or all of the CRP association at the instrument could be the
  missense variant seen through LD. This is an approximation from marginal statistics, not a
  conditional analysis, and it uses panel frequencies. If it is right, the MR variant affects
  expression through one mechanism and "affects" CRP through LD with another: an LD-confounded
  instrument, which is the failure the control was meant to expose, though with the opposite
  sign from the one predicted.
- Whether 358Ala raises IL6R mRNA in other tissues or isoform-specific assays. eQTLGen
  measures total gene expression in blood.
- Soluble versus membrane-bound receptor. No protein data were used in this control.

## 7. Bottom line for the manuscript

The IL6R control is reported as a **pre-specified failure with verified harmonization**. It is
consistent with total IL6R mRNA in blood being a poor proxy for a receptor that exists in
membrane-bound and soluble forms, the balance of which Asp358Ala shifts (Ferreira et al.,
*PLoS Genet* 2013, doi 10.1371/journal.pgen.1003444: 358Ala raises soluble IL-6R). In the
manuscript it sits in the paragraph on mismatch between what an instrument measures and what
a drug does.

The IL6R control fails on its pre-specified criterion and is reported as a fail. The audit
found no pipeline error. It found that the expected sign depended on an assumption about which
variant the instrument would tag, and that assumption was wrong. The manuscript should say:
the control failed on the pre-specified sign; allele harmonization was verified; the instrument
is in weak LD (r² = 0.16) with Asp358Ala; and the result illustrates that a strong cis-eQTL can
be confounded by LD with a functional variant acting through a different mechanism.
