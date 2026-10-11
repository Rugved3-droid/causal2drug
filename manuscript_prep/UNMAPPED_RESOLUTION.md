# UNMAPPED_RESOLUTION — gene-record mapping for the 12 unmapped CV targets

> **Update 2026-10-10: fixes approved and applied.** The same rule was then run on all 60 unmapped symbols: 45 resolved, 15 unresolved (`tables/tableS8_gene_record_resolution.csv`). In the primary set (726 targets) 15 were unmapped and **3 stay unresolved: MIR29B1, MIR92A1 and MIR195**. The table below covers the 12 in the original 630-target set.

**Status: resolved 10 of 12, 2026-10-10. Provisional until you approve; the atlas file is
unchanged.** Covers the 12 unmapped targets of the original 630 set. Targets added by the
extended set will be put through the same rule when that set is approved.

Reproduce (HGNC responses are cached, so a re-run is offline):

```powershell
cd D:\causal2drug
python manuscript_prep\verify\resolve_unmapped.py      # writes unmapped_resolution.csv
python manuscript_prep\verify\gtf_lookup_unmapped.py   # the raw GTF records
```

## The rule, fixed before any eQTLGen state was looked up

1. Ask HGNC for the approved record of the Minikel symbol and take its Ensembl gene ID.
2. If that ID is a `gene` record in the pinned Ensembl GRCh37 release-87 GTF, the target maps
   to it.
3. Otherwise the target stays unresolved and is named in the paper.

No coordinate, symbol or ID was chosen by hand. State is then assigned by the project's
existing definition (`DEFINITIONS.md`, Part 1), unchanged.

## Why they were unmapped

Two causes, both in `deliverables/d1b/gene_coords_unmapped.tsv`:

* **Two GTF records share the symbol** (ACE, ITGB3, TLR9). In each case one record is the
  merged Ensembl/Havana gene and the other is a Havana-only duplicate that overlaps it. HGNC
  names the merged record. The duplicate was not tested by eQTLGen in any of the three cases.
* **The symbol is newer than the 2017 GTF** (ANGPTL8, CCN1, GRK2, GUCY1A1, GUCY1B1, VEGFD,
  FIGNL2). The gene is in the GTF under its earlier name.

## Result for each target

| target | CV phase | cause | fix: Ensembl ID (name in GTF) | eQTLGen | state |
|---|---|---|---|---|---|
| **ACE** | Launched | 2 records | ENSG00000159640 (ACE); duplicate ENSG00000264813 not tested | 37 variants at p < 5e-8 in window, best 1.6 × 10⁻¹⁷ | **INSTRUMENTED** |
| **ITGB3** | Launched | 2 records | ENSG00000259207 (ITGB3); duplicate ENSG00000259753 not tested | 215 variants, best 3.1 × 10⁻³⁴ | **INSTRUMENTED** |
| **GUCY1A1** | Launched | renamed | ENSG00000164116 (GUCY1A3) | 153 variants, best 7.6 × 10⁻³⁷ | **INSTRUMENTED** |
| **GUCY1B1** | Launched | renamed | ENSG00000061918 (GUCY1B3) | tested, 0 variants in window | **TESTED-NULL** |
| TLR9 | Phase I | 2 records | ENSG00000239732 (TLR9); duplicate ENSG00000173366 not tested | 17 variants, best 2.0 × 10⁻¹⁰ | INSTRUMENTED |
| FIGNL2 | Preclinical | renamed | ENSG00000261308 (clone name RP11-923I11.7) | 69 variants, best 3.6 × 10⁻⁵⁰ | INSTRUMENTED |
| GRK2 | Preclinical | renamed | ENSG00000173020 (ADRBK1) | tested, 0 variants in window | TESTED-NULL |
| VEGFD | Phase II | renamed | ENSG00000165197 (FIGF) | not tested | ABSENT |
| ANGPTL8 | Phase I | renamed | ENSG00000130173 (C19orf80) | not tested | ABSENT |
| CCN1 | Preclinical | renamed | ENSG00000142871 (CYR61) | not tested | ABSENT |
| **MIR29B1** | Phase I | — | HGNC's Ensembl ID is not a gene record in GRCh37.87 | — | **UNRESOLVED** |
| **MIR92A1** | Phase I | — | HGNC's Ensembl ID is not a gene record in GRCh37.87 | — | **UNRESOLVED** |

Totals: 5 instrumented, 2 tested-null, 3 absent, 2 unresolved.

For the three absent genes I also checked that neither the current nor the previous symbol
appears among eQTLGen's tested gene symbols; none does.

## The two that stay unresolved

**MIR29B1 and MIR92A1** are microRNA genes. Their current Ensembl IDs do not exist in the
GRCh37 release-87 GTF, so rule 2 fails. Both will be listed by name in the paper as
unclassified. Neither is at Phase III or launched (both Phase I).

A related point for the limitations section: the absent group already contains four other
microRNA targets (MIR122, MIR143, MIR145, MIR21). My understanding is that eQTLGen's expression
platforms do not measure mature microRNAs, which would make "absent" for a microRNA a platform
property and not a finding about blood expression. That is from general knowledge and is
**UNVERIFIED** against the eQTLGen methods; it should be checked before the paper says it.

## What this changes (original 630 set, provisional)

| | as shipped | after the fixes |
|---|---|---|
| Absent | 177 | 180 |
| Tested-null | 55 | 57 |
| Instrumented | 386 | 391 |
| Unmapped | 12 | 2 |
| Launched targets with a state | 48 of 52 | 52 of 52 |
| Launched, absent + tested-null | 21 of 52 (40.4 %) | 22 of 52 (42.3 %) |

**ACE can be instrumented in blood.** That answers the reviewer question raised in Phase 1.
The launched count rises by one because GUCY1B1, a subunit of the soluble guanylate cyclase
that riociguat and vericiguat act on, is tested-null.

## Caveats

* LD clumping was not re-run for the five newly instrumented genes, so the number of
  independent instruments is not reported for them. The state does not depend on clumping
  (`DEFINITIONS.md`).
* The same two failure causes affect the other 48 unmapped targets outside the CV set
  (60 genome-wide). They were not resolved here. Genome-wide comparison figures (for example
  35.0 % across all 2,517 targets) still use the unresolved mapping.
* ADRA2B is mapped (ENSG00000222040) and absent from eQTLGen, but has no GTEx or Human Protein
  Atlas values because that ID is not in those resources. Its state is not in question; its
  expression columns are empty.
