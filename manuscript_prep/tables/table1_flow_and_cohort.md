**Table 1. Derivation of the cardiovascular target set and its composition**

|  | N |  |  |  |
|---|---|---|---|---|
| Rows in the source target-indication table | 29476 |  |  |  |
| Distinct target-indication pairs | 29460 |  |  |  |
| Pairs after removing 4 programmes with conflicting phase labels | 29456 |  |  |  |
| Distinct targets / indications in those pairs | 2517 / 1014 |  |  |  |
|  | Primary set | Sensitivity 1 (original 630) | Sensitivity 2 | Sensitivity 3 (MeSH C14 only) |
| Cardiovascular indications (MeSH terms) | 93 | 44 | 62 | 37 |
| Cardiovascular target-indication pairs | 2128 | 1543 | 1879 | 1364 |
| Cardiovascular targets | 726 | 630 | 663 | 594 |
|   Instrumented | 443 (61.0%) | 391 (62.1%) | 409 (61.7%) | 373 (62.8%) |
|   No eQTL data (lower bound) | 217 (29.9%) | 180 (28.6%) | 191 (28.8%) | 164 (27.6%) |
|     autosomal gene not tested | 188 (25.9%) | 157 (24.9%) | 167 (25.2%) | 144 (24.2%) |
|     X-linked, not analysed by eQTLGen | 27 (3.7%) | 21 (3.3%) | 22 (3.3%) | 18 (3.0%) |
|     mitochondrial, not analysed by eQTLGen | 2 (0.3%) | 2 (0.3%) | 2 (0.3%) | 2 (0.3%) |
|   No qualifying instrument | 63 (8.7%) | 57 (9.0%) | 60 (9.0%) | 55 (9.3%) |
|   Uninstrumentable, upper bound | 280 (38.6%) | 237 (37.6%) | 251 (37.9%) | 219 (36.9%) |
|   Unresolved gene record | 3 (0.4%) | 2 (0.3%) | 3 (0.5%) | 2 (0.3%) |
|   Autosomal-only lower bound | 188 / 697 (27.0%) | 157 / 607 (25.9%) | 167 / 639 (26.1%) | 144 / 574 (25.1%) |
| Highest phase reached in the set |  |  |  |  |
|   Preclinical | 323 | 297 | 300 | 294 |
|   Phase I | 91 | 90 | 93 | 85 |
|   Phase II | 171 | 140 | 155 | 121 |
|   Phase III | 75 | 51 | 60 | 50 |
|   Launched | 66 | 52 | 55 | 44 |
| Launched, approved CV drug verified in ChEMBL: strict rule | 35 | 32 | 35 | 26 |
| Launched, drug ever approved by FDA or EMA for a CV indication | 33 | 30 | 33 | 24 |
| Launched, approved CV drug verified in ChEMBL: permissive rule | 58 | 46 | 49 | 40 |

Primary set: every indication that Minikel et al. tag as cardiovascular, plus indications matching eight lipid/atherosclerosis tokens. Sensitivity 1: the originally analysed set (MeSH branch C14 plus the lipid tokens). Sensitivity 2: sensitivity 1 plus cardiovascular-tagged terms filed under signs and symptoms. Sensitivity 3: MeSH branch C14 only.

Phase is the highest Minikel phase over the target's indications in the same set. Percentages are of all targets in the set.

Instrumented: at least one eQTLGen whole-blood cis-eQTL at P<5×10⁻⁸ within the gene body ±100 kb. No eQTL data: the gene is not among the 19,250 genes eQTLGen tested; its sub-types are autosomal genes not tested, and X-linked and mitochondrial genes, which eQTLGen did not analyse (no target is Y-linked). The autosomal-only lower bound is autosomal no-eQTL-data targets over targets not on X, Y or MT. No qualifying instrument: tested, but no variant meets the criterion. The upper bound is no eQTL data plus no qualifying instrument and is reported only as a bound.

The evidence rows on which the two ChEMBL rules disagree were reviewed and the strict rule was retained. Ever approved by FDA or EMA: the strict rule, and at least one of the target's strict-rule drugs has an approval record in Drugs@FDA or the EMA medicines table; this is the definition of 'approved' used in the text. Two of the 33 rest on an indication that is cardiovascular only in a broad sense: NR3C1 (dexamethasone, prednisolone and prednisone; ChEMBL approved indication myocarditis) and MS4A1 (rituximab; systemic vasculitis). Both are instrumented, so they enter the denominator only. Strict: an approved drug linked directly to the gene (single protein, specific complex or mRNA) has an approved indication that sits only in cardiovascular MeSH branches. Permissive: any ChEMBL link, any approved indication under MeSH C14 or a lipid term.
