# Deliverable 1b — instrument yield census (eQTLGen phase I, whole blood)

Run 2026-08-26/27 UTC. Rules applied unchanged from `deliverables/d1/PREREGISTRATION.md`
(R1 window, R5 clumping) under DECISION 6 (threshold frozen at p < 5e-8) and DECISION 7 /
AMENDMENT R3a (resource changed to eQTLGen phase I). **No MR was run. No control was re-run.**

---

## Headline

**1,576 of 2,517 Minikel targets (62.6 %) have at least one independent cis-eQTL instrument.**

Reproduce:

```bash
cd /d/causal2drug
python src/d1b_scan_eqtlgen.py      # significant file -> hits + significant-gene universe
python src/d1b_tested_universe.py   # FULL file -> 19,250 tested genes
python src/d1b_gene_coords.py       # GRCh37 gene bodies from Ensembl 87 GTF
python src/d1b_clump.py             # LD clump vs 1000G phase3 EUR (resumable)
python src/d1b_clump_delta.py       # reconcile onto the ENSG join
python src/d1b_census.py            # every number below
```

---

## 1. Sources retrieved

All logged in `logs/retrieval_log.csv` (URL, timestamp, bytes, sha256).

| file | bytes | sha256 (first 16) | retrieved UTC |
|---|---|---|---|
| eQTLGen significant cis-eQTL (FDR<0.05) | 322,775,879 | `8d963046d7b74cf3` | 2026-08-26T23:15:04Z |
| eQTLGen **full** cis-eQTL | 3,880,187,026 | `1133a89fe33b795d` | 2026-08-26T23:16:56Z |
| Ensembl GRCh37.87 GTF | 42,400,095 | `56bc520f5e9cbc14` | 2026-08-26T23:23:12Z |
| 1000G phase 3 sample panel | 55,156 | `b4023dc6ee2d62ee` | 2026-08-26T23:19:58Z |

Release: eQTLGen phase I, files dated **2019-12-11** (2-cohort filter), retrieved from
`molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/`. Build **GRCh37**.

**Neither eQTLGen file is position-sorted or tabix-indexed** — the directory contains no `.tbi`.
Range queries are therefore impossible; both required **full download**. The full 3.88 GB file was
not optional: the FDR<0.05 file contains only genes that already have a significant eQTL, so it
cannot distinguish "absent from eQTLGen" from "tested but no instrument".

The GTF was verified against Ensembl's **own published `CHECKSUMS`** file — BSD sum `25146 41407`,
matching exactly. This is the only source in the project so far with a genuine publisher-side
checksum to compare against.

Two independent cross-checks passed:
* full file = **127,341,798 rows**, **19,250 distinct tested genes** — exactly the figure eQTLGen's
  own README states ("19,250 genes that showed expression in blood were tested").
* significant file: 10,507,664 rows, min p = 3.27e-310, **max p = 2.019e-05**, confirming it is a
  strict superset of the p<5e-8 set (7,844,064 rows below threshold).

## 2. Target list

`data/derived/pp_dedup_frozen.tsv` (DECISION 1: dedup by `ti_uid`, four conflicting-outcome
programs dropped) → 29,456 pairs, **2,517 unique targets**.

## 3. Gene identifier mapping

**Method (primary): Minikel symbol → Ensembl gene ID (pinned GRCh37 release-87 GTF) → eQTLGen
`Gene` (ENSG) column.** A symbol maps only if it matches exactly one `gene` feature on a primary
chromosome. Zero matches, or more than one, is UNMAPPED. Nothing is aliased or guessed.

ENSG was chosen over symbol matching because eQTLGen's `GeneSymbol` column is Ensembl v71 (~2013)
nomenclature. I tested the symbol join as a comparator rather than assuming: it **rescued 10** and
**lost 8** targets relative to the ENSG join, so nomenclature drift is real but *not* the main
driver of non-mapping. The ENSG join also disambiguates four symbols that map to >1 ENSG inside
eQTLGen (`F11R`, `IFNAR2`, `KLRK1`, `S1PR2`) but to exactly one gene record in the GTF.

**The two failure modes, kept strictly separate as instructed:**

| category | N | % of 2,517 |
|---|---|---|
| **UNMAPPED** — no unique GRCh37 gene record | **60** | 2.4 % |
| **ABSENT from eQTLGen** — never tested | **652** | 25.9 % |
| tested, but **NO** cis-eQTL at p<5e-8 in window | **229** | 9.1 % |
| **≥1 cis-eQTL at p<5e-8 in window** | **1,576** | 62.6 % |

60 + 652 + 229 + 1,576 = 2,517 — the partition is asserted in code, not by hand.

Of the 60 UNMAPPED: 49 absent from the GTF (newer HGNC symbols post-dating release 87, e.g.
`CCN1`/`CCN2`, `ADGRE1`, `ATP5F1D`, `CGAS`), 11 ambiguous (two gene records, e.g. `ACE`, `CFB`,
`AQP1`). Lists in `gene_coords_unmapped.tsv`, `mapping_methodB_absent.txt`,
`mapping_methodB_tested_no_instrument.txt`.

## 4. Census

LD panel: **1000 Genomes phase 3, GRCh37, EUR, 503 samples**
(`integrated_call_samples_v3.20130502.ALL.panel`, `super_pop == EUR`). GRCh37 throughout, matching
eQTLGen — **no liftover was involved in this census**.

| quantity | N |
|---|---|
| targets with ≥1 cis-eQTL p<5e-8 in gene body ±100 kb | **1,576** |
| targets with ≥1 **independent** instrument (r²<0.01 clumped) | **1,576** |
| targets with NO significant cis-eQTL although tested | 229 |
| targets absent from eQTLGen entirely | 652 |
| targets unmapped | 60 |

**No target was lost at the clumping step**: every target with ≥1 significant hit retained ≥1
instrument. Clumping only reduces the *number* of instruments, never below one.

**Independent instruments per target** (N = 1,576):

| statistic | value |
|---|---|
| median | **3** |
| IQR | [1, 5] |
| min / max | 1 / 30 |
| ≥2 instruments | 1,116 (70.8 %) |
| **≥3 instruments** | **818 (51.9 %)** |
| ≥5 instruments | 452 (28.7 %) |
| ≥10 instruments | 128 (8.1 %) |

Only about half of instrumented targets reach the ≥3 instruments that MR sensitivity analyses
(MR-Egger, weighted median) need.

**Caveat on what "independent" means here.** `r² < 0.01` is inherited unchanged from D1 R5. It is a
*strict independence* criterion, so it retains any SNP near-uncorrelated with those already kept;
in dense 200 kb windows with very strong eQTLGen signals this yields up to 30 instruments at a
single cis locus. Verified directly: among kept SNPs for CBS, max pairwise r² = 0.00967 while the
lead SNP had r² = 1.000 with six of its neighbours — so the computation is correct, and the counts
are a property of the rule. **Thirty near-uncorrelated SNPs at one locus are not thirty independent
natural experiments.** This needs a pre-registered decision before D2 MR; I have not changed it.

## 5. Propagation to pairs

| quantity | N |
|---|---|
| total pairs (frozen set) | 29,456 |
| pairs whose target has ≥1 independent instrument | **19,754 (67.1 %)** |

**DISTINCT TARGET counts** (not pair counts), by highest clinical phase reached (`ccatnum`;
coding verified from data: 1=Preclinical … 5=Launched, no missing values):

| phase | all targets | surviving | retention |
|---|---|---|---|
| any phase | 2,517 | **1,576** | 62.6 % |
| Phase II+ | 1,306 | **831** | 63.6 % |
| Phase III+ | 704 | **431** | 61.2 % |
| launched | 450 | **255** | 56.7 % |

Retention is roughly flat across phases, declining mildly toward launched (62.6 % → 56.7 %).

## 6. Blood-expression penalty by therapeutic area

Area mapping uses **Minikel's own files** (`areas.tsv` × `indic_topl_match.tsv`), not a mapping of
my own — all 1,014 MeSH indications map with zero misses, and pair counts sum to 29,456 exactly.

Distinct targets, surviving / all:

| area | all | surviving | yield | Phase II+ | Phase III+ | launched |
|---|---|---|---|---|---|---|
| infection | 238 | 180 | **75.6 %** | 39/62 | 17/30 | 10/22 |
| hematology | 321 | 236 | 73.5 % | 126/170 | 46/62 | 24/33 |
| dermatology | 455 | 320 | 70.3 % | 153/223 | 55/75 | 35/48 |
| respiratory | 534 | 369 | 69.1 % | 173/256 | 67/103 | 16/29 |
| immune | 805 | 549 | 68.2 % | 221/328 | 97/137 | 57/81 |
| musculoskeletal | 248 | 169 | 68.2 % | 48/80 | 25/38 | 11/15 |
| oncology | 1,573 | 1,070 | 68.0 % | 470/655 | 179/251 | 97/140 |
| digestive | 599 | 398 | 66.4 % | 184/287 | 54/88 | 20/40 |
| ophthalmology | 374 | 246 | 65.8 % | 113/173 | 43/65 | 23/29 |
| congenital | 332 | 212 | 63.9 % | 104/161 | 39/70 | 16/41 |
| metabolic | 511 | 319 | 62.4 % | 116/195 | 48/84 | 29/53 |
| cardiovascular | 594 | 368 | 62.0 % | 136/215 | 55/94 | 23/44 |
| signs/symptoms | 1,051 | 644 | 61.3 % | 232/423 | 99/193 | 54/113 |
| neurology | 534 | 327 | 61.2 % | 89/171 | 33/72 | 13/37 |
| psychiatry | 459 | 278 | 60.6 % | 96/193 | 40/102 | 13/52 |
| other | 547 | 330 | 60.3 % | 128/230 | 57/101 | 27/60 |
| endocrine | 545 | 327 | 60.0 % | 125/226 | 31/58 | 19/37 |

**Stated plainly: no therapeutic area is eliminated.** The yield spread is **60.0 % – 75.6 %, a
range of 15.6 percentage points**. The ordering matches the expected biology — infection,
haematology and immune at the top, neurology, psychiatry and endocrine at the bottom — but the
penalty is far milder than "whole blood only" would suggest. Oncology (68.0 %) and cardiovascular
(62.0 %) sit mid-pack rather than being penalised.

The reason the penalty is muted: targets are shared across many indications, so area-specific
target sets overlap heavily (the per-area target counts sum to far more than 2,517).

**Where the real constraint bites is absolute N at late phases, not the yield rate.** Surviving
distinct targets at *launched* fall to 10 (infection), 11 (musculoskeletal), 13 (neurology),
13 (psychiatry) and 16 (congenital, respiratory). Any leave-one-therapeutic-area-out design at
launched will be operating on double-digit target counts in most areas — that is a Deliverable 4
question and is not answered here.

## 7. Control status (lookup, not a re-run)

| control | in eQTLGen? | hits p<5e-8 in window | independent instruments | best p | best SNP |
|---|---|---|---|---|---|
| **HMGCR** | yes | 277 | **1** | **1.21e-50** | rs6453133 |
| **IL6R** | yes | 226 | **2** | **4.09e-34** | rs10908839 |
| **PCSK9** | **NO** | 0 | — | — | — |

**PCSK9 is absent from eQTLGen phase I entirely.** Verified by direct scan of the full
127,341,798-row file: **zero rows** mention `ENSG00000169174` or the symbol `PCSK9`. Other
PCSK-family genes (PCSK4/5/6/7) are present, so this is not a nomenclature artefact. PCSK9 is not
appreciably expressed in whole blood — it is a hepatocyte-secreted protein.

**Advance read on the D1 re-run:** HMGCR and IL6R will have strong instruments (p = 1e-50 and
4e-34, versus the 1.3e-2 and 6.4e-3 that GTEx liver produced). **The PCSK9→LDL control cannot be
run against this resource at all** — not a weak-instrument failure but an absence. That decision is
yours; I have not substituted a control or a resource.

## 8. Three retrieval bugs found and fixed — two were silent

These were in my own code, found by checking rather than by error messages. They matter because two
produced **wrong numbers with no error**.

1. **Collapsed chunk spans.** 1000G VCFs scatter a cis-window across up to 17 disjoint tabix index
   chunks (structural variants). Taking min→max across them requested **491 MB for a 230 kb
   window** (vs 2.5 MB when the window lands in one chunk). Cost-only — extra records are filtered
   by position — but it hung the run. Fixed by fetching each chunk separately.

2. **Unsound upper bound — the serious one.** I used the tabix linear index as both bounds. It is
   valid only as a **lower** bound: a long record overlapping a later 16 kb interval drags that
   entry backwards, so the "upper bound" can fall below data genuinely belonging to the query.
   Result: **silent under-retrieval**. SLC2A9 recovered **0 of 1,226** hit SNPs; PTGES 1 of 335.
   Before the fix, 3.2 % of targets had <80 % of their hits found in the LD panel and 5.9 % had
   80–95 %; after, **0 % have zero, 96.0 % are ≥95 %**. All results computed before this fix were
   quarantined (`clump_ckpt_BUGGY_truncated/`, `clumped_instruments*_BUGGY.tsv`) and recomputed
   from scratch — they under-counted instruments.

3. **Records split across BGZF blocks.** The incremental reader resumed at block boundaries,
   silently dropping any record straddling one (PPT1: 484 vs 485). Fixed with a carry-over buffer;
   verified to give identical results at 128 KB, 512 KB and 2 MB slice sizes.

**Consequence for Deliverable 1 that you should know:** bug #2 lives in `tabix_http.py`, which
produced the D1 region extracts. D1's **conclusion is unaffected** — the BLOCKER rested on the
eQTLGen-independent permutation file (`p_beta`), which never passed through my retrieval code — but
**D1's extract row counts may be undercounts and should be treated as UNVERIFIED** until re-run.
`tabix_http.py` now carries the fixed bound; `d1b_clump.py` additionally has per-chunk fetching,
early stop, and the carry buffer. Porting those three to `tabix_http.query()` is outstanding.

## 9. Operational notes

* The 1000genomes host caps aggregate throughput at **~450 KB/s regardless of connection count**
  (measured: 1 worker 437 KB/s, 8 workers 482 KB/s; 16 workers causes dropped connections). Total
  LD volume was ~1.8 GB, so this step is bandwidth-bound and took several hours across restarts.
* Merging overlapping windows per chromosome would save only **14.8 %** (428.7 Mb → 365.3 Mb of
  genome span), so the cost is close to irreducible with this LD source.
* The clumper checkpoints **per target** (`clump_ckpt/chr*.tsv`, fsync'd) and resumes. Background
  tasks in this environment were killed three times mid-run; with per-target checkpointing this
  cost at most one target each time. An earlier per-*chromosome* checkpoint design lost 290
  completed targets to the first kill and made zero net progress on chromosomes longer than the run
  window — that is why the granularity changed.
* **Total bytes transferred for the clump is UNVERIFIED**: the byte counter is per-process and the
  run was restarted ~23 times, so no cumulative figure was retained. The final resumed run reported
  1,880,818 bytes, which covers only that run.
* Ensembl's **GRCh37 REST server returned HTTP 500** (HTML error page) for both GET and POST during
  this run, so gene coordinates came from the pinned release-87 GTF instead — more reproducible
  regardless.

## 10. Files

```
deliverables/d1b/
  DELIVERABLE_1B.md                     this report
  census.json                           every number above
  clumped_instruments_final.tsv         1,576 targets: n_hits, n_in_ld_panel, n_independent,
                                        best_p, best_snp, independent_snps (rsIDs)
  clump_ckpt/chr*.tsv                   per-target checkpoints (resumable)
  eqtlgen_tested_universe.tsv           19,250 tested genes
  eqtlgen_sig_gene_universe.tsv         16,987 genes with >=1 FDR<0.05 cis-eQTL
  eqtlgen_hits_by_ensg.tsv              885,846 hit rows, ENSG join (primary)
  eqtlgen_hits_targets.tsv              888,379 hit rows, symbol join (comparator)
  gene_coords_grch37.tsv                2,457 resolved gene bodies
  gene_coords_unmapped.tsv              60 unmapped, with reason
  mapping_methodB_absent.txt            652 absent from eQTLGen
  mapping_methodB_tested_no_instrument.txt   tested but no instrument
  clump_ckpt_BUGGY_truncated/           QUARANTINED pre-fix results, do not use
  clumped_instruments*_BUGGY.tsv        QUARANTINED pre-fix results, do not use
```

**Deliverable 2 has not been started. No MR was computed. No control was re-run.**
