# Causal2Drug — frozen decisions and open items

Decisions here are binding on all downstream work. Frozen items are not revisited.

---

## DECISION 1 — pp.tsv deduplication rule — FROZEN 2026-08-26

**Rule.** Deduplicate `pp.tsv` by `ti_uid`. For the four `ti_uid` values carrying conflicting
highest-phase-reached values, **drop all rows for those programs entirely** — no max, no min, no
first-row-wins:

| dropped ti_uid | indication | conflict observed |
|---|---|---|
| `DAO-D022124` | Hyperammonemia | Phase III vs Preclinical |
| `HRH1-D000086382` | COVID-19 | Launched vs Phase II |
| `IFNAR2-D012141` | Respiratory Tract Infections | Launched vs Phase II |
| `MC1R-D010787` | Photosensitivity Disorders | Phase III vs Phase II |

**Rationale.** Two of the four flip between *Launched* and *Phase II*. Any reconciliation rule
(max / min / first) is an outcome-bearing judgment call made by the analyst on the dependent
variable. Dropping four programs out of ~29,460 costs essentially nothing in power and removes the
judgment call from the study entirely.

**Resulting N** — recomputed from raw, not carried from an earlier step:

```bash
python -c "
import collections
h,*d=open(r'D:\causal2drug\data\raw\pp.tsv',encoding='utf-8').read().rstrip('\n').split('\n')
iu=h.split('\t').index('ti_uid'); by=collections.defaultdict(list)
for l in d: by[l.split('\t')[iu]].append(l)
DROP={'DAO-D022124','HRH1-D000086382','IFNAR2-D012141','MC1R-D010787'}
print(len(d), len(by), len({k:v for k,v in by.items() if k not in DROP}))"
# 29476 29460 29456
```

| stage | N |
|---|---|
| raw data rows | 29,476 |
| after dedup by `ti_uid` | 29,460 |
| **after dropping the 4 conflicting programs — FROZEN ANALYSIS SET** | **29,456** |
| rows removed in total | 20 |
| distinct genes / distinct MeSH indications after | 2,517 / 1,014 |

Frozen set written to `data/derived/pp_dedup_frozen.tsv`
(sha256 `5cb7b8015bdc4756f98e55ac787d04d00c3e8a733fb62762cfdac96ad56c153c`, 29,457 lines incl. header).

**Residual note, flagged not decided.** Twelve of the remaining duplicated `ti_uid` keys collapse
under first-row-wins: nine are byte-identical (no effect), but three — `MAPK14-D000086382`,
`TGFB1-D000086382`, `TNFRSF1A-D000086382` — differ in `oto_acat`/`oto_ccat`/`orphan` only. These
columns are **not** the outcome, so the collapse is outcome-neutral, but if `orphan` or any `oto_*`
field is later used as a *feature*, first-row-wins becomes a silent arbitrary choice for those
three. Revisit only at the point such a feature is introduced.

---

## ADDITION 2 — 1000 Genomes access mode — FROZEN 2026-08-26

No full-release download of 1000 Genomes. LD is obtained by HTTP-range tabix query throughout.
Verified in Deliverable 0: a 228 kb cis-window on chr5 cost 2,422,986 bytes = 0.128 % of the
1.90 GB chromosome file, versus 30.20 GiB for the full release.

---

## ADDITION 3 — `.tbi` availability is its own cascade step — FROZEN 2026-08-26

In the Deliverable 4 cascade, "has full public summary statistics" and "is retrievable by region
query" are **separate counted steps**. A study having full summary statistics does **not** imply it
is queryable.

Measured in Deliverable 0:

| population | n probed | indexed | proportion | 95 % CI |
|---|---|---|---|---|
| random sample of all harmonised studies | 200 | 173 | 86.50 % | [81.8 %, 91.2 %] |
| LDL-trait harmonised studies (census) | 294 | 219 | 74.49 % | — (census) |

The LDL census figure is not a sample estimate and carries no sampling error; it also shows the
rate is trait-dependent, so the global 86.5 % must not be applied blindly per trait.

---

## OPEN ITEM 4 — Open Targets release integrity — DEFERRED 2026-08-26

`latest` is byte-identical to release `26.06` (same parquet part UUIDs). The release-level
`manifest.json` for `26.06` reports `"result": "failure"` (`started_at 2026-06-15T16:11:47Z`);
`25.12` also reports `failure`, `26.03` reports `pending`. Individual dataset directories do carry
`_SUCCESS` markers.

**Deferred — not being chased now.** When the tractability and safety gates are reached, pin to a
specific release whose manifest verifies as successful, rather than `latest`. What the
release-level failure actually covers is currently **UNVERIFIED**.

---

## DECISION 5 — Deliverable 1 pre-registration — FROZEN 2026-08-26, before any result was seen

See `deliverables/d1/PREREGISTRATION.md`. Written and timestamped before the analysis was run.

---

## DECISION 6 — instrument threshold is frozen at p < 5e-8 — FROZEN 2026-08-26

**Rule.** The pre-registered instrument threshold stays at **p < 5 × 10⁻⁸**, unchanged. The
*resource* changes (Decision 7); the *threshold* does not.

**Rationale, recorded verbatim as instructed.** GTEx v10 liver yields 2,603 eGenes of 20,540 genes
tested (~13%). That is an instrument-availability failure of the resource, not of the three
controls. Pre-registered instrument threshold stays at p < 5e-8, unchanged. Explicitly noted and
frozen: PCSK9's best liver cis-eQTL was 4.57e-6 with F=22.26. I have now seen that relaxing to 1e-5
would rescue exactly one control. Because I have seen it, the threshold is not negotiable in this
project.

**Why this matters methodologically.** The threshold is now contaminated as a free parameter: its
consequences for a specific control are known. Any later move of it could not be distinguished from
selecting the value that produces a desired result, regardless of intent. Freezing it removes that
degree of freedom permanently. This applies to every deliverable in the project, not only the
positive controls.

---

## DECISION 7 — exposure resource changed to eQTLGen phase I — FROZEN 2026-08-26

**Old rule (R3, 2026-08-26):** GTEx_v10 liver, eQTL Catalogue R8 `QTD000266`, `quant_method=ge`,
fixed for all pairs.

**New rule (R3, amended 2026-08-26):** **eQTLGen phase I cis-eQTL, whole blood, n ≈ 31,684.**
Fixed for all pairs. No per-pair selection. No maximum over resources or tissues.

**Reason:** instrument power. The old resource produced zero instruments at p < 5e-8 for all three
positive controls (Deliverable 1). The change is to the resource only; the threshold is frozen by
Decision 6.

**Sequencing of the change is the point:** the rule was amended and recorded *before* any re-run of
the controls. Deliverable 1b (instrument yield census) is a lookup against the new resource and is
explicitly **not** a re-run.

**Known consequence, accepted in advance:** eQTLGen is whole blood only. Targets not expressed in
blood will be systematically penalised. Deliverable 1b quantifies that penalty by therapeutic area
before any further pipeline work.

**Known open issue:** eQTLGen phase I is distributed on **GRCh37**, while the GWAS Catalog
harmonised files and the Deliverable 1 chain are **GRCh38**. A build-reconciliation rule must be
pre-registered before the controls are re-run. Not decided here.

---

## DECISION 8 — PCSK9 retired as a control, NOT replaced — FROZEN 2026-08-27

**Status: ABSENCE, not failure.** PCSK9 (`ENSG00000169174`) appears in **zero of 127,341,798 rows**
of the full eQTLGen phase I cis-eQTL file, verified by direct scan of the raw file rather than
inferred from the tested-gene list. PCSK4, PCSK5, PCSK6 and PCSK7 are all present, so this is not a
nomenclature artefact. PCSK9 is a hepatocyte-secreted protein and is not appreciably expressed in
whole blood.

**No substitute control will be chosen.** Selecting a replacement control after seeing which genes
eQTLGen covers would be selection on a known result. The surviving controls are **HMGCR and IL6R
only**, and the reduced control set is itself reported as a limitation.

---

## DECISION 9 — instrument availability is a PRIMARY RESULT — FROZEN 2026-08-27

Instrument availability is not a data-cleaning statistic to be mentioned in a methods footnote. It
is a finding about which drug targets the causal-genetics toolkit can and cannot see, and it is
accumulated deliberately in `INSTRUMENT_AVAILABILITY.md`.

Every target carries one of three states: **absent from the exposure resource**, **tested but
null**, **instrumented**. These are never merged. Where possible each absence is annotated with
whether it tracks tissue-restricted expression or secreted-protein status.

**PCSK9 is the anchor example**: liver-expressed, secreted, the target of one of the most
successful cardiovascular drug classes of the last two decades — and completely invisible to a
whole-blood exposure resource. A pipeline that silently drops it would report a clean result on a
target list that had been truncated in a way correlated with biology.

---

## DECISION 10 — Lane A and Lane B use DIFFERENT exposure resources — FROZEN 2026-08-27

Declared before any cardiology screening has been run.

**Lane A (benchmark, thousands of indications):** eQTLGen phase I whole blood, n ≈ 31,684. One
resource, one uniform rule across all targets. Measured yield: 1,576 / 2,517 targets (62.6 %).

**Lane B (cardiology discovery):** a pre-registered cardiac/hepatic/vascular tissue set, **not**
whole blood — specified in full in `deliverables/d1/PREREGISTRATION.md` §R3b. Whole blood cannot
see cardiac- and hepatic-restricted targets, so a blood-only cardiology screen would return a
truncated target list that does not look truncated.

**Honest power verdict, measured not assumed (see §R3b):** the cardiac/vascular expression
resources are **underpowered**. They instrument 8.8–10.6 % of Minikel targets individually and
16.8 % as a union of the three cardiac/aortic tissues, versus 62.6 % for whole blood — and those
are *upper bounds*. Lane B at expression level will be instrument-limited, and protein-level
exposures are the technically indicated route. That is a finding, not a failure.

**OPEN — CREDENTIAL GATE, requires the user's decision:** the natural protein-level exposure
(UKB-PPP, n ≈ 54,000) is **not** in the GWAS Catalog (verified: PMID 37794186 returns 0 studies)
and its Synapse project `syn51364943` requires credentials. Per the standing rule on credentials I
have stopped rather than worked around it. See §R3b for the alternatives and what remains
UNVERIFIED.

---

## DECISION 11 — IL6R control: FAIL stands, reclassified as a FINDING — FROZEN 2026-08-30

**Verdict unchanged: the IL6R → CRP control FAILED against its pre-registered expectation.** That
stays on the record. What changes is the *interpretation*, and the change is made by the user, not
by the analyst who ran it.

**Mechanism, now understood:** eQTLGen measures **total IL6R transcript**. Higher transcript → more
receptor → more IL-6 signalling → higher CRP. **POSITIVE is mechanistically correct for that
exposure.** The canonical negative result derives from rs2228145 (Asp358Ala), a **shedding** variant
acting at the splicing and protein level, which a total-mRNA exposure structurally cannot
represent — and which was verified absent from eQTLGen's significant cis-eQTL set for IL6R
(rs4129267 present only at p = 1.11e-05, above the frozen 5e-8 threshold).

**Recorded as:** control FAILED against its pre-registered expectation; the expectation was
**mis-specified for the exposure modality**; **the pipeline is not implicated.**

Note on process: the analyst declined to convert the FAIL into a PASS using this argument, because
the argument was constructed after seeing the result. The reclassification is the user's override
of interpretation, with the verdict left intact. The distinction is the point.

---

## DECISION 12 — two failure classes, and Class B is the more dangerous — FROZEN 2026-08-30

Defined in `INSTRUMENT_AVAILABILITY.md`; every target is classified against them from here on.

* **CLASS A — target absent or uninstrumented in the exposure resource.** Anchor: **PCSK9**
  (0 of 127,341,798 eQTLGen rows). The analysis visibly cannot run.
* **CLASS B — target well-instrumented, but the exposure modality cannot represent the mechanism
  the drug acts on** (splicing, shedding, protein activity, localisation). Anchor: **IL6R**
  (F = 134.5, p = 1.16e-59 — a textbook-strong instrument, and the answer is still uninterpretable
  for the drug mechanism).

**Class B is strictly worse than Class A, and is invisible on every existing platform, because the
analysis appears to succeed.** Class A announces itself as a missing row. Class B produces a
confident, well-powered, correctly computed estimate of the wrong quantity.

---

## DECISION 13 — named positive controls are retired entirely — FROZEN 2026-08-30

No replacement controls will be proposed. Enough is now known about eQTLGen's coverage that any
newly chosen named control would be selected on information already in hand, and would validate
nothing.

Validation becomes **rule-based and distributional**: concordance between cis-MR direction and the
direction implied by approved-drug mechanism, across every qualifying target-indication pair, with
the selection rule fixed before computation. Specified and executed in Part C.

---

## GATE D (colocalization) — MARKED NOT TRUSTWORTHY, project-wide — 2026-08-30

Colocalization returned PP3 ≈ 1.000 / PP4 ≈ 1e-11 at HMGCR–LDL-C, one of the most robust
colocalizations in human genetics. No colocalization output from this project may be used, cited or
carried into any downstream gate until the cause is identified and fixed. Status tracked in
`PROTOCOL.md`.
