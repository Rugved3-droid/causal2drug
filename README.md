# causal2drug — cardiovascular drug targets invisible to eQTLGen whole-blood expression instruments

Code, rule documents and derived data for the study "Cardiovascular Drug Targets Invisible to eQTLGen Whole-Blood Expression Instruments: A Target-Level Map of Where Candidate Instruments Exist" (Parmar R, Fahim M).

The study asks, for 726 targets with a cardiovascular drug programme, whether eQTLGen whole
blood supplies a cis-eQTL instrument, and, for the 280 that it does not, whether a candidate
instrument exists in another GTEx v8 tissue or in a published plasma cis-pQTL table. The gap is
specific to eQTLGen: 6 of 10 lipid targets with landmark human genetic evidence have no eQTLGen
instrument (PCSK9, APOC3, ANGPTL3, APOB, NPC1L1, LPA), while GTEx whole blood holds one for PCSK9.

## Where to start

| file | what it is |
|---|---|
| `manuscript_prep/KEY_NUMBERS.md` | the headline numbers and the command that produces each |
| `manuscript_prep/manuscript/NUMBER_TRACE.md` | every number in the paper, traced to an output file and line |
| `manuscript_prep/manuscript/recovery_map.csv` | the recovery map, one row per target (280) |
| `manuscript_prep/analysis/out/per_target_primary.csv` | status of all 726 targets |
| `manuscript_prep/manuscript/RECOVERY_RULE.md` | the recovery rule, with every clarification, amendment and correction, dated |
| `manuscript_prep/manuscript/RESCUE_CONTROLS_PROTOCOL.md`, `RESCUE_CONTROLS_RESULTS.md` | positive controls with candidate instruments |
| `manuscript_prep/manuscript/IL6R_AUDIT.md` | audit of the IL6R control |
| `manuscript_prep/manuscript/APPROVED_35_CHECK.csv` | FDA and EMA approval records for the approved-drug targets |
| `manuscript_prep/drug_links/` | ChEMBL drug–target–indication links and the completed review |
| `manuscript_prep/manuscript/display/`, `manuscript_prep/tables/`, `manuscript_prep/figures/` | tables and figures |
| `audit/AUDIT.md` | a second implementation, written in a separate AI-assisted session from the rule documents alone, without access to the original code or per-target outputs, and what it found |
| `PROTOCOL.md`, `DECISIONS.md`, `manuscript_prep/DECISIONS_LOG.md` | protocol and decision history, including post hoc decisions |
| `archive/` | earlier deliverables and superseded notes, kept for the record; not needed to follow the paper |

## What is not here

Third-party data are not redistributed. To rerun the pipeline, obtain them from their sources
(versions, file names and checksums are in `manuscript_prep/METHODS_FACTS.md` and
`manuscript_prep/manuscript/DATA_ACCESS.md`):

- Minikel et al. target–indication table (`ericminikel/genetic_support`, tag v1.0)
- eQTLGen phase I cis-eQTL files
- Ensembl GRCh37 release 87 annotation
- GTEx v8 single-tissue cis-eQTL archive; GTEx v10 median expression
- UKB-PPP (Sun 2023), deCODE (Ferkingstad 2021) and Fenland (Pietzner 2021) supplementary tables
- GWAS Catalog harmonised summary statistics; 1000 Genomes phase 3
- ChEMBL 37, Human Protein Atlas, HGNC

Scripts read these by absolute path under `D:\causal2drug`; adjust the paths to run elsewhere.
The manuscript text is not in this repository.

## Reproduce

Order of commands: `manuscript_prep/manuscript/README.md`.

## Honest limits

Validation of the Mendelian randomization pipeline is thin: of eight positive controls, three
passed, two failed and three could not be run. Several rules and the reporting tiers were set
after results were seen; each is dated and labelled in the rule files. The protocol documents
were not registered or externally timestamped. The recovery map records candidate instruments,
not validated ones.

## A note on names

The recovery tiers were first called "plausible rescue" and "off-target-tissue rescue". They
were renamed "candidate: CV-relevant tissue or plasma protein" and "candidate: other tissue
only", because the map records candidate instruments and not validated ones. The rule file,
and the decision log were written before the renaming and keep the first names; the
definitions did not change.

## Changes to the plan, on the record

The first protocol (`archive/deliverables/d1/PREREGISTRATION.md`) fixed GTEx liver as the
source of instruments for the positive controls. No control had an instrument there, and the
exposure resource was changed to eQTLGen whole blood in a dated amendment in the same file
("AMENDMENT R3a"), recorded before the controls were re-run but with the liver result known.
The significance threshold was frozen at 5e-8 after one control was seen to pass at a laxer
value (`DECISIONS.md`, decision 6). An earlier colocalization arm was tested against a
criterion written in advance, refuted, and dropped (`manuscript_prep/DEVIATION_CLASS_C_FIGURE.md`).

## Sources and attribution

Derived files contain values from the resources listed above, each of which has its own terms
of use; ChEMBL data are provided under CC BY-SA 3.0. Cite the original resources when using
anything derived from them. The STROBE-MR checklist wording is copyright of the Equator Network
(CC BY 3.0). Code and project documents are under the MIT licence (`LICENSE`).

## Contact

Rugved Parmar, parmarrugved@gmail.com
