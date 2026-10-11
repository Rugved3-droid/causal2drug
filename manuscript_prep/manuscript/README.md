# manuscript/ — draft for Circulation: Genomic and Precision Medicine

Drafted and revised 2026-10-11. Start with `OPEN_ITEMS.md`, then `RESCUE_CONTROLS_RESULTS.md` and `IL6R_AUDIT.md`.

| file | what it is |
|---|---|
| `manuscript.docx`, `manuscript.md` | the draft: title page, structured abstract, text, references, legends, Tables 1–3, Figures 1–3 |
| `supplement.docx`, `supplement.md` | Supplemental Methods, Tables S1–S23 (long ones in part), Figures S1–S2, STROBE-MR checklist |
| `cover_letter.md`, `cover_letter.docx` | cover letter |
| `WORD_COUNT.md` | counts against the journal limits |
| `NUMBER_TRACE.md` | every number in the manuscript → output file and line → command |
| `OPEN_ITEMS.md` | decisions, unverified items, known weaknesses, title options |
| `RESCUE_CONTROLS_PROTOCOL.md`, `RESCUE_CONTROLS_RESULTS.md`, `rescue_controls/` | five positive controls using candidate instruments: protocol (written first), results, script and raw rows |
| `IL6R_AUDIT.md`, `il6r_audit/` | audit of the IL6R control: variant identity, harmonization, LD with rs2228145 |
| `revision_sources/` | GTEx v8 tissue sample sizes, official STROBE-MR checklist, FDA and EMA drug records, the Zenodo preprint text |
| `RECOVERY_RULE.md` | the pre-specified recovery rule and its dated amendments (8 and 9 are post hoc) |
| `DATA_ACCESS.md` | which resources were obtained, and how |
| `recovery_map.csv`, `recovery_map_summary.txt` | the recovery map (280 targets) and its summary counts, tags [R01]–[R31] |
| `display/` | Table 2, Table 3, Tables S11–S21, STROBE-MR checklist, `recovery_tables.xlsx`, Figure 3 (PDF, SVG, PNG) |
| `recovery_out/` | per-tissue and per-study detail, ID-change evidence, the first-run summary kept as a record |
| `recovery_scripts/` | `build_recovery_map.py`, `make_recovery_outputs.py`, `fetch_id_change_records.py`, `inspect_xlsx.py` |
| `draft/` | `manuscript_src.md`, `supplement_src.md`, `build_manuscript.py`, `fetch_refs.py`, `refs_crossref.json`, `make_number_trace.py` |
| `PRIOR_ART.md`, `prior_art_search/` | the prior-art search and its verdict |
| `recovery_data/` | downloaded resources (1.6 GB; not in the zip) |

Figures 1 and 2 and Tables 1 and S1–S10 are in `../figures/` and `../tables/`.

Rebuild, from `D:\causal2drug`:

```
python manuscript_prep\analysis\classify.py > manuscript_prep\analysis\classify_output.txt
python manuscript_prep\tables\make_tables.py
python manuscript_prep\figures\make_figures.py
python -I manuscript_prep\manuscript\recovery_scripts\build_recovery_map.py > manuscript_prep\manuscript\recovery_map_summary.txt
python manuscript_prep\manuscript\rescue_controls\run_rescue_controls.py
python manuscript_prep\manuscript\il6r_audit\il6r_audit.py > manuscript_prep\manuscript\il6r_audit\il6r_audit_output.txt
python manuscript_prep\analysis\phase_ci.py > manuscript_prep\analysis\phase_ci_output.txt
python manuscript_prep\manuscript\revision_sources\fetch_drug_status.py
python manuscript_prep\manuscript\recovery_scripts\make_recovery_outputs.py
python manuscript_prep\manuscript\draft\make_revision_tables.py
python manuscript_prep\manuscript\draft\build_manuscript.py
python manuscript_prep\manuscript\draft\make_number_trace.py
```
