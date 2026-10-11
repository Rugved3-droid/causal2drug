"""AUDIT: dump the pQTL workbook sheets named in RECOVERY_RULE.md to TSV, unchanged.

Run:  python audit/a02_dump_pqtl_sheets.py
Writes audit/out/pqtl/<study>_<sheet>.tsv (every row, every column, as read by openpyxl)
"""
import csv, os, openpyxl

D = r"D:\causal2drug\manuscript_prep\manuscript\recovery_data"
OUT = r"D:\causal2drug\audit\out\pqtl"
os.makedirs(OUT, exist_ok=True)
JOBS = [("ukbppp", "41586_2023_6592_MOESM3_ESM.xlsx", ["ST3", "ST9", "ST10"]),
        ("decode", "ferkingstad2021_MOESM4.xlsx", ["ST01", "ST02"]),
        ("fenland", "pietzner2021_supplementary_tables.xlsx", ["Header_ST2", "ST2", "Header_ST3", "ST3"])]
for study, fn, sheets in JOBS:
    wb = openpyxl.load_workbook(os.path.join(D, fn), read_only=True, data_only=True)
    for s in sheets:
        n = 0
        with open(os.path.join(OUT, f"{study}_{s}.tsv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            for row in wb[s].iter_rows(values_only=True):
                w.writerow(["" if v is None else v for v in row]); n += 1
        print(study, s, n, "rows", flush=True)
