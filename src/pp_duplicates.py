#!/usr/bin/env python
"""Duplicate analysis of Minikel pp.tsv, recomputed from the raw file every run.

Reports, in order:
  A. exact duplicate full lines
  B. duplicated ti_uid keys, every row printed in full
  C. duplicated (gene, indication_mesh_id) pairs
  D. key cardinality
  E. classification of duplicated keys into identical / outcome-conflicting /
     annotation-only-conflicting

Usage:  python src/pp_duplicates.py [path-to-pp.tsv]
"""
import collections
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 else r"D:\causal2drug\data\raw\pp.tsv"

OUTCOME_COLS = ["hcat", "acat", "ccat", "hcatnum", "acatnum", "ccatnum",
                "succ_p_1", "succ_1_2", "succ_2_3", "succ_3_a", "year_launch"]

raw = open(PATH, "rb").read()
print("file:", PATH)
print("bytes:", len(raw))
print("ends_with_newline:", raw.endswith(b"\n"))
print("CR_bytes_present:", b"\r" in raw)

lines = raw.decode("utf-8").split("\n")
if lines and lines[-1] == "":
    lines = lines[:-1]
header = lines[0].split("\t")
data = lines[1:]
print("total_lines_incl_header:", len(lines))
print("n_data_rows:", len(data))
print("n_cols_header:", len(header))
print("field_count_histogram:",
      dict(collections.Counter(len(l.split("\t")) for l in data)))

i_uid = header.index("ti_uid")
i_gene = header.index("gene")
i_mesh = header.index("indication_mesh_id")

print("\n=== A. EXACT DUPLICATE FULL LINES ===")
c_full = collections.Counter(data)
dup_full = {k: v for k, v in c_full.items() if v > 1}
print("n_distinct_lines_duplicated:", len(dup_full))
print("n_excess_rows:", sum(v - 1 for v in dup_full.values()))
for k, v in sorted(dup_full.items()):
    print("count=%d | %s" % (v, k))

print("\n=== B. DUPLICATE ti_uid ===")
by_uid = collections.defaultdict(list)
for l in data:
    by_uid[l.split("\t")[i_uid]].append(l)
dup_uid = {k: v for k, v in by_uid.items() if len(v) > 1}
print("n_distinct_ti_uid:", len(by_uid))
print("n_ti_uid_duplicated:", len(dup_uid))
for k in sorted(dup_uid):
    print("--- ti_uid=%s count=%d ---" % (k, len(dup_uid[k])))
    for l in dup_uid[k]:
        print("   " + l)

print("\n=== C. DUPLICATE (gene, indication_mesh_id) PAIRS ===")
by_gi = collections.defaultdict(list)
for l in data:
    f = l.split("\t")
    by_gi[(f[i_gene], f[i_mesh])].append(l)
dup_gi = {k: v for k, v in by_gi.items() if len(v) > 1}
print("n_distinct_gene_indication_pairs:", len(by_gi))
print("n_pairs_duplicated:", len(dup_gi))
for k in sorted(dup_gi):
    print("--- %s / %s count=%d ---" % (k[0], k[1], len(dup_gi[k])))
    for l in dup_gi[k]:
        print("   " + l)

print("\n=== D. KEY CARDINALITY ===")
print("distinct_genes:", len({l.split("\t")[i_gene] for l in data}))
print("distinct_indication_mesh_id:", len({l.split("\t")[i_mesh] for l in data}))

print("\n=== E. CLASSIFICATION OF DUPLICATED ti_uid ===")
identical, conflict_outcome, conflict_annot = [], [], []
for k in sorted(dup_uid):
    rows = [l.split("\t") for l in dup_uid[k]]
    if all(r == rows[0] for r in rows):
        identical.append(k)
        continue
    diff = [c for c in header if len({r[header.index(c)] for r in rows}) > 1]
    (conflict_outcome if any(c in OUTCOME_COLS for c in diff)
     else conflict_annot).append((k, diff))
print("identical rows                       :", len(identical), identical)
print("CONFLICT in outcome-bearing columns  :", len(conflict_outcome))
for k, d in conflict_outcome:
    print("     ", k, "-> differs in", d)
print("conflict only in oto_/orphan columns :", len(conflict_annot))
for k, d in conflict_annot:
    print("     ", k, "-> differs in", d)

print("\nrows if dedup by identical-line :", len(set(data)))
print("rows if dedup by ti_uid         :", len(by_uid))
