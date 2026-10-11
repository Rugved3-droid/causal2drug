#!/usr/bin/env bash
# AUDIT: one pass over the raw eQTLGen full cis-eQTL file (127,341,798 rows, ~4 min).
# Needs audit/out/gene_coords_audit.tsv from a01. Run from D:\causal2drug in Git Bash:
#   bash audit/a00_eqtlgen_scan.sh
# Writes audit/out/eqtlgen_tested_universe_audit.tsv  (every distinct Gene ID + its row count)
#        audit/out/eqtlgen_target_rows.tsv            (every row, any p, for a primary-set target's Ensembl ID)
#        audit/out/eqtlgen_nrows.txt
set -euo pipefail
cut -f3 audit/out/gene_coords_audit.tsv | tail -n +2 | grep ENSG > audit/out/target_ensg.txt
gzip -dc data/raw/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | gawk -F'\t' -v OUTD=audit/out '
  BEGIN { while ((getline l < (OUTD"/target_ensg.txt")) > 0) T[l] = 1 }
  NR == 1 { print > (OUTD"/eqtlgen_target_rows.tsv"); next }
  { n[$8]++; if ($8 in T) print > (OUTD"/eqtlgen_target_rows.tsv") }
  END { for (g in n) print g"\t"n[g] > (OUTD"/eqtlgen_tested_universe_audit.tsv"); print NR-1 > (OUTD"/eqtlgen_nrows.txt") }'
