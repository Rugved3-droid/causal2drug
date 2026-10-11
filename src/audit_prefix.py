#!/usr/bin/env python
"""Part B3: audit which D0 / D1 / D1b numbers were computed BEFORE the tabix
upper-bound fix, and recompute every pre-fix number that used tabix.

The defect: the tabix linear index was used as an upper bound as well as a lower bound.
It is valid only as a lower bound, so regions could be silently truncated.
Numbers that never passed through tabix (gzip streaming, HEAD probes, GTF parsing,
checksums) are unaffected by construction and are listed as such rather than re-run.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T  # noqa: E402

GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
EQ = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
      "QTS000015/QTD000266/QTD000266.all.tsv.gz")
KG38 = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
        "1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

CASES = [
    dict(label="D0 GWAS GCST90239658 chr5 HMGCR window", prev=3323,
         url=GW + "GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         chrom="5", start=75236329, end=75464241),
    dict(label="D0 eQTL QTD000266 chr5 HMGCR window (all genes)", prev=8990,
         url=EQ, chrom="5", start=75236329, end=75464241, gene=None),
    dict(label="D0 1000G GRCh38 chr5 HMGCR window", prev=5644,
         url=KG38, chrom="chr5", start=75236329, end=75464241),
    dict(label="D1 HMGCR GWAS GCST90239658", prev=3323,
         url=GW + "GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         chrom="5", start=75236329, end=75464241),
    dict(label="D1 HMGCR eQTL (ENSG00000113161 only)", prev=592,
         url=EQ, chrom="5", start=75236329, end=75464241, gene="ENSG00000113161"),
    dict(label="D1 PCSK9 GWAS GCST90239658", prev=3925,
         url=GW + "GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         chrom="1", start=54939445, end=55164852),
    dict(label="D1 PCSK9 eQTL (ENSG00000169174 only)", prev=846,
         url=EQ, chrom="1", start=54939445, end=55164852, gene="ENSG00000169174"),
    dict(label="D1 IL6R GWAS GCST90179146", prev=581,
         url=GW + "GCST90179001-GCST90180000/GCST90179146/harmonised/GCST90179146.h.tsv.gz",
         chrom="1", start=154305141, end=154569450),
    dict(label="D1 IL6R eQTL (ENSG00000160712 only)", prev=839,
         url=EQ, chrom="1", start=154305141, end=154569450, gene="ENSG00000160712"),
]


def main():
    out = []
    print("%-48s %10s %10s %8s %6s" % ("query", "before", "after", "match", "spans"))
    print("-" * 88)
    for c in CASES:
        rows, st = T.query(c["url"], c["chrom"], c["start"], c["end"])
        if c.get("gene"):
            rows = [r for r in rows if len(r) > 16 and r[16] == c["gene"]]
        n = len(rows)
        ok = n == c["prev"]
        print("%-48s %10d %10d %8s %6s"
              % (c["label"][:48], c["prev"], n, "YES" if ok else "*** NO ***",
                 st.get("n_spans")))
        out.append(dict(label=c["label"], before=c["prev"], after=n, match=ok,
                        n_spans=st.get("n_spans"), bytes=st.get("bytes")))
    print("-" * 88)
    n_ok = sum(1 for r in out if r["match"])
    print("%d/%d pre-fix numbers reproduce unchanged" % (n_ok, len(out)))
    json.dump(out, open(r"D:\causal2drug\deliverables\d1c\prefix_audit.json", "w"),
              indent=2)
    return 0 if n_ok == len(out) else 1


if __name__ == "__main__":
    os.makedirs(r"D:\causal2drug\deliverables\d1c", exist_ok=True)
    sys.exit(main())
