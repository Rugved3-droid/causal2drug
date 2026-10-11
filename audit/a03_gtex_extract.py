"""AUDIT: pull GTEx v8 rows for the primary-set targets out of the tar, unchanged.

Keeps a signif_variant_gene_pairs row when the gene's Ensembl ID (version stripped) is a target's
GRCh37 ID or the ID HGNC currently records for the target symbol (needed for RECOVERY_RULE
amendment 8), or - for reporting only - when nothing else. No p or distance filter here.

Run:  python audit/a03_gtex_extract.py
Writes audit/out/gtex_egenes_all.tsv   (every gene in every tissue's egenes file: coordinates)
       audit/out/gtex_signif_targets.tsv
"""
import csv, gzip, io, os, tarfile

ROOT = r"D:\causal2drug"
TAR = os.path.join(ROOT, "manuscript_prep", "manuscript", "recovery_data", "GTEx_Analysis_v8_eQTL.tar")
OUT = os.path.join(ROOT, "audit", "out")

ids = set()
for r in csv.DictReader(open(os.path.join(OUT, "gene_coords_audit.tsv"), encoding="utf-8"), delimiter="\t"):
    for k in ("ensembl_id", "hgnc_current_ensembl_id"):
        if r[k]:
            ids.add(r[k])
print("target ids (GRCh37 + HGNC current):", len(ids))

eg = open(os.path.join(OUT, "gtex_egenes_all.tsv"), "w", newline="", encoding="utf-8")
sg = open(os.path.join(OUT, "gtex_signif_targets.tsv"), "w", newline="", encoding="utf-8")
eg.write("tissue\tgene_id\tgene_name\tgene_chr\tgene_start\tgene_end\tstrand\n")
sg.write("tissue\tvariant_id\tgene_id\tpval_nominal\tslope\tslope_se\n")
n_e = n_s = 0
with tarfile.open(TAR) as tf:
    for m in tf:
        base = os.path.basename(m.name)
        if base.endswith(".v8.egenes.txt.gz"):
            tissue = base[:-len(".v8.egenes.txt.gz")]
            with gzip.open(tf.extractfile(m), "rt") as fh:
                h = fh.readline().rstrip("\n").split("\t")
                ix = [h.index(c) for c in ("gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand")]
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    eg.write(tissue + "\t" + "\t".join(f[i] for i in ix) + "\n"); n_e += 1
        elif base.endswith(".v8.signif_variant_gene_pairs.txt.gz"):
            tissue = base[:-len(".v8.signif_variant_gene_pairs.txt.gz")]
            with gzip.open(tf.extractfile(m), "rt") as fh:
                h = fh.readline().rstrip("\n").split("\t")
                iv, ig, ip, ib, ise = (h.index(c) for c in ("variant_id", "gene_id", "pval_nominal", "slope", "slope_se"))
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[ig].split(".")[0] in ids:
                        sg.write("\t".join([tissue, f[iv], f[ig], f[ip], f[ib], f[ise]]) + "\n"); n_s += 1
            print(tissue, "done; signif rows kept so far", n_s, flush=True)
print("egenes rows", n_e, "signif rows", n_s)
