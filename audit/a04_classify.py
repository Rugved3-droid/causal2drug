"""AUDIT step 2b: eQTLGen status of every primary-set target, from DEFINITIONS.md Part 1.

  unresolved gene record                                  -> UNRESOLVED
  Ensembl ID in no row of the full eQTLGen cis file       -> NO_EQTL_DATA  (autosomal / X / MT by GTF chromosome)
  >=1 row p < 5e-8, same chromosome, gene body +-100 kb   -> INSTRUMENTED
  otherwise                                               -> NO_QUALIFYING_INSTRUMENT

Inputs are my own scan of the raw eQTLGen file (see AUDIT.md for the gawk command).
Run:  python audit/a04_classify.py
Writes audit/out/audit_per_target.csv
"""
import csv, collections, os

OUT = r"D:\causal2drug\audit\out"
FLANK, P = 100_000, 5e-8

tested = {l.split("\t")[0] for l in open(os.path.join(OUT, "eqtlgen_tested_universe_audit.tsv"))}
nrows = open(os.path.join(OUT, "eqtlgen_nrows.txt")).read().strip()
print(f"[B01] eQTLGen full file: {nrows} data rows, {len(tested)} distinct genes tested")

coords = {r["symbol"]: r for r in csv.DictReader(open(os.path.join(OUT, "gene_coords_audit.tsv"), encoding="utf-8"), delimiter="\t")}
by_id = {r["ensembl_id"]: r for r in coords.values() if r["ensembl_id"]}

nhit = collections.Counter(); best = {}; chroms = collections.Counter()
with open(os.path.join(OUT, "eqtlgen_target_rows.tsv"), encoding="utf-8") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    iP, iC, iPos, iG, iS, iZ = (h.index(c) for c in ("Pvalue", "SNPChr", "SNPPos", "Gene", "SNP", "Zscore"))
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c = by_id[f[iG]]
        chroms[f[iC]] += 1
        if float(f[iP]) >= P or f[iC] != c["chrom"]:
            continue
        if int(c["start"]) - FLANK <= int(f[iPos]) <= int(c["end"]) + FLANK:
            nhit[f[iG]] += 1
            z = abs(float(f[iZ]))
            if f[iG] not in best or z > best[f[iG]][0]:
                best[f[iG]] = (z, f[iS], f[iP])
print(f"[B02] SNP chromosomes present in target rows: {sorted(chroms, key=lambda x: (len(x), x))}")
print(f"[B02] chromosomes among ALL tested genes cannot be read from the universe file; "
      f"X/MT targets tested: {sum(1 for r in coords.values() if r['chrom'] in ('X','MT','Y') and r['ensembl_id'] in tested)}")

rows = []
for s, r in sorted(coords.items()):
    eid = r["ensembl_id"]
    if not eid:
        st, sub = "UNRESOLVED", ""
    elif eid not in tested:
        st = "NO_EQTL_DATA"
        sub = {"X": "x_linked", "MT": "mitochondrial", "Y": "y_linked"}.get(r["chrom"], "autosomal")
    elif nhit[eid]:
        st, sub = "INSTRUMENTED", ""
    else:
        st, sub = "NO_QUALIFYING_INSTRUMENT", ""
    rows.append(dict(symbol=s, phase=r["phase"], ensembl_id=eid, chrom=r["chrom"], start=r["start"], end=r["end"],
                     status=st, subtype=sub, n_qualifying=nhit[eid] if eid else "",
                     lead_snp=best.get(eid, ("", "", ""))[1], lead_p=best.get(eid, ("", "", ""))[2], how=r["how"]))
with open(os.path.join(OUT, "audit_per_target.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

c = collections.Counter(r["status"] for r in rows); sub = collections.Counter(r["subtype"] for r in rows if r["status"] == "NO_EQTL_DATA")
n = len(rows)
print(f"[B03] targets {n}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
print(f"[B04] NO_EQTL_DATA {c['NO_EQTL_DATA']} = autosomal {sub['autosomal']} + X-linked {sub['x_linked']} + mitochondrial {sub['mitochondrial']}")
print(f"[B05] NO_QUALIFYING_INSTRUMENT {c['NO_QUALIFYING_INSTRUMENT']}; upper bound (sum) {c['NO_EQTL_DATA']+c['NO_QUALIFYING_INSTRUMENT']}")
print(f"[B06] lower bound {c['NO_EQTL_DATA']}/{n} = {100*c['NO_EQTL_DATA']/n:.1f}% ; upper {100*(c['NO_EQTL_DATA']+c['NO_QUALIFYING_INSTRUMENT'])/n:.1f}%")
nonsex = sum(1 for r in rows if r["chrom"] not in ("X", "Y", "MT") and r["status"] != "UNRESOLVED")
nonsex_all = sum(1 for r in rows if r["chrom"] not in ("X", "Y", "MT"))
print(f"[B07] autosomal-only lower bound: {sub['autosomal']}/{nonsex_all} (denominator incl. 3 unresolved) = {100*sub['autosomal']/nonsex_all:.1f}% ; "
      f"/{nonsex} (excl. unresolved) = {100*sub['autosomal']/nonsex:.1f}%")
la = [r for r in rows if r["phase"] == "5"]
print(f"[B08] launched {len(la)}: " + ", ".join(f"{k} {v}" for k, v in sorted(collections.Counter(r['status'] for r in la).items())))
