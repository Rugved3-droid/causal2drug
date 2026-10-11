#!/usr/bin/env python
"""Recovery map: for every primary-set target with no instrument in eQTLGen whole blood, is
there a cis instrument in GTEx v8 (any tissue) or in a published plasma cis-pQTL table?

Implements manuscript_prep/manuscript/RECOVERY_RULE.md, including its 2026-10-11 amendments.
Instrument = variant with p < 5e-8 within gene body +-100,000 bp, in the build of the resource.

    python -I manuscript_prep/manuscript/recovery_scripts/build_recovery_map.py \
        > manuscript_prep/manuscript/recovery_map_summary.txt

Reads   manuscript_prep/analysis/out/per_target_primary.csv
        manuscript_prep/manuscript/recovery_data/{GTEx tar, UKB-PPP, deCODE, Fenland workbooks}
        data/raw/Homo_sapiens.GRCh37.87.gtf.gz   (only to establish the Fenland build)
Writes  manuscript_prep/manuscript/recovery_map.csv
        manuscript_prep/manuscript/recovery_out/*.csv|tsv
"""
import collections
import csv
import gzip
import json
import os
import re
import sys
import tarfile

import openpyxl

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
MP = os.path.dirname(MS)
ROOT = os.path.dirname(MP)
DATA = os.path.join(MS, "recovery_data")
OUT = os.path.join(MS, "recovery_out")
os.makedirs(OUT, exist_ok=True)

P_MAX = 5e-8
LOG10_MIN = 7.30102999566398           # -log10(5e-8)
WINDOW = 100000
CV_TISSUES = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Aorta",
              "Artery_Coronary", "Artery_Tibial", "Liver", "Kidney_Cortex", "Adrenal_Gland"]
BLOOD = "Whole_Blood"
NAMED = ["AGTR1", "ANGPTL3", "APOB", "APOC3", "AVPR2", "F2", "PCSK9", "PLG", "REN", "SLC5A1",
         "IL1B", "ITGA2B", "NPC1L1"]
PHASES = ["Preclinical", "Phase I", "Phase II", "Phase III", "Launched"]
CATS = ["fixable-tissue", "fixable-protein", "fixable-both", "locked"]


def clean(v):
    return "" if v is None else str(v).strip().lstrip("'").strip()


def chrom(v):
    c = clean(v).replace("chr", "")
    return {"23": "X", "24": "Y"}.get(c, c)


def dist(pos, a, b):
    return 0 if a <= pos <= b else min(abs(pos - a), abs(pos - b))


def tokens(v):
    return [t for t in re.split(r"[|;,/\s]+", clean(v)) if t]


# ------------------------------------------------------------------ scope
rows = list(csv.DictReader(open(os.path.join(MP, "analysis", "out", "per_target_primary.csv"),
                                encoding="utf-8")))
scope = [r for r in rows if r["status"] in ("NO_EQTL_DATA", "NO_QUALIFYING_INSTRUMENT")]
unres = sorted(r["gene"] for r in rows if r["status"] == "UNRESOLVED")
by_ens = {r["ensembl_id"]: r for r in scope}
by_sym = {r["gene"]: r for r in scope}
assert len(by_ens) == len(scope) == len(by_sym)
print("[R01] primary-set targets: %d; in scope (no eQTL data or no qualifying instrument): %d;"
      " out of scope, unresolved gene record: %d (%s)" % (len(rows), len(scope), len(unres),
                                                         ", ".join(unres)))
assert all(g in by_sym for g in NAMED), [g for g in NAMED if g not in by_sym]

# ------------------------------------------------------------------ documented ID changes
# Amendment 8 (2026-10-11). HGNC's current Ensembl gene ID for the target symbol, where it
# differs from the GRCh37 ID used in the main analysis. Cached by fetch_id_change_records.py
# for the targets that were in no GTEx file. The alternative ID is used ONLY if the GRCh37 ID
# is in no GTEx egenes file and the checks printed under [R02b] pass.
IDC = os.path.join(OUT, "id_change")
hgnc_now = {}
for r in scope:
    f = os.path.join(IDC, "hgnc_%s.json" % r["gene"])
    if os.path.exists(f):
        d = json.load(open(f, encoding="utf-8"))["response"]["docs"]
        if d and d[0].get("ensembl_gene_id") and d[0]["ensembl_gene_id"] != r["ensembl_id"]:
            hgnc_now[r["gene"]] = (d[0]["hgnc_id"], d[0]["ensembl_gene_id"])
alt_to_gene = {e: g for g, (_, e) in hgnc_now.items()}
WANT = set(by_ens) | set(alt_to_gene)

# ------------------------------------------------------------------ (a) GTEx v8
EG = os.path.join(OUT, "gtex_target_egenes.tsv")
PR = os.path.join(OUT, "gtex_target_signif_pairs.tsv")
V26 = os.path.join(OUT, "gencode_v26_gene_bodies_from_gtex_egenes.tsv")
if not (os.path.exists(EG) and os.path.exists(PR) and os.path.exists(V26)):
    v26 = {}
    with tarfile.open(os.path.join(DATA, "GTEx_Analysis_v8_eQTL.tar")) as tf, \
            open(EG, "w", newline="", encoding="utf-8") as feg, \
            open(PR, "w", newline="", encoding="utf-8") as fpr:
        weg, wpr = csv.writer(feg, delimiter="\t"), csv.writer(fpr, delimiter="\t")
        weg.writerow(["tissue", "ensembl_id", "gene_id", "gene_name", "chr", "gene_start",
                      "gene_end", "qval", "pval_nominal_lead"])
        wpr.writerow(["tissue", "ensembl_id", "variant_id", "tss_distance", "pval_nominal"])
        ntis = collections.Counter()
        for m in tf:
            base = os.path.basename(m.name)
            if not base.endswith(".txt.gz"):
                continue
            tissue = base.split(".v8.")[0]
            with gzip.open(tf.extractfile(m), "rt", encoding="utf-8") as fh:
                h = fh.readline().rstrip("\n").split("\t")
                ix = {c: i for i, c in enumerate(h)}
                if ".egenes." in base:
                    ntis["egenes"] += 1
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        gid = f[ix["gene_id"]]
                        if "_PAR_Y" in gid:
                            continue
                        e = gid.split(".")[0]
                        v26.setdefault(e, (f[ix["gene_name"]], f[ix["gene_chr"]],
                                           f[ix["gene_start"]], f[ix["gene_end"]]))
                        if e in WANT:
                            weg.writerow([tissue, e, gid, f[ix["gene_name"]], f[ix["gene_chr"]],
                                          f[ix["gene_start"]], f[ix["gene_end"]], f[ix["qval"]],
                                          f[ix["pval_nominal"]]])
                elif ".signif_variant_gene_pairs." in base:
                    ntis["pairs"] += 1
                    gi, vi, ti, pi = (ix["gene_id"], ix["variant_id"], ix["tss_distance"],
                                      ix["pval_nominal"])
                    for line in fh:
                        f = line.split("\t", pi + 1)
                        if "_PAR_Y" in f[gi]:
                            continue
                        e = f[gi].split(".")[0]
                        if e in WANT:
                            wpr.writerow([tissue, e, f[vi], f[ti], f[pi]])
            sys.stderr.write("scanned %s\n" % base)
    assert ntis["egenes"] == 49 and ntis["pairs"] == 49, ntis
    with open(V26, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ensembl_id", "gene_name", "chr", "gene_start", "gene_end"])
        for e in sorted(v26):
            w.writerow([e] + list(v26[e]))

v26 = {r["ensembl_id"]: (r["gene_name"], chrom(r["chr"]), int(r["gene_start"]),
                         int(r["gene_end"]))
       for r in csv.DictReader(open(V26, encoding="utf-8"), delimiter="\t")}
tested = collections.defaultdict(dict)          # ens -> tissue -> (chr, start, end)
tissues = set()
for r in csv.DictReader(open(EG, encoding="utf-8"), delimiter="\t"):
    tested[r["ensembl_id"]][r["tissue"]] = (chrom(r["chr"]), int(r["gene_start"]),
                                            int(r["gene_end"]))
gt = collections.defaultdict(dict)              # ens -> tissue -> dict(n, lead, p, n_sig_any)
for r in csv.DictReader(open(PR, encoding="utf-8"), delimiter="\t"):
    tissues.add(r["tissue"])
    c, a, b = tested[r["ensembl_id"]][r["tissue"]]
    vc, vp = r["variant_id"].split("_")[:2]
    d = gt[r["ensembl_id"]].setdefault(r["tissue"], dict(n=0, lead="", p=1.0, n_sig_any=0,
                                                         n_p5e8_any=0))
    d["n_sig_any"] += 1
    p = float(r["pval_nominal"])
    if p < P_MAX:
        d["n_p5e8_any"] += 1
        if chrom(vc) == c and dist(int(vp), a, b) <= WINDOW:
            d["n"] += 1
            if p < d["p"]:
                d["p"], d["lead"] = p, r["variant_id"]
# apply amendment 8
g37 = {}
with gzip.open(os.path.join(ROOT, "data", "raw", "Homo_sapiens.GRCh37.87.gtf.gz"), "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] == "gene" and f[0] in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
            nm = re.search('gene_name "([^"]+)"', f[8]).group(1)
            gi = re.search('gene_id "([^"]+)"', f[8]).group(1)
            g37.setdefault(nm, []).append((f[0], int(f[3]), int(f[4])))
            g37["id:" + gi] = (f[0], int(f[3]), int(f[4]), f[6])
gtex_id = {r["ensembl_id"]: r["ensembl_id"] for r in scope}     # target ID -> ID used in GTEx
idrows = []
for r in scope:
    g, e = r["gene"], r["ensembl_id"]
    if tested.get(e):
        continue
    row = dict(gene=g, grch37_ensembl_id=e, hgnc_id="", hgnc_current_ensembl_id="",
               documented_id_change=False, new_id_in_gtex_v8=False, gtex_gene_name="",
               grch37_coordinates="", gtex_v26_coordinates="", ensembl_grch38_coordinates="",
               same_chromosome="", same_strand="", gene_length_grch37="", gene_length_gtex_v26="",
               gtex_equals_ensembl_grch38="",
               share_of_gtex_gene_body_inside_ensembl_grch38_record="",
               symbol_match_applied=False)
    c37 = g37.get("id:" + e)
    if c37:
        row["grch37_coordinates"] = "%s:%d-%d (%s)" % c37
        row["gene_length_grch37"] = c37[2] - c37[1] + 1
    if g in hgnc_now:
        hid, ne = hgnc_now[g]
        row.update(hgnc_id=hid, hgnc_current_ensembl_id=ne, documented_id_change=True)
        x = json.load(open(os.path.join(IDC, "ensembl38_%s.json" % ne), encoding="utf-8"))
        row["ensembl_grch38_coordinates"] = "%s:%d-%d (%s)" % (
            x["seq_region_name"], x["start"], x["end"], "+" if x["strand"] == 1 else "-")
        if ne in v26:
            nm, c, a, b = v26[ne]
            row.update(new_id_in_gtex_v8=True, gtex_gene_name=nm,
                       gtex_v26_coordinates="%s:%d-%d" % (c, a, b),
                       gene_length_gtex_v26=b - a + 1,
                       same_chromosome=bool(c37) and c == c37[0] == x["seq_region_name"],
                       same_strand=bool(c37) and c37[3] == ("+" if x["strand"] == 1 else "-"),
                       gtex_equals_ensembl_grch38=(c, a, b) == (x["seq_region_name"],
                                                               x["start"], x["end"]))
            ov = max(0, min(b, x["end"]) - max(a, x["start"]) + 1)
            row["share_of_gtex_gene_body_inside_ensembl_grch38_record"] = round(
                ov / float(b - a + 1), 4)
            # the check that gates the match: same symbol, same chromosome in all three
            # records, and >= 95% of the GTEx gene body inside the current GRCh38 record of
            # the HGNC-documented ID. Strand and length versus GRCh37 are reported, not gated.
            ok = (nm == g and row["same_chromosome"] and ov / float(b - a + 1) >= 0.95)
            row["symbol_match_applied"] = ok
            if ok:
                gtex_id[e] = ne
    idrows.append(row)
with open(os.path.join(OUT, "gtex_id_change_resolution.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(idrows[0].keys()))
    w.writeheader()
    w.writerows(idrows)
print("[R02b] targets in no GTEx egenes file under their GRCh37 Ensembl ID: %d; HGNC documents "
      "a different current ID for %d (%s); that ID is in GTEx v8 for %d; symbol match applied "
      "(same symbol and chromosome; >= 95%% of the GTEx gene body inside the current Ensembl "
      "GRCh38 record of that ID) for %d"
      % (len(idrows), sum(r["documented_id_change"] for r in idrows),
         ", ".join(r["gene"] for r in idrows if r["documented_id_change"]),
         sum(r["new_id_in_gtex_v8"] for r in idrows),
         sum(r["symbol_match_applied"] for r in idrows)))
for r in idrows:
    if r["new_id_in_gtex_v8"]:
        print("        %-7s %s -> %s | GRCh37 %s, length %s | GTEx v26 %s, length %s | Ensembl "
              "GRCh38 %s | share of GTEx body inside it %s | same strand as GRCh37: %s | "
              "applied: %s"
              % (r["gene"], r["grch37_ensembl_id"], r["hgnc_current_ensembl_id"],
                 r["grch37_coordinates"], r["gene_length_grch37"], r["gtex_v26_coordinates"],
                 r["gene_length_gtex_v26"], r["ensembl_grch38_coordinates"],
                 r["share_of_gtex_gene_body_inside_ensembl_grch38_record"], r["same_strand"],
                 r["symbol_match_applied"]))


ALL_TISSUES = sorted({t for e in tested for t in tested[e]} | tissues)
xg = sum(1 for e, v in v26.items() if v[1] == "X")
print("[R02] GTEx v8: tissues with target rows %d; genes in any egenes file %d, of which on "
      "chromosome X %d (GTEx v8 covers X)" % (len(ALL_TISSUES), len(v26), xg))

with open(os.path.join(OUT, "recovery_gtex_detail.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "ensembl_id", "gtex_ensembl_id", "tissue", "cv_relevant_tissue",
                "tested_in_tissue",
                "n_significant_pairs", "n_pairs_p5e8_any_distance",
                "n_pairs_p5e8_within_100kb", "instrument", "lead_variant_b38", "lead_p"])
    for r in scope:
        e = gtex_id[r["ensembl_id"]]
        for t in sorted(tested.get(e, {})):
            d = gt.get(e, {}).get(t, dict(n=0, lead="", p="", n_sig_any=0, n_p5e8_any=0))
            w.writerow([r["gene"], r["ensembl_id"], e, t, t in CV_TISSUES, True, d["n_sig_any"],
                        d["n_p5e8_any"], d["n"], d["n"] > 0, d["lead"],
                        ("%.3g" % d["p"]) if d["n"] else ""])

# ------------------------------------------------------------------ (b) pQTL tables


def sheet(path, name, header_has):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    it = wb[name].iter_rows(values_only=True)
    for row in it:
        h = [clean(c).replace("\n", " ") for c in row]
        if header_has in h:
            break
    else:
        raise SystemExit("header %r not found in %s" % (header_has, name))
    out = [dict(zip(h, row)) for row in it if any(c is not None for c in row)]
    wb.close()
    return out


pq = collections.defaultdict(dict)   # gene -> study -> dict
pdetail = []


def record(gene, study, **kw):
    pq[gene][study] = kw


def best(cands):
    """cands: (p, variant, distance or None, extra). Lowest-p qualifying, else lowest-p."""
    ok = sorted((c for c in cands if c[2] is not None and c[2] <= WINDOW), key=lambda c: c[0])
    return (ok[0], True) if ok else ((sorted(cands, key=lambda c: c[0])[0], False)
                                     if cands else (None, False))


def settle(gene, study, assayed, key, cands, n_listed_cis):
    """cands are cis-labelled rows with p < 5e-8."""
    b, ok = best(cands)
    unver = (not ok) and any(c[2] is None for c in cands)
    record(gene, study, assayed=assayed, match=key, n_cis_rows=n_listed_cis,
           instrument=ok, unverified=unver,
           lead=b[1] if b else "", p=b[0] if b else "", dist=("" if not b or b[2] is None
                                                             else b[2]),
           extra=b[3] if b else "")
    for c in cands:
        pdetail.append([gene, study, key, c[1], c[0], "" if c[2] is None else c[2],
                        c[2] is not None and c[2] <= WINDOW, c[3]])


# --- UKB-PPP (Sun 2023): ST3 proteins (GRCh38 gene coordinates), ST9 discovery, ST10 combined
UKB = os.path.join(DATA, "41586_2023_6592_MOESM3_ESM.xlsx")
u3 = sheet(UKB, "ST3", "UKBPPP ProteinID")
uprot = collections.defaultdict(list)          # symbol -> [(protein id, chr, start, end)]
multi = set()
for r in u3:
    syms = tokens(r["Gene symbol"].replace("_", "|")) if "_" in clean(r["Gene symbol"]) \
        else tokens(r["Gene symbol"])
    if len(syms) > 1:
        multi.add(clean(r["Gene symbol"]))
    cs, ss, es = tokens(r["Gene CHROM"]), tokens(r["Gene start"]), tokens(r["Gene end"])
    for i, s in enumerate(syms):
        try:
            j = i if len(ss) == len(syms) else 0
            co = (chrom(cs[j if len(cs) == len(syms) else 0]), int(float(ss[j])),
                  int(float(es[j])))
        except (ValueError, IndexError):
            co = None
        uprot[s].append((clean(r["UKBPPP ProteinID"]), co))
print("[R03] UKB-PPP ST3: %d protein rows, %d distinct gene symbols; multi-gene entries: %s"
      % (len(u3), len(uprot), "; ".join(sorted(multi)) or "none"))
ucis = collections.defaultdict(list)
uX = 0
for sh in ("ST9", "ST10"):
    rr = sheet(UKB, sh, "UKBPPP ProteinID")
    nc = 0
    for r in rr:
        if clean(r["cis/trans"]).lower() != "cis":
            continue
        nc += 1
        uX += chrom(r["CHROM"]) == "X"
        ucis[clean(r["UKBPPP ProteinID"])].append(
            (sh, chrom(r["CHROM"]), int(float(r["GENPOS (hg38)"])),
             float(r["log10(p)" if sh == "ST10" else "log10(p) (discovery)"]),
             clean(r["rsID"]) or clean(r[[k for k in r if k.startswith("Variant ID")][0]]),
             clean(r["cis gene"])))
    print("[R04] UKB-PPP %s: %d rows, %d labelled cis" % (sh, len(rr), nc))
print("[R05] UKB-PPP cis rows on chromosome X (both sheets): %d" % uX)
for r in scope:
    g = r["gene"]
    if g not in uprot:
        record(g, "UKB-PPP", assayed=False, match="", n_cis_rows=0, instrument=False,
               unverified=False, lead="", p="", dist="", extra="")
        continue
    cands, n = [], 0
    for pid, co in uprot[g]:
        for sh, c, pos, lp, var, cg in ucis.get(pid, []):
            n += 1
            if lp <= LOG10_MIN:
                continue
            d = dist(pos, co[1], co[2]) if (co and co[0] == c) else None
            cands.append((10 ** -lp if lp < 300 else 0.0, "%s (chr%s:%d)" % (var, c, pos), d,
                          sh))
    settle(g, "UKB-PPP", True, "symbol", cands, n)

# --- deCODE (Ferkingstad 2021): ST01 aptamers, ST02 sentinel + secondary pQTLs (GRCh38)
DEC = os.path.join(DATA, "ferkingstad2021_MOESM4.xlsx")
d1 = sheet(DEC, "ST01", "SeqId")
dapt = collections.defaultdict(set)            # gene -> {(seqid, matchkey)}
dinc = collections.Counter()
for r in d1:
    inc = clean(r["Included in analysis"])
    dinc[inc] += 1
    if inc.lower() != "yes":
        continue
    for s in tokens(r["Gene"]):
        if s in by_sym:
            dapt[s].add((clean(r["SeqId"]), "symbol"))
    for e in tokens(r["Ensembl.Gene.ID"]):
        if e in by_ens:
            dapt[by_ens[e]["gene"]].add((clean(r["SeqId"]), "ensembl_id"))
print("[R06] deCODE ST01: %d aptamer rows; 'Included in analysis': %s" % (len(d1), dict(dinc)))
d2 = sheet(DEC, "ST02", "SeqId")
dcis = collections.defaultdict(list)
dX = 0
for r in d2:
    if clean(r["cis/ trans"]).lower() != "cis":
        continue
    dX += chrom(r["chr (var.)"]) == "X"
    dcis[clean(r["SeqId"])].append((chrom(r["chr (var.)"]), int(float(r["pos (var.)"])),
                                    float(r["-Log10(P) (adj.)"]), clean(r["variant"]),
                                    clean(r["distTSS"])))
print("[R07] deCODE ST02: %d rows, %d labelled cis, %d of those on chromosome X"
      % (len(d2), sum(len(v) for v in dcis.values()), dX))
for r in scope:
    g, e = r["gene"], r["ensembl_id"]
    if g not in dapt:
        record(g, "deCODE", assayed=False, match="", n_cis_rows=0, instrument=False,
               unverified=False, lead="", p="", dist="", extra="")
        continue
    co = v26.get(gtex_id[e])
    cands, n = [], 0
    for sid in sorted({s for s, _ in dapt[g]}):
        for c, pos, lp, var, dtss in dcis.get(sid, []):
            n += 1
            if lp <= LOG10_MIN:
                continue
            d = dist(pos, co[2], co[3]) if (co and co[1] == c) else None
            cands.append((10 ** -lp if lp < 300 else 0.0, "%s (chr%s:%d)" % (var, c, pos), d,
                          "aptamer %s; distTSS %s" % (sid, dtss)))
    settle(g, "deCODE", True, "+".join(sorted({k for _, k in dapt[g]})), cands, n)

# --- Fenland (Pietzner 2021): ST3 aptamers with gene coordinates, ST2 pQTLs
FEN = os.path.join(DATA, "pietzner2021_supplementary_tables.xlsx")
f3 = sheet(FEN, "ST3", "SomaScan.ID")
n38 = collections.defaultdict(list)
for e, (nm, c, a, b) in v26.items():
    n38[nm].append((c, a, b))
fapt = collections.defaultdict(list)           # symbol -> [(somascan id, chr, start, end)]
fsid = collections.defaultdict(list)           # somascan id -> [(symbol, (chr, start, end))]
cmp_ = collections.Counter()
for r in f3:
    syms = tokens(r["Protein encoding gene"])
    cs, ss, es = tokens(r["Chr"]), tokens(r["start"]), tokens(r["end"])
    for i, s in enumerate(syms):
        try:
            j = i if len(ss) == len(syms) else 0
            co = (chrom(cs[j if len(cs) == len(syms) else 0]), int(float(ss[j])),
                  int(float(es[j])))
        except (ValueError, IndexError):
            co = None
        fapt[s].append((clean(r["SomaScan.ID"]), co))
        fsid[clean(r["SomaScan.ID"])].append((s, co))
        if co and len(syms) == 1:
            cmp_["compared"] += 1
            cmp_["exact GRCh37.87"] += co in g37.get(s, [])
            cmp_["exact GENCODE v26 (GRCh38)"] += co in n38.get(s, [])
print("[R08] Fenland ST3: %d aptamer rows, %d gene symbols. Build check on single-gene "
      "aptamers with coordinates: %d compared; start and end identical to Ensembl GRCh37.87 "
      "for %d; identical to GENCODE v26 (GRCh38) for %d"
      % (len(f3), len(fapt), cmp_["compared"], cmp_["exact GRCh37.87"],
         cmp_["exact GENCODE v26 (GRCh38)"]))
FEN_BUILD_OK = cmp_["exact GRCh37.87"] >= 0.9 * cmp_["compared"]
print("[R08] Fenland build taken as GRCh37: %s" % FEN_BUILD_OK)
f2 = sheet(FEN, "ST2", "SomaScan.id")
fcis = collections.defaultdict(list)           # somascan id -> rows
fmatch = collections.defaultdict(set)
fX = 0
flab = collections.Counter(clean(r["cis.trans"]).lower() for r in f2)
fblank = [(clean(r["HGNC.symbol.protein"]), clean(r["cis.trans"])) for r in f2
          if not clean(r["P.value"])]
print("[R09] Fenland ST2 cis.trans labels: %s; rows with a blank p-value (not counted): %s; "
      "in-scope targets among them: %s"
      % (dict(flab), fblank, [g for g, _ in fblank if g in by_sym] or "none"))


def fen_p(v):
    s = clean(v)
    if not s:                      # blank p-value in the table: not estimated, not counted
        return float("nan")
    return 0.0 if s.startswith("<") else float(s)


for r in f2:
    if clean(r["cis.trans"]).lower() != "cis":
        continue
    sid = clean(r["SomaScan.id"])
    fX += chrom(r["Chr"]) == "X"
    fcis[sid].append((chrom(r["Chr"]), int(float(r["Position"])), fen_p(r["P.value"]),
                      clean(r["rsID"]) or clean(r["MarkerName"]), clean(r["signal.locus"])))
    for s in tokens(r["HGNC.symbol.protein"]):
        if s in by_sym:
            fmatch[s].add((sid, "symbol"))
    for e in tokens(r["ENSEMBL.id.protein"]):
        if e in by_ens:
            fmatch[by_ens[e]["gene"]].add((sid, "ensembl_id"))
print("[R09] Fenland ST2: %d rows, %d labelled cis, %d of those on chromosome X"
      % (len(f2), sum(len(v) for v in fcis.values()), fX))
for r in scope:
    g = r["gene"]
    apt = {s: co for s, co in fapt.get(g, [])}
    for sid, _ in fmatch.get(g, []):
        own = [co for s_, co in fsid.get(sid, []) if co]
        apt.setdefault(sid, own[0] if len(own) == 1 else None)
    if not apt:
        record(g, "Fenland", assayed=False, match="", n_cis_rows=0, instrument=False,
               unverified=False, lead="", p="", dist="", extra="")
        continue
    keys = set(k for _, k in fmatch.get(g, [])) | ({"symbol"} if g in fapt else set())
    cands, n = [], 0
    for sid in sorted(apt):
        co = apt[sid]
        for c, pos, p, var, rank in fcis.get(sid, []):
            n += 1
            if not p < P_MAX:
                continue
            d = dist(pos, co[1], co[2]) if (FEN_BUILD_OK and co and co[0] == c) else None
            cands.append((p, "%s (chr%s:%d)" % (var, c, pos), d,
                          "aptamer %s; signal rank %s" % (sid, rank)))
    settle(g, "Fenland", True, "+".join(sorted(keys)), cands, n)

with open(os.path.join(OUT, "recovery_pqtl_detail.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "study", "matched_by", "variant_and_position_in_table_build", "p",
                "bp_from_gene_body", "within_100kb", "note"])
    w.writerows(sorted(pdetail, key=lambda x: (x[0], x[1], x[4])))

# ------------------------------------------------------------------ categories
STUDIES = ["UKB-PPP", "deCODE", "Fenland"]
out = []
for r in scope:
    g, e0 = r["gene"], r["ensembl_id"]
    e = gtex_id[e0]
    T0 = tested.get(e0, {})
    before = None
    if e != e0:                      # category under the rule before amendment 8
        pin0 = [s for s in STUDIES if pq[g][s]["instrument"]]
        before = "fixable-protein" if pin0 else "locked"
    T = tested.get(e, {})
    inst = sorted(t for t, d in gt.get(e, {}).items() if d["n"] > 0)
    nonblood = [t for t in inst if t != BLOOD]
    cvt = [t for t in CV_TISSUES if t in inst]
    lead_t = min(nonblood, key=lambda t: gt[e][t]["p"]) if nonblood else ""
    lead_cv = min(cvt, key=lambda t: gt[e][t]["p"]) if cvt else ""
    pin = [s for s in STUDIES if pq[g][s]["instrument"]]
    unv = [s for s in STUDIES if pq[g][s]["unverified"]]
    cat = ("fixable-both" if nonblood and pin else "fixable-tissue" if nonblood
           else "fixable-protein" if pin else "locked")
    assayed = [s for s in STUDIES if pq[g][s]["assayed"]]
    tier = ("plausible rescue" if (cvt or pin) else "off-target-tissue rescue" if nonblood
            else "locked: never assayed" if not assayed else "locked: assayed, no instrument")
    flags = []
    if e != e0:
        flags.append("GTEx matched by symbol after documented Ensembl ID change (%s)" % e)
    if r["chromosome"] == "X":
        flags.append("X")
    if cat == "locked" and BLOOD in inst:
        flags.append("GTEx whole blood only")
    if unv and not pin:
        flags.append("category depends on UNVERIFIED cell (%s)" % ", ".join(unv))
    if not T:
        flags.append("gene not in any GTEx v8 egenes file")
    row = collections.OrderedDict(
        gene=g, ensembl_id=e0, gtex_ensembl_id=e if T else "", chromosome=r["chromosome"],
        x_linked=r["chromosome"] == "X",
        eqtlgen_status=r["status"],
        eqtlgen_subtype=r["no_eqtl_data_subtype"] or r["no_instrument_subtype"],
        phase_primary=r["phase_PRIMARY"], indication_areas=r["indication_areas"],
        chembl_verified_strict=r["chembl_verified_strict"], named_13=g in NAMED,
        category=cat, reporting_tier=tier,
        category_before_id_amendment=before or cat, flags="; ".join(flags),
        gtex_tissues_with_instrument_excl_blood_with_p="; ".join(
            "%s (%.2g)" % (t, gt[e][t]["p"])
            for t in sorted(nonblood, key=lambda t: gt[e][t]["p"])),
        gtex_n_tissues_tested=len(T),
        gtex_n_tissues_with_instrument_excl_blood=len(nonblood),
        gtex_tissues_with_instrument_excl_blood="; ".join(nonblood),
        gtex_cv_relevant_tissues_with_instrument="; ".join(cvt),
        gtex_cv_tissue_fixable=bool(cvt),
        gtex_cv_tissues_tested="; ".join(t for t in CV_TISSUES if t in T),
        gtex_best_tissue=lead_t,
        gtex_best_variant_b38=gt[e][lead_t]["lead"] if lead_t else "",
        gtex_best_p=("%.3g" % gt[e][lead_t]["p"]) if lead_t else "",
        gtex_best_cv_tissue=lead_cv,
        gtex_best_cv_variant_b38=gt[e][lead_cv]["lead"] if lead_cv else "",
        gtex_best_cv_p=("%.3g" % gt[e][lead_cv]["p"]) if lead_cv else "",
        gtex_whole_blood_tested=BLOOD in T, gtex_whole_blood_instrument=BLOOD in inst,
        gtex_whole_blood_p=("%.3g" % gt[e][BLOOD]["p"]) if BLOOD in inst else "")
    for s in STUDIES:
        k, d = s.lower().replace("-", ""), pq[g][s]
        row[k + "_assayed"] = d["assayed"]
        row[k + "_matched_by"] = d["match"]
        row[k + "_n_cis_rows_listed"] = d["n_cis_rows"]
        row[k + "_instrument"] = d["instrument"]
        row[k + "_cell"] = ("instrument" if d["instrument"] else
                            "cis-pQTL reported; distance UNVERIFIED" if d["unverified"] else
                            "cis-pQTL reported; outside +-100 kb" if d["lead"] else
                            "assayed; no cis-pQTL in table" if d["assayed"] else "not assayed")
        row[k + "_lead_variant"] = d["lead"]
        row[k + "_p"] = ("" if not d["lead"] else "<1e-300" if d["p"] == 0
                         else "%.3g" % d["p"])
        row[k + "_bp_from_gene_body"] = d["dist"]
        row[k + "_note"] = d["extra"]
    row["pqtl_studies_with_instrument"] = "; ".join(pin)
    row["pqtl_studies_assayed"] = "; ".join(s for s in STUDIES if pq[g][s]["assayed"])
    out.append(row)
TIERS = ["plausible rescue", "off-target-tissue rescue", "locked: never assayed",
         "locked: assayed, no instrument"]
out.sort(key=lambda x: (TIERS.index(x["reporting_tier"]), x["gene"]))
with open(os.path.join(MS, "recovery_map.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
src = {r["gene"]: r for r in scope}

# ------------------------------------------------------------------ summaries


KEY = "category"


def table(tag, title, groups):
    print("\n[%s] %s" % (tag, title))
    print("  %-46s %5s  %14s %15s %12s %6s   %s" % ("", "n", *CATS, "fixable (any)"))
    for name, rs in groups:
        c = collections.Counter(x[KEY] for x in rs)
        n = len(rs)
        fx = n - c["locked"]
        if KEY != "category":
            raise SystemExit("table() is for categories; use ttable() for tiers")
        x = sum(1 for q in rs if q["x_linked"])
        xl = sum(1 for q in rs if q["x_linked"] and q["category"] == "locked")
        print("  %-46s %5d  %14d %15d %12d %6d   %d (%s)   [X-linked: %d, of which locked %d]"
              % (name[:46], n, *(c[k] for k in CATS), n - c["locked"] if KEY == "category"
                 else fx,
                 ("%.1f%%" % (100.0 * fx / n)) if n else "-", x, xl))


cc = collections.Counter((r["no_eqtl_data_subtype"], r["chromosome"]) for r in scope
                         if r["chromosome"] in ("X", "Y", "MT"))
print("[R09b] CONFLICT CHECK, sub-type label vs chromosome for in-scope non-autosomal "
      "targets: %s" % dict(cc))
table("R10", "1. Overall", [("all in-scope targets", out)])
table("R11", "2. By eQTLGen status", [
    ("no eQTL data: autosomal, not tested",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA"
      and x["chromosome"] not in ("X", "Y", "MT")]),
    ("no eQTL data: X-linked, not analysed",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA" and x["x_linked"]]),
    ("no eQTL data: mitochondrial genome (MT)",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA" and x["chromosome"] == "MT"]),
    ("no qualifying instrument",
     [x for x in out if x["eqtlgen_status"] == "NO_QUALIFYING_INSTRUMENT"])])
table("R12", "3. By highest clinical phase in the primary set",
      [(p, [x for x in out if x["phase_primary"] == p]) for p in PHASES])
areas = [k[len("area_phase: "):] for k in rows[0] if k.startswith("area_phase: ")]
table("R13", "4. By indication area (a target is counted in each area it appears in)",
      [(a, [x for x in out if src[x["gene"]]["area_phase: " + a]]) for a in areas])
table("R14", "Strict-verified launched targets (PROVISIONAL pending manual review)",
      [("strict-verified, Launched in primary set",
        [x for x in out if x["chembl_verified_strict"] == "True"
         and x["phase_primary"] == "Launched"])])

print("\n[R15] 5. The 13 named targets")
for g in NAMED:
    x = next(q for q in out if q["gene"] == g)
    print("  %-8s %-16s %s%s" % (g, x["category"], ("[" + x["flags"] + "] ") if x["flags"]
                                 else "", ""))
    print("           tier: %s" % x["reporting_tier"])
    print("           GTEx non-blood tissues with instrument (%d), lead p: %s"
          % (x["gtex_n_tissues_with_instrument_excl_blood"],
             x["gtex_tissues_with_instrument_excl_blood_with_p"] or "none"))
    print("           of which CV-relevant: %s; GTEx whole blood instrument: %s"
          % (x["gtex_cv_relevant_tissues_with_instrument"] or "none",
             x["gtex_whole_blood_instrument"]))
    for s in STUDIES:
        k = s.lower().replace("-", "")
        print("           %-8s %s%s" % (s, x[k + "_cell"],
              (": %s, p=%s, %s bp from gene body" % (x[k + "_lead_variant"], x[k + "_p"],
                                                     x[k + "_bp_from_gene_body"]))
              if x[k + "_lead_variant"] else ""))



def ttable(tag, title, groups):
    print("\n[%s] %s" % (tag, title))
    print("  %-46s %5s  %16s %24s %21s %30s" % ("", "n", *TIERS))
    for name, rs in groups:
        c = collections.Counter(x["reporting_tier"] for x in rs)
        n = len(rs)
        print("  %-46s %5d  %16s %24s %21s %30s   [locked, total %d (%s)]"
              % (name[:46], n, *(("%d (%.1f%%)" % (c[k], 100.0 * c[k] / n)) if n else "0"
                                 for k in TIERS), c[TIERS[2]] + c[TIERS[3]],
                 ("%.1f%%" % (100.0 * (c[TIERS[2]] + c[TIERS[3]]) / n)) if n else "-"))


ttable("R23", "REPORTING TIERS, overall and by eQTLGen status", [
    ("all in-scope targets", out),
    ("no eQTL data: autosomal, not tested",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA"
      and x["chromosome"] not in ("X", "Y", "MT")]),
    ("no eQTL data: X-linked, not analysed",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA" and x["x_linked"]]),
    ("no eQTL data: mitochondrial genome (MT)",
     [x for x in out if x["eqtlgen_status"] == "NO_EQTL_DATA" and x["chromosome"] == "MT"]),
    ("no qualifying instrument",
     [x for x in out if x["eqtlgen_status"] == "NO_QUALIFYING_INSTRUMENT"])])
ttable("R24", "REPORTING TIERS by highest clinical phase in the primary set",
       [(p, [x for x in out if x["phase_primary"] == p]) for p in PHASES])
ttable("R25", "REPORTING TIERS by indication area (a target is counted in each area it "
       "appears in)",
       [(a, [x for x in out if src[x["gene"]]["area_phase: " + a]]) for a in areas])
ttable("R26", "REPORTING TIERS, the 13 named targets and the strict-verified launched targets "
       "(PROVISIONAL)",
       [("the 13 named targets", [x for x in out if x["named_13"]]),
        ("strict-verified, Launched in primary set",
         [x for x in out if x["chembl_verified_strict"] == "True"
          and x["phase_primary"] == "Launched"])])
pr = [x for x in out if x["reporting_tier"] == "plausible rescue"]
print("\n[R27] plausible rescue (%d): through a CV-relevant tissue only %d; through a plasma "
      "cis-pQTL only %d; through both %d"
      % (len(pr),
         sum(1 for x in pr if x["gtex_cv_tissue_fixable"] and not x["pqtl_studies_with_instrument"]),
         sum(1 for x in pr if not x["gtex_cv_tissue_fixable"] and x["pqtl_studies_with_instrument"]),
         sum(1 for x in pr if x["gtex_cv_tissue_fixable"] and x["pqtl_studies_with_instrument"])))
ot = [x for x in out if x["reporting_tier"] == "off-target-tissue rescue"]
print("[R28] off-target-tissue rescues (%d), with every tissue that has an instrument:"
      % len(ot))
for x in ot:
    print("        %-9s %-11s %s" % (x["gene"], x["phase_primary"],
                                     x["gtex_tissues_with_instrument_excl_blood"]))
tc = collections.Counter(t for x in ot
                         for t in x["gtex_tissues_with_instrument_excl_blood"].split("; "))
print("[R29] off-target-tissue rescues: tissue frequency: %s"
      % "; ".join("%s %d" % kv for kv in tc.most_common()))
print("[R29] off-target-tissue rescues resting on a single tissue: %d; on testis only: %d"
      % (sum(1 for x in ot if x["gtex_n_tissues_with_instrument_excl_blood"] == 1),
         sum(1 for x in ot if x["gtex_tissues_with_instrument_excl_blood"] == "Testis")))
with open(os.path.join(OUT, "off_target_tissue_rescues.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "phase_primary", "indication_areas", "eqtlgen_status",
                "n_tissues", "tissues_with_instrument_and_lead_p", "pqtl_studies_assayed"])
    for x in ot:
        w.writerow([x["gene"], x["phase_primary"], x["indication_areas"], x["eqtlgen_status"],
                    x["gtex_n_tissues_with_instrument_excl_blood"],
                    x["gtex_tissues_with_instrument_excl_blood_with_p"],
                    x["pqtl_studies_assayed"]])
for k in TIERS[2:]:
    print("[R30] %s (%d): %s" % (k, sum(1 for x in out if x["reporting_tier"] == k),
                                 ", ".join("%s%s" % (x["gene"], "(X)" if x["x_linked"] else "")
                                           for x in out if x["reporting_tier"] == k)))
print("[R31] amendment 8, before -> after: %s"
      % "; ".join("%s: %s -> %s (%s)" % (x["gene"], x["category_before_id_amendment"],
                                         x["category"], x["reporting_tier"])
                  for x in out if "documented Ensembl ID change" in x["flags"]))

fixt = [x for x in out if x["gtex_n_tissues_with_instrument_excl_blood"] > 0]
cvf = [x for x in fixt if x["gtex_cv_tissue_fixable"]]
print("\n[R16] 6. GTEx: targets with an instrument in >=1 non-blood tissue: %d; through a "
      "CV-relevant tissue: %d; only through another tissue: %d"
      % (len(fixt), len(cvf), len(fixt) - len(cvf)))
for t in CV_TISSUES:
    print("        %-24s tested %3d   instrument %3d"
          % (t, sum(1 for x in out if t in x["gtex_cv_tissues_tested"].split("; ")),
             sum(1 for x in out if t in x["gtex_cv_relevant_tissues_with_instrument"].split("; "))))
print("[R17] GTEx: targets in no egenes file (not tested in any tissue): %d; tested in >=1 "
      "tissue: %d" % (sum(1 for x in out if x["gtex_n_tissues_tested"] == 0),
                      sum(1 for x in out if x["gtex_n_tissues_tested"] > 0)))
print("[R18] GTEx Whole Blood: tested %d; instrument %d; locked targets flagged 'GTEx whole "
      "blood only': %d"
      % (sum(1 for x in out if x["gtex_whole_blood_tested"]),
         sum(1 for x in out if x["gtex_whole_blood_instrument"]),
         sum(1 for x in out if "GTEx whole blood only" in x["flags"])))
print("[R19] pQTL tables, per study (of %d in-scope targets):" % len(out))
for s in STUDIES:
    k = s.lower().replace("-", "")
    c = collections.Counter(x[k + "_cell"] for x in out)
    print("        %-8s %s" % (s, "; ".join("%s %d" % kv for kv in sorted(c.items()))))
print("[R20] pQTL: assayed in >=1 study %d; instrument in >=1 study %d; in all three %d; "
      "targets flagged UNVERIFIED-dependent %d"
      % (sum(1 for x in out if x["pqtl_studies_assayed"]),
         sum(1 for x in out if x["pqtl_studies_with_instrument"]),
         sum(1 for x in out if x["pqtl_studies_with_instrument"].count(";") == 2),
         sum(1 for x in out if "UNVERIFIED" in x["flags"])))
lk = [x for x in out if x["category"] == "locked"]
print("[R21] locked (relative to GTEx v8 non-blood tissues and the UKB-PPP, deCODE and "
      "Fenland published cis-pQTL tables): %d. Of these, assayed in no pQTL study: %d; not in "
      "any GTEx egenes file: %d"
      % (len(lk), sum(1 for x in lk if not x["pqtl_studies_assayed"]),
         sum(1 for x in lk if x["gtex_n_tissues_tested"] == 0)))
print("[R22] locked, Launched in primary set: %s"
      % ", ".join(sorted("%s%s" % (x["gene"], "*" if x["chembl_verified_strict"] == "True"
                                   else "") for x in lk if x["phase_primary"] == "Launched")))
print("       (* = strict-verified)")
