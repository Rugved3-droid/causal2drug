"""AUDIT step 2a (part 1): CV target set and gene-record mapping, from rules alone.

Rules used:
  DECISIONS.md  DECISION 1        dedup pp.tsv by ti_uid, drop 4 conflicting programmes
  CV_SET_RULE.md                  areas contains 'cardiovascular' OR term contains a lipid token
  DEFINITIONS.md Part 1           one GRCh37 gene record per Minikel symbol (Ensembl 87 GTF)
  UNMAPPED_RESOLUTION.md          HGNC approved record -> Ensembl ID -> must be a GTF gene record
  DEFINITIONS.md Part 2           phase(target, S) = max ccatnum over pairs in S

Run:  python audit/a01_cv_set_and_coords.py
Writes audit/out/cv_terms.csv, cv_targets.csv, gene_coords_audit.tsv
"""
import csv, gzip, collections, json, os, re, sys, urllib.request

ROOT = r"D:\causal2drug"
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "audit", "out")
os.makedirs(OUT, exist_ok=True)

DROP = {"DAO-D022124", "HRH1-D000086382", "IFNAR2-D012141", "MC1R-D010787"}
TOKENS = ["hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
          "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler"]

# ---- pairs -----------------------------------------------------------------
rows = list(csv.DictReader(open(os.path.join(RAW, "pp.tsv"), encoding="utf-8"), delimiter="\t"))
seen, pairs = set(), []
for r in rows:
    if r["ti_uid"] in DROP or r["ti_uid"] in seen:
        continue
    seen.add(r["ti_uid"])
    pairs.append(r)
print(f"[A01] raw rows {len(rows)}; after dedup+drop {len(pairs)}; "
      f"genes {len({p['gene'] for p in pairs})}; indications {len({p['indication_mesh_id'] for p in pairs})}")

# ---- CV terms --------------------------------------------------------------
indic = list(csv.DictReader(open(os.path.join(RAW, "minikel_indic.tsv"), encoding="utf-8"), delimiter="\t"))
cv_terms = {}
for r in indic:
    areas = [a.strip().lower() for a in r["areas"].split(",")]
    by_area = "cardiovascular" in areas
    tok = [t for t in TOKENS if t in r["indication_mesh_term"].lower()]
    if by_area or tok:
        cv_terms[r["indication_mesh_id"]] = dict(term=r["indication_mesh_term"], areas=r["areas"],
                                                 by_area=by_area, tokens=";".join(tok))
cv_pairs = [p for p in pairs if p["indication_mesh_id"] in cv_terms]
used_terms = {p["indication_mesh_id"] for p in cv_pairs}
phase = collections.defaultdict(int)
for p in cv_pairs:
    phase[p["gene"]] = max(phase[p["gene"]], int(p["ccatnum"] or 0))
targets = sorted(phase)
print(f"[A02] CV terms by rule {len(cv_terms)} (with >=1 pair: {len(used_terms)}); "
      f"token-only terms {sum(1 for v in cv_terms.values() if not v['by_area'])}; "
      f"pairs {len(cv_pairs)}; TARGETS {len(targets)}")
print(f"[A03] launched (phase 5) targets in primary set: {sum(1 for t in targets if phase[t]==5)}")

npairs = collections.Counter(p["indication_mesh_id"] for p in cv_pairs)
with open(os.path.join(OUT, "cv_terms.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["mesh_id", "term", "areas", "by_area", "tokens", "pairs"])
    for k, v in sorted(cv_terms.items()):
        w.writerow([k, v["term"], v["areas"], v["by_area"], v["tokens"], npairs[k]])

# ---- GTF gene records ------------------------------------------------------
by_name, by_id = collections.defaultdict(list), {}
with gzip.open(os.path.join(RAW, "Homo_sapiens.GRCh37.87.gtf.gz"), "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        rec = dict(ensembl_id=gid, chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], gtf_name=name)
        by_name[name].append(rec); by_id[gid] = rec
print(f"[A04] GTF gene records {len(by_id)}")

# ---- HGNC (complete set, fetched once and cached in audit/out) --------------
hg_path = os.path.join(OUT, "hgnc_complete_set.txt")
if not os.path.exists(hg_path):
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    urllib.request.urlretrieve(url, hg_path)
hgnc = {}
for r in csv.DictReader(open(hg_path, encoding="utf-8"), delimiter="\t"):
    if r["status"] == "Approved":
        hgnc[r["symbol"]] = r
print(f"[A05] HGNC approved records {len(hgnc)}")

coords, how = {}, {}
for t in targets:
    recs = by_name.get(t, [])
    if len(recs) == 1:
        coords[t], how[t] = recs[0], "gtf_unique_symbol"
        continue
    h = hgnc.get(t)
    eid = h["ensembl_gene_id"] if h else ""
    if eid and eid in by_id:
        coords[t], how[t] = by_id[eid], f"hgnc_fix({len(recs)}_gtf_records)"
    else:
        how[t] = f"UNRESOLVED({len(recs)}_gtf_records;hgnc_id={eid or 'none'})"
c = collections.Counter(v.split("(")[0] for v in how.values())
print(f"[A06] mapping: {dict(c)}")
print(f"[A06] unresolved: {sorted(t for t in targets if t not in coords)}")
print(f"[A06] hgnc-fixed: {sorted(t for t in targets if how[t].startswith('hgnc_fix'))}")

with open(os.path.join(OUT, "gene_coords_audit.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["symbol", "phase", "ensembl_id", "chrom", "start", "end", "strand", "gtf_name", "how",
                "hgnc_current_ensembl_id"])
    for t in targets:
        c_ = coords.get(t)
        h = hgnc.get(t)
        w.writerow([t, phase[t]] + ([c_["ensembl_id"], c_["chrom"], c_["start"], c_["end"], c_["strand"],
                                     c_["gtf_name"]] if c_ else [""] * 6)
                   + [how[t], h["ensembl_gene_id"] if h else ""])
