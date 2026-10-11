"""AUDIT step 2c: recovery map for the uninstrumentable targets, from RECOVERY_RULE.md alone
(original rule + clarifications 1-7 + post-hoc amendments 8 and 9).

Instrument everywhere: variant p < 5e-8 inside gene body +-100,000 bp, in the resource's own build.

Run:  python audit/a05_recovery_map.py
Writes audit/out/audit_recovery_map.csv
"""
import csv, collections, json, math, os, urllib.request

OUT = r"D:\causal2drug\audit\out"
PQ = os.path.join(OUT, "pqtl")
FLANK, P, LOGP = 100_000, 5e-8, -math.log10(5e-8)
CV_TISSUES = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Aorta", "Artery_Coronary",
              "Artery_Tibial", "Liver", "Kidney_Cortex", "Adrenal_Gland"]
BLOOD = "Whole_Blood"


def table(path, skip, delim="\t"):
    r = list(csv.reader(open(path, encoding="utf-8"), delimiter=delim))
    h = [c.replace("\n", " ").strip() for c in r[skip]]
    return h, [x for x in r[skip + 1:] if any(x)]


def clean(s):
    return str(s).strip().lstrip("'").strip()


# ---- scope -----------------------------------------------------------------
targets = [r for r in csv.DictReader(open(os.path.join(OUT, "audit_per_target.csv"), encoding="utf-8"))]
scope = [r for r in targets if r["status"] in ("NO_EQTL_DATA", "NO_QUALIFYING_INSTRUMENT")]
hg_cur = {r["symbol"]: r["hgnc_current_ensembl_id"]
          for r in csv.DictReader(open(os.path.join(OUT, "gene_coords_audit.tsv"), encoding="utf-8"), delimiter="\t")}
print(f"[C01] in scope: {len(scope)} (no eQTL data {sum(r['status']=='NO_EQTL_DATA' for r in scope)}, "
      f"no qualifying instrument {sum(r['status']=='NO_QUALIFYING_INSTRUMENT' for r in scope)}); "
      f"out of scope unresolved: {[r['symbol'] for r in targets if r['status']=='UNRESOLVED']}")

# ---- GTEx ------------------------------------------------------------------
egenes = collections.defaultdict(dict)        # id -> tissue -> (name, chr, start, end, strand)
tissues = set()
with open(os.path.join(OUT, "gtex_egenes_all.tsv"), encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        t, gid, name, c, s, e, strand = line.rstrip("\n").split("\t")
        egenes[gid.split(".")[0]][t] = (name, c, int(s), int(e), strand); tissues.add(t)
print(f"[C02] GTEx tissues {len(tissues)}; genes in >=1 egenes file {len(egenes)}")
signif = collections.defaultdict(list)        # id -> [(tissue, variant, p, slope, se)]
with open(os.path.join(OUT, "gtex_signif_targets.tsv"), encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        t, v, gid, p, b, se = line.rstrip("\n").split("\t")
        signif[gid.split(".")[0]].append((t, v, float(p), b, se))


def ensembl_lookup(eid):
    cache = os.path.join(OUT, f"ensembl_lookup_{eid}.json")
    if not os.path.exists(cache):
        req = urllib.request.Request(f"https://rest.ensembl.org/lookup/id/{eid}?content-type=application/json")
        try:
            data = urllib.request.urlopen(req, timeout=60).read().decode()
        except Exception as ex:
            data = json.dumps({"error": str(ex)})
        open(cache, "w").write(data)
    return json.load(open(cache))


gtex_id, a8 = {}, []
for r in scope:
    s, eid = r["symbol"], r["ensembl_id"]
    if eid in egenes:
        gtex_id[s] = (eid, "grch37_id")
        continue
    cur = hg_cur.get(s, "")
    note = "not_in_gtex"
    if cur and cur != eid:                                     # (i) documented ID change
        if cur in egenes:
            recs = list(egenes[cur].values())
            same_symbol = all(x[0] == s for x in recs)             # (ii)
            g_chr = recs[0][1].replace("chr", ""); g_s = min(x[2] for x in recs); g_e = max(x[3] for x in recs)
            ens = ensembl_lookup(cur)
            ok_chr = (not ens.get("error")) and str(ens.get("seq_region_name")) == g_chr == r["chrom"]   # (iii)
            frac = 0.0
            if not ens.get("error"):
                ov = max(0, min(g_e, ens["end"]) - max(g_s, ens["start"]) + 1)
                frac = ov / (g_e - g_s + 1)                    # (iv)
            passed = same_symbol and ok_chr and frac >= 0.95
            note = (f"id_change {eid}->{cur}; same_symbol={same_symbol}; chr_ok={ok_chr}; inside={frac:.4f}; "
                    f"strand37={r.get('strand','?')}/gtex={recs[0][4]}; {'APPLIED' if passed else 'NOT applied'}")
            if passed:
                gtex_id[s] = (cur, "amendment8")
        else:
            note = f"id_change {eid}->{cur}; new ID not in GTEx v8"
    a8.append((s, note))
print(f"[C03] targets whose GRCh37 ID is in no GTEx egenes file: {len(a8)}")
for s, n in a8:
    print(f"        {s}: {n}")

gt = {}
for r in scope:
    s = r["symbol"]
    d = dict(tested=[], inst={}, id_used="", id_how="")
    if s in gtex_id:
        gid, how = gtex_id[s]
        d["id_used"], d["id_how"] = gid, how
        d["tested"] = sorted(egenes[gid])
        for t, v, p, b, se in signif.get(gid, []):
            if p >= P or t not in egenes[gid]:
                continue
            _, c, gs, ge, _ = egenes[gid][t]
            vc, vp = v.split("_")[0], int(v.split("_")[1])
            if vc == c and gs - FLANK <= vp <= ge + FLANK:
                if t not in d["inst"] or p < d["inst"][t][1]:
                    d["inst"][t] = (v, p, b, se)
    gt[s] = d

# ---- gene body lookup for deCODE (GENCODE v26 as carried in GTEx egenes) ----
def gtex_body(eid):
    if eid not in egenes:
        return None
    recs = list(egenes[eid].values())
    return recs[0][1].replace("chr", ""), recs[0][2], recs[0][3]


def best(cur, cand):
    return cand if cur is None or cand["p"] < cur["p"] else cur


# ---- UKB-PPP ---------------------------------------------------------------
h3, d3 = table(os.path.join(PQ, "ukbppp_ST3.tsv"), 2)
ukb_prot = collections.defaultdict(list)       # symbol -> [(proteinID, chr, start, end)]
ukb_pid = {}
for x in d3:
    r = dict(zip(h3, x))
    syms = [clean(a) for a in r["Gene symbol"].split(";")]
    chs = str(r["Gene CHROM"]).split(";"); sts = str(r["Gene start"]).split(";"); ens = str(r["Gene end"]).split(";")
    for i, sy in enumerate(syms):
        j = i if i < len(chs) else 0
        try:
            rec = (r["UKBPPP ProteinID"], chs[j], int(float(sts[j])), int(float(ens[j])))
        except ValueError:
            rec = (r["UKBPPP ProteinID"], chs[j], None, None)
        ukb_prot[sy].append(rec)
        ukb_pid.setdefault(r["UKBPPP ProteinID"], []).append((sy,) + rec[1:])
ukb_rows = collections.defaultdict(list)
for sheet in ("ST9", "ST10"):
    h, d = table(os.path.join(PQ, f"ukbppp_{sheet}.tsv"), 4)
    cb = [c for c in h if c.startswith("BETA")][0]; cs = [c for c in h if c.startswith("SE")][0]
    cl = [c for c in h if c.startswith("log10(p)")][0]; cf = [c for c in h if c.startswith("A1FREQ")][0]
    for x in d:
        r = dict(zip(h, x))
        if clean(r["cis/trans"]).lower() != "cis":
            continue
        # FIX after Step 3 (error in my first run): ST9/ST10 code chromosome X as 23, ST3 as X
        ukb_rows[r["UKBPPP ProteinID"]].append(dict(sheet=sheet, var=r[h[0]], chr="X" if str(r["CHROM"]) == "23" else str(r["CHROM"]), pos=int(float(r["GENPOS (hg38)"])),
                                                    logp=float(r[cl]), beta=r[cb], se=r[cs], a1freq=r[cf], rsid=r["rsID"]))

# ---- deCODE ----------------------------------------------------------------
h1, d1 = table(os.path.join(PQ, "decode_ST01.tsv"), 2)
dec_by_sym, dec_by_id, dec_seq = collections.defaultdict(list), collections.defaultdict(list), {}
for x in d1:
    r = dict(zip(h1, x))
    inc = clean(r["Included in analysis"])
    syms = [clean(a) for a in r["Gene"].split("|")]; ids = [clean(a) for a in r["Ensembl.Gene.ID"].split("|")]
    dec_seq[r["SeqId"]] = (syms, ids, inc)
    for a in syms:
        dec_by_sym[a].append((r["SeqId"], inc))
    for a in ids:
        if a.startswith("ENSG"):
            dec_by_id[a].append((r["SeqId"], inc))
h2, d2 = table(os.path.join(PQ, "decode_ST02.tsv"), 2)
dec_rows = collections.defaultdict(list)
for x in d2:
    r = dict(zip(h2, x))
    if clean(r["cis/ trans"]).lower() != "cis":
        continue
    dec_rows[r["SeqId"]].append(dict(var=r["variant"], chr=r["chr (var.)"].replace("chr", ""), pos=int(float(r["pos (var.)"])),
                                     logp=float(r["-Log10(P) (adj.)"]), gene=clean(r["gene  (prot.)"]) if "gene  (prot.)" in r else ""))

# ---- Fenland ---------------------------------------------------------------
hf3, f3 = table(os.path.join(PQ, "fenland_ST3.tsv"), 0)
fen_by_sym, fen_body = collections.defaultdict(list), {}
for x in f3:
    r = dict(zip(hf3, x))
    syms = [clean(a) for a in r["Protein encoding gene"].replace(",", "|").split("|") if clean(a)]
    chs = str(r["Chr"]).split("|"); sts = str(r["start"]).split("|"); ens = str(r["end"]).split("|")
    for i, sy in enumerate(syms):
        j = i if i < len(chs) else 0
        try:
            body = ("X" if chs[j] == "23" else chs[j], int(float(sts[j])), int(float(ens[j])))
        except (ValueError, IndexError):
            body = None
        fen_by_sym[sy].append(r["SomaScan.ID"]); fen_body[(r["SomaScan.ID"], sy)] = body
        if len(syms) == 1:
            fen_body[(r["SomaScan.ID"], None)] = body      # single-gene aptamer: body regardless of symbol spelling
# build check (clarification 4): ST3 coordinates vs Ensembl GRCh37.87
import gzip, re
g37 = {}
with gzip.open(r"D:\causal2drug\data\raw\Homo_sapiens.GRCh37.87.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] == "gene":
            g37.setdefault(re.search(r'gene_name "([^"]+)"', f[8]).group(1), []).append((f[0], int(f[3]), int(f[4])))
n37 = n_cmp = 0
for (sid, sy), body in fen_body.items():
    if sy and body and len(g37.get(sy, [])) == 1 and len(fen_by_sym[sy]) >= 1:
        n_cmp += 1; n37 += (g37[sy][0][1], g37[sy][0][2]) == (body[1], body[2])
print(f"[C04] Fenland ST3 gene start/end identical to Ensembl GRCh37.87: {n37} of {n_cmp} comparable aptamer-genes -> build GRCh37")
hf2, f2 = table(os.path.join(PQ, "fenland_ST2.tsv"), 0)
fen_rows = []
fen_ids_in_st2 = set()
for x in f2:
    r = dict(zip(hf2, x))
    fen_ids_in_st2.update(clean(a) for a in r["ENSEMBL.id.protein"].split("|"))
    if clean(r["cis.trans"]).lower() != "cis":
        continue
    try:
        p = float(r["P.value"])
    except ValueError:
        continue
    fen_rows.append(dict(sid=r["SomaScan.id"], syms=[clean(a) for a in r["HGNC.symbol.protein"].split("|")],
                         ids=[clean(a) for a in r["ENSEMBL.id.protein"].split("|")], var=r["MarkerName"], rsid=r["rsID"],
                         chr="X" if str(r["Chr"]) == "23" else str(r["Chr"]), pos=int(float(r["Position"])), p=p,
                         ea=r["EA"], nea=r["NEA"], eaf=r["EAF"], beta=r["Effect"], se=r["SE"]))

# ---- per target ------------------------------------------------------------
rows = []
for r in scope:
    s, eid = r["symbol"], r["ensembl_id"]
    g = gt[s]
    nonblood = sorted(t for t in g["inst"] if t != BLOOD)
    cv = [t for t in CV_TISSUES if t in g["inst"]]
    best_cv = min(cv, key=lambda t: g["inst"][t][1]) if cv else ""

    # UKB-PPP (symbol only)
    u_assayed = s in ukb_prot; u_best = None; u_unver = False
    for pid, c, gs, ge in ukb_prot.get(s, []):
        for q in ukb_rows.get(pid, []):
            if q["logp"] <= LOGP:
                continue
            if gs is None:
                u_unver = True; continue
            if q["chr"] == c and gs - FLANK <= q["pos"] <= ge + FLANK:
                u_best = best(u_best, dict(q, p=-q["logp"]))
    # deCODE (symbol or Ensembl ID); gene body from GTEx egenes
    seqs = {a for a, _ in dec_by_sym.get(s, [])} | {a for a, _ in dec_by_id.get(eid, [])}
    d_listed = bool(seqs); d_included = any(dec_seq[a][2] == "Yes" for a in seqs)
    body = gtex_body(g["id_used"]) if g["id_used"] else None
    d_best = None; d_unver = False
    for a in seqs:
        for q in dec_rows.get(a, []):
            if q["logp"] <= LOGP:
                continue
            if body is None:
                d_unver = True; continue
            if q["chr"] == body[0] and body[1] - FLANK <= q["pos"] <= body[2] + FLANK:
                d_best = best(d_best, dict(q, p=-q["logp"]))
    # Fenland (symbol or Ensembl ID); gene body from ST3
    f_assayed = s in fen_by_sym or eid in fen_ids_in_st2; f_best = None; f_unver = False
    for q in fen_rows:
        if s not in q["syms"] and eid not in q["ids"]:
            continue
        if q["p"] >= P:
            continue
        # FIX after Step 3 (error in my first run): a row matched by Ensembl ID whose ST3 symbol is an older
        # name (VEGFD is FIGF in ST3) still has an ST3 gene body; take it from the single-gene aptamer
        b = fen_body.get((q["sid"], s)) or fen_body.get((q["sid"], None))
        if b is None:
            f_unver = True; continue
        if q["chr"] == b[0] and b[1] - FLANK <= q["pos"] <= b[2] + FLANK:
            f_best = best(f_best, q)

    pq = bool(u_best or d_best or f_best)
    cat = ("fixable-both" if nonblood and pq else "fixable-tissue" if nonblood else "fixable-protein" if pq else "locked")
    assayed_any = u_assayed or d_listed or f_assayed
    assayed_any_strict = u_assayed or d_included or f_assayed
    if cv or pq:
        tier = "plausible rescue"
    elif nonblood:
        tier = "off-target-tissue rescue"
    else:
        tier = "locked: assayed, no instrument" if assayed_any else "locked: never assayed"
    rows.append(dict(
        symbol=s, ensembl_id=eid, chrom=r["chrom"], phase=r["phase"], eqtlgen_status=r["status"], subtype=r["subtype"],
        x_flag=r["chrom"] == "X", gtex_id_used=g["id_used"], gtex_id_how=g["id_how"], gtex_n_tissues_tested=len(g["tested"]),
        gtex_nonblood_tissues=";".join(nonblood), gtex_n_nonblood=len(nonblood), gtex_cv_tissues=";".join(cv),
        gtex_whole_blood=BLOOD in g["inst"], gtex_best_cv_tissue=best_cv,
        gtex_best_cv_variant=g["inst"][best_cv][0] if best_cv else "", gtex_best_cv_p=g["inst"][best_cv][1] if best_cv else "",
        gtex_best_cv_slope=g["inst"][best_cv][2] if best_cv else "", gtex_best_cv_slope_se=g["inst"][best_cv][3] if best_cv else "",
        ukb_assayed=u_assayed, ukb_instrument=bool(u_best), ukb_lead=u_best["var"] if u_best else "",
        ukb_sheet=u_best["sheet"] if u_best else "", ukb_log10p=u_best["logp"] if u_best else "",
        ukb_beta=u_best["beta"] if u_best else "", ukb_se=u_best["se"] if u_best else "", ukb_pos38=u_best["pos"] if u_best else "",
        decode_listed=d_listed, decode_included=d_included, decode_instrument=bool(d_best),
        decode_lead=d_best["var"] if d_best else "", decode_log10p=d_best["logp"] if d_best else "",
        fenland_assayed=f_assayed, fenland_instrument=bool(f_best), fenland_lead=f_best["var"] if f_best else "",
        fenland_p=f_best["p"] if f_best else "",
        unverified="|".join(k for k, v in (("ukb", u_unver), ("decode", d_unver), ("fenland", f_unver)) if v),
        assayed_any=assayed_any, assayed_any_decode_included_only=assayed_any_strict, category=cat, tier=tier))

with open(os.path.join(OUT, "audit_recovery_map.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

C = collections.Counter
print(f"[C05] categories: {dict(sorted(C(r['category'] for r in rows).items()))}")
print(f"[C06] tiers: {dict(sorted(C(r['tier'] for r in rows).items()))}")
cvp = sum(1 for r in rows if r["gtex_cv_tissues"]); pq_ = sum(1 for r in rows if r["ukb_instrument"] or r["decode_instrument"] or r["fenland_instrument"])
print(f"[C07] CV-tissue instrument {cvp}; plasma cis-pQTL instrument {pq_} (UKB {sum(r['ukb_instrument'] for r in rows)}, "
      f"deCODE {sum(r['decode_instrument'] for r in rows)}, Fenland {sum(r['fenland_instrument'] for r in rows)})")
print(f"[C08] locked under 'assayed = deCODE Included==Yes only': "
      f"{dict(C(('assayed' if r['assayed_any_decode_included_only'] else 'never') for r in rows if r['category']=='locked'))}")
print(f"[C09] tier x eQTLGen status: {dict(sorted(C((r['tier'], r['eqtlgen_status'], r['subtype']) for r in rows).items()))}")
print(f"[C10] GTEx whole-blood-only (locked, flagged): {[r['symbol'] for r in rows if r['category']=='locked' and r['gtex_whole_blood']]}")
print(f"[C11] rows with an UNVERIFIED cell: {[(r['symbol'], r['unverified']) for r in rows if r['unverified']]}")
print(f"[C12] locked: {sorted(r['symbol'] for r in rows if r['category']=='locked')}")
