#!/usr/bin/env python
"""Deliverable 1 -- positive-control chain, run strictly under
deliverables/d1/PREREGISTRATION.md. No parameter in this file may be changed in
response to a result.

Per pair: retrieve trait GWAS cis-region -> retrieve matched cis-eQTL region ->
verify genome build -> harmonize alleles (incl. explicit strand-ambiguous handling)
-> retrieve 1000G EUR-unrelated genotypes by range query -> clump -> MR -> coloc
-> write provenance. Every intermediate is written unmodified.
"""
import csv
import datetime
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T  # noqa: E402
import requests  # noqa: E402

OUT = r"D:\causal2drug\deliverables\d1"
PED = r"D:\causal2drug\data\raw\20130606_g1k_3202_samples_ped_population.txt"

# ---- pre-registered constants (R1, R5, R8) ----
FLANK = 100_000
P_INSTRUMENT = 5e-8
R2_CLUMP = 0.01
F_MIN = 10.0
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_QUANT = 0.15 ** 2

GWAS_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
EQTL_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
            "QTS000015/QTD000266/QTD000266.all.tsv.gz")
KG_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
          "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
          "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

GWAS_COLS = ("chromosome base_pair_location effect_allele other_allele beta standard_error "
             "effect_allele_frequency p_value variant_id n N_studies hm_coordinate_conversion "
             "hm_code rsid").split()
EQTL_COLS = ("molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue "
             "beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid").split()

PAIRS = [
    dict(name="HMGCR_LDL", gene="HMGCR", ensg="ENSG00000113161",
         trait="LDL cholesterol", efo="EFO_0004611", gwas_acc="GCST90239658",
         gwas_rel="./GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         expect="POSITIVE"),
    dict(name="PCSK9_LDL", gene="PCSK9", ensg="ENSG00000169174",
         trait="LDL cholesterol", efo="EFO_0004611", gwas_acc="GCST90239658",
         gwas_rel="./GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
         expect="POSITIVE"),
    dict(name="IL6R_CRP", gene="IL6R", ensg="ENSG00000160712",
         trait="C-reactive protein", efo="EFO_0004458", gwas_acc="GCST90179146",
         gwas_rel="./GCST90179001-GCST90180000/GCST90179146/harmonised/GCST90179146.h.tsv.gz",
         expect="NEGATIVE"),
]

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def rc(a):
    return "".join(COMP.get(b, "N") for b in reversed(a))


def is_snp(a):
    return len(a) == 1 and a in COMP


def write_tsv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def gwas_header(url):
    """Harmonised GWAS Catalog files do NOT share one column layout (observed: 14 cols
    with n/N_studies, 13 cols with z). Read the header from the file itself."""
    lines, st = T.header_lines(url, nbytes=1 << 16)
    cols = lines[0].rstrip(chr(13)).split(chr(9))
    need = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "p_value"]
    missing = [c for c in need if c not in cols]
    if missing:
        raise RuntimeError("GWAS file missing required columns %r: %s" % (missing, url))
    return cols, st


# ------------------------- gene coordinates (R1) -------------------------
def ensembl_gene(symbol):
    r = requests.get("https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s" % symbol,
                     headers={"Content-Type": "application/json"}, timeout=120)
    r.raise_for_status()
    d = r.json()
    return {"symbol": d["display_name"], "ensembl_id": d["id"], "chrom": str(d["seq_region_name"]),
            "start": int(d["start"]), "end": int(d["end"]), "strand": d["strand"],
            "assembly": d["assembly_name"], "retrieved_utc": now(),
            "source": "https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s" % symbol}


# ------------------------- harmonization (R6) -------------------------
def harmonize(gwas_rows, eqtl_rows, gwas_cols):
    gi = {c: i for i, c in enumerate(gwas_cols)}
    qi = {c: i for i, c in enumerate(EQTL_COLS)}
    gmap = {}
    for r in gwas_rows:
        gmap.setdefault(int(r[gi["base_pair_location"]]), []).append(r)

    counts = dict(eqtl_variants=len(eqtl_rows), gwas_variants=len(gwas_rows),
                  position_matched=0, aligned_same=0, aligned_flipped=0,
                  strand_flipped=0, allele_mismatch=0, indel_or_multiallelic=0,
                  ambiguous_total=0, ambiguous_resolved=0, ambiguous_unresolved=0,
                  missing_stats=0, harmonized=0)
    out = []
    for q in eqtl_rows:
        pos = int(q[qi["position"]])
        ref, alt = q[qi["ref"]], q[qi["alt"]]
        cands = gmap.get(pos)
        if not cands:
            continue
        counts["position_matched"] += 1
        if not (is_snp(ref) and is_snp(alt)):
            counts["indel_or_multiallelic"] += 1
            continue
        resolved = None
        for g in cands:
            ea, oa = g[gi["effect_allele"]], g[gi["other_allele"]]
            if not (is_snp(ea) and is_snp(oa)):
                continue
            if (ea, oa) == (alt, ref):
                resolved = (g, 1.0, "aligned_same"); break
            if (ea, oa) == (ref, alt):
                resolved = (g, -1.0, "aligned_flipped"); break
            rea, roa = rc(ea), rc(oa)
            if (rea, roa) == (alt, ref):
                resolved = (g, 1.0, "strand_flipped"); break
            if (rea, roa) == (ref, alt):
                resolved = (g, -1.0, "strand_flipped"); break
        if resolved is None:
            if all(not (is_snp(g[gi["effect_allele"]]) and is_snp(g[gi["other_allele"]]))
                   for g in cands):
                counts["indel_or_multiallelic"] += 1
            else:
                counts["allele_mismatch"] += 1
            continue
        g, sign, how = resolved

        # eQTL alt-allele frequency from ac/an
        try:
            af_alt_eqtl = float(q[qi["ac"]]) / float(q[qi["an"]])
        except (ValueError, ZeroDivisionError):
            af_alt_eqtl = float("nan")
        try:
            af_eff_gwas = float(g[gi["effect_allele_frequency"]])
        except ValueError:
            af_eff_gwas = float("nan")
        # frequency of the ALT allele as implied by the GWAS, after alignment
        af_alt_gwas = af_eff_gwas if sign > 0 else (1.0 - af_eff_gwas)

        ambiguous = (ref, alt) in AMBIG
        amb_note = ""
        if ambiguous:
            counts["ambiguous_total"] += 1
            ok = (not math.isnan(af_alt_eqtl) and not math.isnan(af_alt_gwas)
                  and (af_alt_eqtl < 0.42 or af_alt_eqtl > 0.58)
                  and (af_alt_gwas < 0.42 or af_alt_gwas > 0.58)
                  and ((af_alt_eqtl < 0.5) == (af_alt_gwas < 0.5)))
            if ok:
                counts["ambiguous_resolved"] += 1
                amb_note = "resolved_by_freq"
            else:
                counts["ambiguous_unresolved"] += 1
                continue

        try:
            b_g = float(g[gi["beta"]]); se_g = float(g[gi["standard_error"]])
            p_g = float(g[gi["p_value"]])
            b_q = float(q[qi["beta"]]); se_q = float(q[qi["se"]]); p_q = float(q[qi["pvalue"]])
        except ValueError:
            counts["missing_stats"] += 1
            continue
        if se_g <= 0 or se_q <= 0:
            counts["missing_stats"] += 1
            continue

        counts[how] = counts.get(how, 0) + 1
        counts["harmonized"] += 1
        out.append(dict(
            chrom=q[qi["chromosome"]], pos=pos, variant_eqtl=q[qi["variant"]],
            rsid_eqtl=q[qi["rsid"]], rsid_gwas=(g[gi["rsid"]] if "rsid" in gi else "NA"),
            eqtl_ref=ref, eqtl_alt=alt,
            gwas_effect_allele=g[gi["effect_allele"]], gwas_other_allele=g[gi["other_allele"]],
            harmonization=how, sign_applied=sign,
            ambiguous="TRUE" if ambiguous else "FALSE", ambiguous_note=amb_note,
            effect_allele_harmonized=alt, other_allele_harmonized=ref,
            beta_eqtl=b_q, se_eqtl=se_q, p_eqtl=p_q,
            beta_gwas_raw=b_g, beta_gwas_harmonized=sign * b_g, se_gwas=se_g, p_gwas=p_g,
            af_alt_eqtl=af_alt_eqtl, af_alt_gwas=af_alt_gwas,
            n_gwas=(g[gi["n"]] if "n" in gi else "NA"),
            an_eqtl=q[qi["an"]], maf_eqtl=q[qi["maf"]],
            F_eqtl=(b_q / se_q) ** 2))
    return out, counts


# ------------------------- 1000G LD (R5) -------------------------
def eur_unrelated_ids():
    ids = []
    with open(PED) as fh:
        head = fh.readline().split()
        i_s, i_f, i_m, i_sp = (head.index("SampleID"), head.index("FatherID"),
                               head.index("MotherID"), head.index("Superpopulation"))
        for l in fh:
            f = l.split()
            if f[i_sp] == "EUR" and f[i_f] == "0" and f[i_m] == "0":
                ids.append(f[i_s])
    return ids


def kg_dosages(chrom, start, end, keep_ids):
    url = KG_URL % chrom
    hdr_lines, hstats = T.header_lines(url, nbytes=1 << 18)
    chrom_line = next(l for l in hdr_lines if l.startswith("#CHROM"))
    samples = chrom_line.rstrip("\r").split("\t")[9:]
    idx = [i for i, s in enumerate(samples) if s in set(keep_ids)]
    rows, stats = T.query(url, "chr" + str(chrom), start, end)
    stats["bytes"] = stats.get("bytes", 0) + hstats.get("bytes", 0)
    dos = {}
    for r in rows:
        ref, alt = r[3], r[4]
        if not (is_snp(ref) and is_snp(alt)):
            continue
        gts = r[9:]
        v = []
        for i in idx:
            g = gts[i]
            v.append((g[0] == "1") + (g[2] == "1") if len(g) >= 3 else 0)
        dos[(int(r[1]), ref, alt)] = v
    return dos, stats, len(idx)


def r2(a, b):
    n = len(a)
    ma = sum(a) / n
    mb = sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va == 0 or vb == 0:
        return 0.0
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    return (cov * cov) / (va * vb)


def clump(cands, dos):
    """cands sorted by eQTL p ascending. Greedy, drop r2 >= R2_CLUMP."""
    kept, dropped = [], []
    for c in cands:
        k = (c["pos"], c["eqtl_ref"], c["eqtl_alt"])
        gk = dos.get(k)
        conflict = None
        for s in kept:
            ks = (s["pos"], s["eqtl_ref"], s["eqtl_alt"])
            gs = dos.get(ks)
            if gk is None or gs is None:
                continue
            v = r2(gk, gs)
            if v >= R2_CLUMP:
                conflict = (s["variant_eqtl"], v)
                break
        if conflict:
            dropped.append((c["variant_eqtl"], conflict[0], conflict[1]))
        else:
            kept.append(c)
    return kept, dropped


# ------------------------- MR (R7) -------------------------
def norm_sf(z):
    return 0.5 * math.erfc(abs(z) / math.sqrt(2.0)) * 2.0


def mr(insts):
    if len(insts) == 1:
        i = insts[0]
        b = i["beta_gwas_harmonized"] / i["beta_eqtl"]
        se = abs(i["se_gwas"] / i["beta_eqtl"])
        return dict(method="Wald ratio", n_instruments=1, beta=b, se=se,
                    z=b / se, p=norm_sf(b / se))
    num = sum(i["beta_gwas_harmonized"] * i["beta_eqtl"] / i["se_gwas"] ** 2 for i in insts)
    den = sum(i["beta_eqtl"] ** 2 / i["se_gwas"] ** 2 for i in insts)
    b = num / den
    se = math.sqrt(1.0 / den)
    return dict(method="IVW fixed effects", n_instruments=len(insts), beta=b, se=se,
                z=b / se, p=norm_sf(b / se))


# ------------------------- coloc (R8) -------------------------
def labf(beta, se, W=W_QUANT):
    V = se * se
    r = W / (W + V)
    z = beta / se
    return 0.5 * (math.log(1 - r) + r * z * z)


def lse(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def ldiff(a, b):
    return a + math.log1p(-math.exp(b - a)) if b < a else float("-inf")


def coloc(rows):
    l1 = [labf(r["beta_eqtl"], r["se_eqtl"]) for r in rows]
    l2 = [labf(r["beta_gwas_harmonized"], r["se_gwas"]) for r in rows]
    s1, s2 = lse(l1), lse(l2)
    s12 = lse([a + b for a, b in zip(l1, l2)])
    lH0 = 0.0
    lH1 = math.log(P1) + s1
    lH2 = math.log(P2) + s2
    lH3 = math.log(P1) + math.log(P2) + ldiff(s1 + s2, s12)
    lH4 = math.log(P12) + s12
    all_l = [lH0, lH1, lH2, lH3, lH4]
    d = lse(all_l)
    pp = [math.exp(x - d) for x in all_l]
    return dict(n_variants=len(rows), PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4])


# ------------------------- driver -------------------------
def run_pair(p, gene_cache):
    log = {"pair": p["name"], "started_utc": now(), "preregistration":
           "deliverables/d1/PREREGISTRATION.md", "bytes_transferred": 0}
    g = gene_cache.setdefault(p["gene"], ensembl_gene(p["gene"]))
    log["gene"] = g
    if g["assembly"] != "GRCh38":
        raise RuntimeError("gene coords not GRCh38: %s" % g["assembly"])
    chrom, s, e = g["chrom"], g["start"] - FLANK, g["end"] + FLANK
    log["window"] = {"chrom": chrom, "start": s, "end": e, "flank_bp": FLANK}

    # --- GWAS region ---
    t0 = time.time()
    gurl = GWAS_BASE + p["gwas_rel"].lstrip("./")
    gcols, hstats = gwas_header(gurl)
    grows, gstats = T.query(gurl, chrom, s, e)
    gstats["bytes"] += hstats.get("bytes", 0)
    gstats["columns"] = gcols
    gstats["retrieved_utc"] = now(); gstats["seconds"] = round(time.time() - t0, 2)
    gstats["accession"] = p["gwas_acc"]
    log["gwas_retrieval"] = gstats
    log["bytes_transferred"] += gstats["bytes"]
    write_tsv(os.path.join(OUT, "extracts", "%s_GWAS_%s_%s_%d_%d.tsv"
                           % (p["name"], p["gwas_acc"], chrom, s, e)), gcols, grows)

    # --- eQTL region ---
    t0 = time.time()
    qall, qstats = T.query(EQTL_URL, chrom, s, e)
    qrows = [r for r in qall if r[16] == p["ensg"]]
    qstats["retrieved_utc"] = now(); qstats["seconds"] = round(time.time() - t0, 2)
    qstats["dataset"] = "QTD000266 GTEx_v10 liver ge (eQTL Catalogue R8)"
    qstats["rows_all_traits"] = len(qall); qstats["rows_this_gene"] = len(qrows)
    log["eqtl_retrieval"] = qstats
    log["bytes_transferred"] += qstats["bytes"]
    write_tsv(os.path.join(OUT, "extracts", "%s_EQTL_QTD000266_%s_%s_%d_%d.tsv"
                           % (p["name"], p["ensg"], chrom, s, e)), EQTL_COLS, qrows)

    # --- build check ---
    ok = 0
    for r in qrows[:200]:
        v = r[5].split("_")
        if len(v) >= 4 and v[1] == r[2] and v[2] == r[3] and v[3] == r[4]:
            ok += 1
    log["build_check"] = {
        "eqtl_variant_id_matches_pos_ref_alt": "%d/%d" % (ok, min(200, len(qrows))),
        "gwas_file_is_harmonised_GRCh38": p["gwas_rel"].endswith(".h.tsv.gz"),
        "gene_assembly": g["assembly"],
        "kg_panel_build": "GRCh38 (1kGP high coverage phased panel)"}

    # --- harmonize ---
    merged, counts = harmonize(grows, qrows, gcols)
    log["harmonization"] = counts
    mcols = list(merged[0].keys()) if merged else []
    write_tsv(os.path.join(OUT, "harmonized", "%s_harmonized_merge.tsv" % p["name"]),
              mcols, [[m[c] for c in mcols] for m in merged])

    # --- LD panel ---
    t0 = time.time()
    dos, kstats, n_eur = kg_dosages(chrom, s, e, gene_cache["_eur"])
    kstats["retrieved_utc"] = now(); kstats["seconds"] = round(time.time() - t0, 2)
    kstats["n_eur_unrelated_samples"] = n_eur
    kstats["n_biallelic_snps_in_panel_window"] = len(dos)
    log["ld_retrieval"] = kstats
    log["bytes_transferred"] += kstats["bytes"]

    # --- instruments / clumping (R5) ---
    cands = sorted([m for m in merged if m["p_eqtl"] < P_INSTRUMENT],
                   key=lambda m: m["p_eqtl"])
    kept, dropped = clump(cands, dos)
    best = min(merged, key=lambda m: m["p_eqtl"]) if merged else None
    log["instruments"] = {
        "threshold": P_INSTRUMENT, "n_candidates_before_clumping": len(cands),
        "n_after_clumping": len(kept), "clump_r2": R2_CLUMP,
        "n_dropped_by_ld": len(dropped),
        "min_eqtl_p_in_window": best["p_eqtl"] if best else None,
        "min_eqtl_p_variant": best["variant_eqtl"] if best else None,
        "min_eqtl_p_rsid": best["rsid_eqtl"] if best else None,
        "F_of_best_variant": round(best["F_eqtl"], 3) if best else None}

    # --- MR (R7) ---
    if kept:
        Fs = [i["F_eqtl"] for i in kept]
        log["instrument_strength"] = {"F_min": min(Fs), "F_mean": sum(Fs) / len(Fs),
                                      "F_criterion": F_MIN, "passes": min(Fs) >= F_MIN}
        log["mr"] = mr(kept)
        log["mr"]["direction"] = "POSITIVE" if log["mr"]["beta"] > 0 else "NEGATIVE"
        log["mr"]["expected_direction"] = p["expect"]
    else:
        log["instrument_strength"] = None
        log["mr"] = {"status": "NOT COMPUTABLE - zero instruments at p<%g" % P_INSTRUMENT}

    # --- coloc (R8) ---
    log["coloc"] = coloc(merged) if merged else {"status": "no harmonized variants"}

    # --- pass criteria (R10) ---
    log["pass_criteria"] = {
        "1_eGene_FDR005": p.get("egene"),
        "2_ge1_instrument": len(kept) >= 1,
        "3_F_ge_10": bool(log["instrument_strength"] and log["instrument_strength"]["passes"]),
        "4_MR_p_lt_005": bool(log["mr"].get("p") is not None and log["mr"]["p"] < 0.05),
        "5_direction_matches": (log["mr"].get("direction") == p["expect"]
                                if log["mr"].get("direction") else False)}
    log["VERDICT"] = "PASS" if all(log["pass_criteria"].values()) else "FAIL"
    log["finished_utc"] = now()
    with open(os.path.join(OUT, "provenance", "%s_provenance.json" % p["name"]), "w") as fh:
        json.dump(log, fh, indent=2)
    return log


def main():
    for d in ("extracts", "harmonized", "provenance"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    cache = {"_eur": eur_unrelated_ids()}
    print("EUR unrelated samples in 1000G panel:", len(cache["_eur"]))
    egene = json.load(open(os.path.join(OUT, "egene_status.json")))
    summary = []
    for p in PAIRS:
        p["egene"] = egene[p["ensg"]]["is_egene_fdr005"]
        print("\n" + "=" * 70 + "\n" + p["name"] + "\n" + "=" * 70)
        log = run_pair(p, cache)
        print(json.dumps({k: log[k] for k in
                          ("harmonization", "instruments", "mr", "coloc",
                           "pass_criteria", "VERDICT")}, indent=2)[:2600])
        summary.append(log)
    json.dump(summary, open(os.path.join(OUT, "d1_summary.json"), "w"), indent=2)
    print("\nwrote", os.path.join(OUT, "d1_summary.json"))


if __name__ == "__main__":
    main()
