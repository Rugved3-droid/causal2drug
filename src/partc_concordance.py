#!/usr/bin/env python
"""PART C: rule-based distributional validation (PREREGISTRATION R13).

For every qualifying target-indication pair, run the frozen cis-MR (R12.2 lead-SNP Wald
ratio) and compare the sign of the estimate against the direction implied by an approved
drug's mechanism.

Caches the tabix index per outcome GWAS (155 distinct files) and checkpoints per pair, so
an interruption costs at most one pair.
"""
import collections
import concurrent.futures as cf
import csv
import threading
import gzip
import json
import math
import os
import sys
import time

import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_http as T                 # noqa: E402
from harmonize import harmonize_one    # noqa: E402
from pyliftover import LiftOver        # noqa: E402

D1B = r"D:\causal2drug\deliverables\d1b"
OUT = r"D:\causal2drug\deliverables\d1f"
RAW = r"D:\causal2drug\data\raw"
CKPT = OUT + r"\ckpt"
GW = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
FLANK = 100_000
COLS = ["gene", "mesh_id", "mesh_term", "gwas_acc", "drug_direction", "n_independent",
        "lead_rsid", "lead_p_exp", "beta_exp", "se_exp", "F",
        "beta_out_harm", "se_out", "harmonization", "n_overlap",
        "mr_beta", "mr_se", "mr_p", "mr_direction", "expected_direction", "concordant",
        "note"]


def smr(z, af, n):
    v = 2.0 * af * (1.0 - af) * (n + z * z)
    if v <= 0:
        return None, None
    se = 1.0 / math.sqrt(v)
    return z * se, se


def load_af(rsids):
    af = {}
    with gzip.open(RAW + r"\eqtlgen_SNP_AF.txt.gz", "rt", encoding="utf-8") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_s, i_f, i_a, i_b = (h.index("SNP"), h.index("AlleleB_all"),
                              h.index("AlleleA"), h.index("AlleleB"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_s] in rsids:
                try:
                    af[f[i_s]] = (float(f[i_f]), f[i_a], f[i_b])
                except ValueError:
                    pass
    return af


def rel_path(acc, harm_index):
    return harm_index[acc]


def main():
    os.makedirs(CKPT, exist_ok=True)
    pairs = list(csv.DictReader(open(OUT + r"\partc_pairs.tsv", encoding="utf-8"),
                                delimiter="\t"))
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D1B + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    genes = {p["gene"] for p in pairs}

    # harmonised path per accession
    import re
    harm_index = {}
    for l in open(RAW + r"\harmonised_list.txt"):
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m:
            harm_index[m.group(1)] = l.strip().lstrip("./")

    # exposure: eQTLGen hits inside +-100 kb, per target
    hits = collections.defaultdict(list)
    with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = r["MinikelSymbol"]
            if s not in genes:
                continue
            c = coords[s]
            if r["SNPChr"] != c["chrom"]:
                continue
            if int(c["start"]) - FLANK <= int(r["SNPPos"]) <= int(c["end"]) + FLANK:
                hits[s].append(r)
    need = {r["SNP"] for v in hits.values() for r in v}
    print("targets %d ; exposure variants %d ; looking up allele frequencies..."
          % (len(hits), len(need)), flush=True)
    af = load_af(need)
    print("  AF found for %d/%d" % (len(af), len(need)), flush=True)

    lo = LiftOver(RAW + r"\hg19ToHg38.over.chain.gz")
    exp = {}
    for s, rows in hits.items():
        out = []
        for r in rows:
            a = af.get(r["SNP"])
            if not a:
                continue
            afb, alA, alB = a
            if r["AssessedAllele"] == alB:
                fa = afb
            elif r["AssessedAllele"] == alA:
                fa = 1.0 - afb
            else:
                continue
            z = float(r["Zscore"])
            b, se = smr(z, fa, float(r["NrSamples"]))
            if b is None:
                continue
            res = lo.convert_coordinate("chr" + r["SNPChr"], int(r["SNPPos"]) - 1)
            if not res:
                continue
            out.append(dict(rsid=r["SNP"], chrom=res[0][0].replace("chr", ""),
                            pos38=res[0][1] + 1, assessed=r["AssessedAllele"],
                            other=r["OtherAllele"], af=fa, beta=b, se=se,
                            p=float(r["Pvalue"]), F=z * z))
        if out:
            exp[s] = out
    print("targets with usable exposure after AF + liftover: %d" % len(exp), flush=True)

    by_acc = collections.defaultdict(list)
    for p in pairs:
        by_acc[p["gwas_acc"]].append(p)

    done = {}
    for f in os.listdir(CKPT):
        for r in csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                delimiter="\t"):
            done[(r["gene"], r["mesh_id"], r["gwas_acc"])] = r
    print("resuming: %d pairs already done" % len(done), flush=True)

    skipped = {}
    t0 = time.time()
    for acc in sorted(by_acc):
        cpath = os.path.join(CKPT, "%s.tsv" % acc)
        todo = [p for p in by_acc[acc]
                if (p["gene"], p["mesh_id"], acc) not in done]
        if not todo:
            continue
        if acc not in harm_index:
            continue
        url = GW + harm_index[acc]
        sess = requests.Session()
        try:
            lines, _ = T.header_lines(url, sess=sess, nbytes=1 << 16)
            gcols = lines[0].rstrip("\r").split("\t")
            gi = {c: i for i, c in enumerate(gcols)}
            base_ok = all(k in gi for k in ("chromosome", "base_pair_location",
                                             "effect_allele", "other_allele",
                                             "standard_error"))
            # case-control harmonised files report odds_ratio instead of beta; the
            # accompanying standard_error is already on the log-OR scale, so
            # beta = ln(OR) is the standard conversion, not a modelling choice.
            use_or = ("beta" not in gi) and ("odds_ratio" in gi)
            if not base_ok or not ("beta" in gi or use_or):
                # distinguish the two causes: a plain-gzip file yields an EMPTY header
                # here, which would otherwise be mislabelled as a missing-column problem
                probe = sess.get(url, headers={"Range": "bytes=0-17"}, timeout=60).content
                is_bgzf, _ = T.is_bgzf(probe)
                has_tbi = sess.head(url + ".tbi", timeout=60).status_code == 200
                skipped[acc] = ("plain_gzip_no_tbi_not_range_queryable"
                                if (not is_bgzf or not has_tbi)
                                else "bgzf_indexed_but_lacks_beta_and_odds_ratio")
                continue
            ri = sess.get(url + ".tbi", timeout=600)
            ri.raise_for_status()
            tbi = T.Tbi(ri.content)
        except Exception as e:
            skipped[acc] = "not_bgzf_or_no_tbi:%s" % type(e).__name__
            continue

        fresh = not os.path.exists(cpath)
        fh = open(cpath, "a", newline="", encoding="utf-8")
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        if fresh:
            w.writeheader()
        lock = threading.Lock()

        def do_pair(p):
            s = p["gene"]
            rec = {c: "" for c in COLS}
            rec.update(dict(gene=s, mesh_id=p["mesh_id"], mesh_term=p["mesh_term"],
                            gwas_acc=acc, drug_direction=p["drug_direction"],
                            n_independent=p["n_independent"]))
            ev = exp.get(s)
            if not ev:
                rec["note"] = "no_usable_exposure"
                return rec
            chrom = ev[0]["chrom"]
            lo_p = min(e["pos38"] for e in ev) - 1000
            hi_p = max(e["pos38"] for e in ev) + 1000
            try:
                use = chrom if chrom in tbi.names else ("chr" + chrom)
                ch = tbi.chunks(use, lo_p - 1, hi_p) if use in tbi.names else []
                if not ch:
                    rec["note"] = "no_gwas_chunks"
                    return rec
                grows, _st = T.query(url, chrom, lo_p, hi_p, sess=sess)
            except Exception as e:
                rec["note"] = "gwas_query_failed:%s" % type(e).__name__
                return rec
            gmap = collections.defaultdict(list)
            for g in grows:
                try:
                    gmap[int(g[gi["base_pair_location"]])].append(g)
                except (ValueError, IndexError):
                    pass
            merged = []
            for e in ev:
                for g in gmap.get(e["pos38"], []):
                    try:
                        so = float(g[gi["standard_error"]])
                        if use_or:
                            orv = float(g[gi["odds_ratio"]])
                            if orv <= 0:
                                continue
                            bo = math.log(orv)
                        else:
                            bo = float(g[gi["beta"]])
                    except ValueError:
                        continue
                    try:
                        afo = float(g[gi["effect_allele_frequency"]])
                    except (ValueError, KeyError):
                        afo = None
                    h = harmonize_one(e["other"], e["assessed"], g[gi["effect_allele"]],
                                      g[gi["other_allele"]], bo,
                                      af_alt_exp=e["af"], af_effect_out=afo)
                    if h["sign"] is None or so <= 0:
                        continue
                    merged.append((e, h["beta_out_harmonized"], so, h["status"]))
                    break
            rec["n_overlap"] = len(merged)
            if not merged:
                rec["note"] = "no_overlap_after_harmonization"
                return rec
            e, bo, so, hstat = min(merged, key=lambda m: m[0]["p"])
            beta = bo / e["beta"]
            se = abs(so / e["beta"])
            z = beta / se
            pval = math.erfc(abs(z) / math.sqrt(2.0))
            exp_dir = "POSITIVE" if p["drug_direction"] == "DOWN" else "NEGATIVE"
            got = "POSITIVE" if beta > 0 else "NEGATIVE"
            rec.update(dict(lead_rsid=e["rsid"], lead_p_exp="%.6g" % e["p"],
                            beta_exp="%.6g" % e["beta"], se_exp="%.6g" % e["se"],
                            F="%.4g" % e["F"], beta_out_harm="%.6g" % bo,
                            se_out="%.6g" % so, harmonization=hstat,
                            mr_beta="%.6g" % beta, mr_se="%.6g" % se,
                            mr_p="%.6g" % pval, mr_direction=got,
                            expected_direction=exp_dir,
                            concordant="YES" if got == exp_dir else "NO"))
            return rec

        with cf.ThreadPoolExecutor(max_workers=6) as ex:
            for rec in ex.map(do_pair, todo):
                with lock:
                    w.writerow(rec); fh.flush()
        os.fsync(fh.fileno())
        fh.close()
        n_done = sum(1 for f in os.listdir(CKPT)
                     for _ in csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                             delimiter="\t"))
        print("  %-14s +%-4d  cum=%-5d  %.0fs" % (acc, len(todo), n_done, time.time() - t0),
              flush=True)

    # merge
    allrows = []
    for f in sorted(os.listdir(CKPT)):
        allrows += list(csv.DictReader(open(os.path.join(CKPT, f), encoding="utf-8"),
                                       delimiter="\t"))
    with open(OUT + r"\partc_results.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(allrows)
    print("wrote %d rows to partc_results.tsv" % len(allrows))
    npair = collections.Counter(p["gwas_acc"] for p in pairs)
    with open(OUT + chr(92) + "partc_skipped_accessions.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w2 = csv.writer(fh, delimiter=chr(9), lineterminator=chr(10))
        w2.writerow(["gwas_acc", "reason", "n_pairs_lost"])
        for a, why in sorted(skipped.items()):
            w2.writerow([a, why, npair[a]])
    print("skipped accessions: %d  (pairs lost: %d)"
          % (len(skipped), sum(npair[a] for a in skipped)))


if __name__ == "__main__":
    main()
