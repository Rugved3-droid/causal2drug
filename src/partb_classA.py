#!/usr/bin/env python
"""PART B: is CLASS A (absent from the exposure resource) enriched for tissue-restricted
and secreted proteins?

Reference: Human Protein Atlas `proteinatlas.tsv` -- chosen over GTEx-derived tau because
it supplies BOTH a curated tissue-specificity class AND a curated secretome annotation
from one consistent source, so the two hypotheses are tested against the same reference
rather than two ad-hoc thresholds. GTEx v10 median TPM is used as an independent
cross-check on the one quantity that is most mechanistically direct: whole-blood
expression, which is what eQTLGen can actually see.

States are recomputed from the Deliverable 1b artefacts, not carried from its summary.
"""
import csv
import gzip
import json
import os
import sys
import zipfile

from scipy import stats

D1B = r"D:\causal2drug\deliverables\d1b"
RAW = r"D:\causal2drug\data\raw"
OUT = r"D:\causal2drug\deliverables\d1e"
FLANK = 100_000

RESTRICTED = {"Tissue enriched", "Group enriched", "Tissue enhanced"}
NOT_RESTRICTED = {"Low tissue specificity"}


def target_states():
    """Recompute ABSENT / TESTED-NULL / INSTRUMENTED / UNMAPPED from raw artefacts."""
    frozen = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
    targets = sorted({r["gene"] for r in csv.DictReader(
        open(frozen, encoding="utf-8"), delimiter="\t")})
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D1B + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    tested = {r["Gene"] for r in csv.DictReader(
        open(D1B + r"\eqtlgen_tested_universe.tsv", encoding="utf-8"), delimiter="\t")}
    hits = set()
    with open(D1B + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = coords.get(r["MinikelSymbol"])
            if not c or r["SNPChr"] != c["chrom"]:
                continue
            if int(c["start"]) - FLANK <= int(r["SNPPos"]) <= int(c["end"]) + FLANK:
                hits.add(r["MinikelSymbol"])
    st = {}
    for s in targets:
        if s not in coords:
            st[s] = ("UNMAPPED", None)
        elif coords[s]["ensembl_id"] not in tested:
            st[s] = ("ABSENT", coords[s]["ensembl_id"])
        elif s in hits:
            st[s] = ("INSTRUMENTED", coords[s]["ensembl_id"])
        else:
            st[s] = ("TESTED_NULL", coords[s]["ensembl_id"])
    return st


def load_hpa():
    z = zipfile.ZipFile(RAW + r"\proteinatlas.tsv.zip")
    out = {}
    with z.open("proteinatlas.tsv") as fh:
        rd = csv.DictReader((l.decode("utf-8") for l in fh), delimiter="\t")
        for r in rd:
            e = r.get("Ensembl", "").strip()
            if e:
                out[e] = r
    return out


def load_gtex_blood():
    out = {}
    with gzip.open(RAW + r"\gtex_v10_median_tpm.gct.gz", "rt", encoding="utf-8") as fh:
        fh.readline(); fh.readline()
        hdr = fh.readline().rstrip("\n").split("\t")
        ib = hdr.index("Whole_Blood")
        tis = [i for i in range(2, len(hdr))]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = f[0].split(".")[0]
            try:
                vals = [float(f[i]) for i in tis]
                out[g] = (float(f[ib]), max(vals))
            except (ValueError, IndexError):
                pass
    return out


def fisher(a, b, c, d, label_rows, label_cols):
    """2x2: [[a,b],[c,d]]. Returns dict with OR, CI, p."""
    table = [[a, b], [c, d]]
    orr, p = stats.fisher_exact(table)
    # Woolf 95% CI on log OR with Haldane-Anscombe correction if any zero
    import math
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a + .5, b + .5, c + .5, d + .5)
    lor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return dict(table=table, rows=label_rows, cols=label_cols,
                odds_ratio=orr, ci95=[math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se)],
                p_fisher=p)


def show(name, res):
    (a, b), (c, d) = res["table"]
    print("\n%s" % name)
    print("            %-22s %-22s" % (res["cols"][0], res["cols"][1]))
    print("  %-9s %-22d %-22d" % (res["rows"][0], a, b))
    print("  %-9s %-22d %-22d" % (res["rows"][1], c, d))
    print("  OR = %.3f  95%% CI [%.3f, %.3f]   Fisher exact p = %.3g"
          % (res["odds_ratio"], res["ci95"][0], res["ci95"][1], res["p_fisher"]))


def main():
    os.makedirs(OUT, exist_ok=True)
    st = target_states()
    hpa = load_hpa()
    gtex = load_gtex_blood()
    print("Minikel targets: %d" % len(st))
    from collections import Counter
    print("states:", dict(Counter(v[0] for v in st.values())))
    print("HPA entries: %d ; GTEx genes: %d" % (len(hpa), len(gtex)))

    rows = []
    for sym, (state, ensg) in sorted(st.items()):
        h = hpa.get(ensg) if ensg else None
        g = gtex.get(ensg) if ensg else None
        spec = (h.get("RNA tissue specificity") or "").strip() if h else ""
        secl = (h.get("Secretome location") or "").strip() if h else ""
        pcls = (h.get("Protein class") or "") if h else ""
        rows.append(dict(
            symbol=sym, state=state, ensembl=ensg or "",
            in_hpa=bool(h), rna_tissue_specificity=spec,
            tissue_restricted=("YES" if spec in RESTRICTED else
                               "NO" if spec in NOT_RESTRICTED else "NA"),
            secretome_location=secl,
            secreted=("YES" if (secl or "secreted" in pcls.lower()) else "NO") if h else "NA",
            blood_tpm=(g[0] if g else ""), max_tissue_tpm=(g[1] if g else "")))
    with open(OUT + r"\classA_annotation.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    mapped = [r for r in rows if r["state"] != "UNMAPPED"]
    inhpa = [r for r in mapped if r["in_hpa"]]
    print("of %d mapped targets, %d (%.1f%%) join to HPA"
          % (len(mapped), len(inhpa), 100.0 * len(inhpa) / len(mapped)))

    def counts(states, field, val):
        sel = [r for r in inhpa if r["state"] in states and r[field] != "NA"]
        yes = sum(1 for r in sel if r[field] == val)
        return yes, len(sel) - yes

    results = {}
    # H1: tissue-restricted
    a, b = counts(["ABSENT"], "tissue_restricted", "YES")
    c, d = counts(["INSTRUMENTED"], "tissue_restricted", "YES")
    results["absent_vs_instrumented_tissue_restricted"] = fisher(
        a, b, c, d, ["ABSENT", "INSTR"], ["tissue-restricted", "low-specificity"])
    show("H1  ABSENT vs INSTRUMENTED  x  tissue-restricted",
         results["absent_vs_instrumented_tissue_restricted"])

    a, b = counts(["TESTED_NULL"], "tissue_restricted", "YES")
    results["testednull_vs_instrumented_tissue_restricted"] = fisher(
        a, b, c, d, ["TESTNUL", "INSTR"], ["tissue-restricted", "low-specificity"])
    show("H1b TESTED-NULL vs INSTRUMENTED  x  tissue-restricted",
         results["testednull_vs_instrumented_tissue_restricted"])

    # H2: secreted
    a, b = counts(["ABSENT"], "secreted", "YES")
    c, d = counts(["INSTRUMENTED"], "secreted", "YES")
    results["absent_vs_instrumented_secreted"] = fisher(
        a, b, c, d, ["ABSENT", "INSTR"], ["secreted", "not secreted"])
    show("H2  ABSENT vs INSTRUMENTED  x  secreted", results["absent_vs_instrumented_secreted"])

    a, b = counts(["TESTED_NULL"], "secreted", "YES")
    results["testednull_vs_instrumented_secreted"] = fisher(
        a, b, c, d, ["TESTNUL", "INSTR"], ["secreted", "not secreted"])
    show("H2b TESTED-NULL vs INSTRUMENTED  x  secreted",
         results["testednull_vs_instrumented_secreted"])

    # H3: whole-blood expression (GTEx), the mechanistically direct quantity
    def tpms(state):
        return [float(r["blood_tpm"]) for r in rows
                if r["state"] == state and r["blood_tpm"] != ""]
    ab, ins, tn = tpms("ABSENT"), tpms("INSTRUMENTED"), tpms("TESTED_NULL")
    import statistics as sstat
    u1, p1 = stats.mannwhitneyu(ab, ins, alternative="two-sided")
    u2, p2 = stats.mannwhitneyu(tn, ins, alternative="two-sided")
    print("\nH3  GTEx v10 Whole_Blood median TPM (independent cross-check)")
    for lab, v in (("ABSENT", ab), ("TESTED_NULL", tn), ("INSTRUMENTED", ins)):
        v2 = sorted(v)
        print("  %-13s n=%-5d median=%-9.3f  Q1=%-9.3f Q3=%-9.3f  %% with TPM<1: %.1f%%"
              % (lab, len(v), sstat.median(v), v2[len(v2) // 4], v2[3 * len(v2) // 4],
                 100.0 * sum(1 for x in v if x < 1) / len(v)))
    print("  Mann-Whitney ABSENT vs INSTRUMENTED      : U=%.0f p=%.3g" % (u1, p1))
    print("  Mann-Whitney TESTED_NULL vs INSTRUMENTED : U=%.0f p=%.3g" % (u2, p2))
    results["blood_tpm"] = dict(
        median_absent=sstat.median(ab), median_tested_null=sstat.median(tn),
        median_instrumented=sstat.median(ins),
        mw_absent_vs_instr_p=p1, mw_testednull_vs_instr_p=p2,
        n_absent=len(ab), n_tested_null=len(tn), n_instrumented=len(ins))

    json.dump(results, open(OUT + r"\classA_tests.json", "w"), indent=2, default=str)
    print("\nwrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
