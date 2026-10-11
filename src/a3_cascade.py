#!/usr/bin/env python
"""ANALYSIS 3 -- verify the denominator behind '57/103 CV pairs lost to unindexed
accessions'.

Traces every drop from 630 CV targets to the pair set that claim was computed on, and
then asks whether 103 is actually the right denominator for an infrastructure-attrition
rate, or whether a different denominator is correct.
"""
import collections
import csv
import re

RAW = r"D:\causal2drug\data\raw"
D1B = r"D:\causal2drug\deliverables\d1b"
D1F = r"D:\causal2drug\deliverables\d1f"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")


def tsv(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def main():
    areas = {r["topl"]: r["area"] for r in tsv(RAW + r"\minikel_areas.tsv")}
    topl = {r["indication_mesh_id"]: r["topl"]
            for r in tsv(RAW + r"\minikel_indic_topl_match.tsv")}

    def is_cv(mid, term):
        return (areas.get(topl.get(mid)) == "cardiovascular"
                or any(k in term.lower() for k in TOK))

    frozen = tsv(FROZEN)
    cv_pairs = [r for r in frozen if is_cv(r["indication_mesh_id"],
                                           r["indication_mesh_term"])]
    cv_targets = {r["gene"] for r in cv_pairs}
    print("=" * 78)
    print("CASCADE: CV pairs -> the denominator used for infrastructure attrition")
    print("=" * 78)
    print("  %-62s %6d" % ("CV-EXTENDED target-indication pairs (frozen set)", len(cv_pairs)))
    print("  %-62s %6d" % ("...distinct CV targets", len(cv_targets)))

    # step: instrumented in eQTLGen
    clump = {r["symbol"] for r in tsv(D1B + r"\clumped_instruments_final.tsv")}
    s1 = [r for r in cv_pairs if r["gene"] in clump]
    print("  %-62s %6d" % ("pairs whose target is INSTRUMENTED in eQTLGen", len(s1)))
    print("  %-62s %6d" % ("...distinct targets", len({r['gene'] for r in s1})))

    # step: target has an unambiguous approved-drug direction (Part C R13.2)
    import json
    direction = json.load(open(D1F + r"\target_direction.json"))
    coords = {r["symbol"]: r["ensembl_id"] for r in tsv(D1B + r"\gene_coords_grch37.tsv")}
    s2 = [r for r in s1 if coords.get(r["gene"]) in direction]
    print("  %-62s %6d" % ("...AND target has a consistent approved-drug direction", len(s2)))
    print("  %-62s %6d" % ("...distinct targets", len({r['gene'] for r in s2})))

    # step: indication exact-matches a EUR-only harmonised GWAS trait (R13.4)
    pairs_pc = tsv(D1F + r"\partc_pairs.tsv")
    pc_keys = {(r["gene"], r["mesh_id"]) for r in pairs_pc}
    s3 = [r for r in s2 if (r["gene"], r["indication_mesh_id"]) in pc_keys]
    print("  %-62s %6d" % ("...AND indication exact-matches a EUR harmonised GWAS trait",
                           len(s3)))
    print("  %-62s %6d" % ("...distinct targets", len({r['gene'] for r in s3})))

    cv_pc = [p for p in pairs_pc if is_cv(p["mesh_id"], p["mesh_term"])]
    print("\n  cross-check against partc_pairs.tsv directly: %d CV pairs" % len(cv_pc))

    # who was lost, and why
    res = tsv(D1F + r"\partc_results.tsv")
    got = {(r["gene"], r["mesh_id"], r["gwas_acc"]) for r in res}
    lost = [p for p in cv_pc if (p["gene"], p["mesh_id"], p["gwas_acc"]) not in got]
    skipped = {r["gwas_acc"]: r["reason"] for r in tsv(D1F + r"\partc_skipped_accessions.tsv")}
    reasons = collections.Counter(skipped.get(p["gwas_acc"], "NOT_IN_SKIPPED_LIST")
                                 for p in lost)
    print("\n" + "=" * 78)
    print("THE 103 DENOMINATOR, AUDITED")
    print("=" * 78)
    print("  CV pairs entering the Part C run          : %d" % len(cv_pc))
    print("  CV pairs with a result row written        : %d" % (len(cv_pc) - len(lost)))
    print("  CV pairs LOST                             : %d" % len(lost))
    for k, v in reasons.most_common():
        print("      reason: %-42s %d" % (k, v))

    # Is 103 the right denominator? It is conditioned on Part C's mechanism + trait rules.
    # The infrastructure question is not conditioned on those. Recompute on the broader set.
    print("\n" + "=" * 78)
    print("IS 103 THE RIGHT DENOMINATOR?")
    print("=" * 78)
    print("  103 is conditioned on TWO Part C requirements that have nothing to do with")
    print("  file indexing: (a) the target must carry an unambiguous approved-drug")
    print("  mechanism direction, and (b) the indication must exact-match a GWAS trait.")
    print("  Those are Part C eligibility rules, not infrastructure. Quoting an")
    print("  infrastructure-attrition rate over them mixes two different filters.")

    # Correct denominator: CV pairs whose indication maps to a EUR-only harmonised GWAS
    # at all -- i.e. every CV pair for which an outcome GWAS could in principle be used.
    harm = set()
    for l in open(RAW + r"\harmonised_list.txt"):
        m = re.search(r"/(GCST\d+)/harmonised/", l)
        if m:
            harm.add(m.group(1))
    ANC = ["African", "East Asian", "South Asian", "Hispanic", "Latin", "Native",
           "Asian unspecified", "Oceanian", "Greater Middle Eastern", "Aboriginal",
           "Sub-Saharan", "NR"]

    def eur_only(s):
        return ("European" in s) and not any(a in s for a in ANC)

    best = {}
    with open(RAW + r"\gwas-catalog-download-studies-v1.0.3.1.txt",
              encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        h = next(rd)
        it, ia, isz = (h.index("MAPPED_TRAIT"), h.index("STUDY ACCESSION"),
                       h.index("INITIAL SAMPLE SIZE"))
        for f in rd:
            if len(f) <= it or f[ia] not in harm or not eur_only(f[isz]):
                continue
            n = max([int(x.replace(",", "")) for x in re.findall(r"[\d,]{4,}", f[isz])] or [0])
            for t in f[it].split(","):
                t = t.strip().lower()
                if t and (t not in best or n > best[t][1]):
                    best[t] = (f[ia], n)
    cv_mappable = [r for r in cv_pairs
                   if r["indication_mesh_term"].strip().lower() in best]
    print("\n  BROADER denominator -- every CV pair whose indication maps to a EUR-only")
    print("  harmonised GWAS (no mechanism-direction requirement):")
    print("    CV pairs mappable to an outcome GWAS        : %d" % len(cv_mappable))
    accs = {r["gene"] + "|" + r["indication_mesh_id"]:
            best[r["indication_mesh_term"].strip().lower()][0] for r in cv_mappable}
    lost_broad = [k for k, a in accs.items() if a in skipped]
    print("    ...whose GWAS accession is unindexed        : %d" % len(lost_broad))
    print("    infrastructure attrition rate               : %.1f%%"
          % (100 * len(lost_broad) / len(cv_mappable) if cv_mappable else 0))
    print("\n  Narrow (Part C-conditioned) rate: %d/%d = %.1f%%"
          % (len(lost), len(cv_pc), 100 * len(lost) / len(cv_pc)))


if __name__ == "__main__":
    main()
