#!/usr/bin/env python
"""Deliverable 1b steps 4-7: instrument-yield census, pair propagation,
therapeutic-area breakdown, and control lookup.

Everything is recomputed from files on disk; no count is carried from a previous step.
"""
import collections
import csv
import json
import statistics

D = r"D:\causal2drug\deliverables\d1b"
R = r"D:\causal2drug\data\raw"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
FLANK = 100_000

PHASES = [("any phase", lambda c: True),
          ("Phase II+", lambda c: c >= 3),
          ("Phase III+", lambda c: c >= 4),
          ("launched", lambda c: c >= 5)]


def tsv(path, delim="\t"):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def main():
    out = {}
    frozen = tsv(FROZEN)
    targets = sorted({r["gene"] for r in frozen})
    coords = {r["symbol"]: r for r in tsv(D + r"\gene_coords_grch37.tsv")}
    tested = tsv(D + r"\eqtlgen_tested_universe.tsv")
    tested_ensg = {r["Gene"] for r in tested}
    clump = {r["symbol"]: r for r in tsv(D + r"\clumped_instruments_final.tsv")}

    # ---------- step 2/3: mapping ----------
    unmapped = [s for s in targets if s not in coords]
    mapped = [s for s in targets if s in coords]
    absent = [s for s in mapped if coords[s]["ensembl_id"] not in tested_ensg]
    intested = [s for s in mapped if coords[s]["ensembl_id"] in tested_ensg]

    hits = collections.defaultdict(list)
    for r in tsv(D + r"\eqtlgen_hits_by_ensg.tsv"):
        c = coords.get(r["MinikelSymbol"])
        if not c or r["SNPChr"] != c["chrom"]:
            continue
        p = int(r["SNPPos"])
        if int(c["start"]) - FLANK <= p <= int(c["end"]) + FLANK:
            hits[r["MinikelSymbol"]].append(r)

    with_hit = sorted(hits)
    no_hit = [s for s in intested if s not in hits]
    out["mapping"] = {
        "minikel_unique_targets": len(targets),
        "UNMAPPED_no_unique_GRCh37_gene_record": len(unmapped),
        "mapped_to_GRCh37_ENSG": len(mapped),
        "ABSENT_from_eQTLGen_tested_universe": len(absent),
        "present_in_eQTLGen_tested_universe": len(intested),
        "eQTLGen_tested_genes_total": len(tested_ensg),
    }

    # ---------- step 4: census ----------
    indep = {s: int(clump[s]["n_independent"]) for s in with_hit if s in clump}
    missing_clump = [s for s in with_hit if s not in clump]
    ge1 = sorted([s for s, n in indep.items() if n >= 1])
    counts = sorted(indep[s] for s in ge1)
    out["census"] = {
        "targets_with_ge1_cis_eQTL_p_lt_5e-8_in_window": len(with_hit),
        "targets_missing_from_clump_output": len(missing_clump),
        "targets_with_ge1_INDEPENDENT_instrument": len(ge1),
        "targets_with_NO_significant_cis_eQTL_though_tested": len(no_hit),
        "targets_ABSENT_from_eQTLGen": len(absent),
        "targets_UNMAPPED": len(unmapped),
        "instruments_per_target": {
            "median": statistics.median(counts) if counts else None,
            "q1": counts[len(counts) // 4] if counts else None,
            "q3": counts[3 * len(counts) // 4] if counts else None,
            "min": counts[0] if counts else None,
            "max": counts[-1] if counts else None,
            "n_with_ge3": sum(1 for c in counts if c >= 3),
            "n_with_ge2": sum(1 for c in counts if c >= 2),
            "n_eq1": sum(1 for c in counts if c == 1),
        },
    }
    assert (len(unmapped) + len(absent) + len(no_hit) + len(with_hit)
            == len(targets)), "census partition does not sum to N targets"
    out["census"]["partition_sums_to_N"] = True

    # ---------- step 5: propagate to pairs ----------
    areas = {r["topl"]: r["area"] for r in tsv(R + r"\minikel_areas.tsv")}
    topl = {r["indication_mesh_id"]: r["topl"]
            for r in tsv(R + r"\minikel_indic_topl_match.tsv")}
    for r in frozen:
        r["area"] = areas.get(topl.get(r["indication_mesh_id"]), "UNMAPPED")
    ok = set(ge1)
    surv = [r for r in frozen if r["gene"] in ok]
    out["pairs"] = {"total_pairs_frozen": len(frozen),
                    "pairs_whose_target_has_ge1_independent_instrument": len(surv)}

    def phase_table(rows):
        t = {}
        for name, fn in PHASES:
            t[name] = len({r["gene"] for r in rows if r["ccatnum"].isdigit()
                           and fn(int(r["ccatnum"]))})
        return t

    out["distinct_targets_by_phase"] = {
        "all_frozen_pairs": phase_table(frozen),
        "surviving_pairs": phase_table(surv)}

    by_area_all = collections.defaultdict(list)
    by_area_surv = collections.defaultdict(list)
    for r in frozen:
        by_area_all[r["area"]].append(r)
    for r in surv:
        by_area_surv[r["area"]].append(r)
    ta = {}
    for a in sorted(by_area_all, key=lambda x: -len(by_area_all[x])):
        allt = {r["gene"] for r in by_area_all[a]}
        survt = {r["gene"] for r in by_area_surv.get(a, [])}
        ta[a] = {"pairs_all": len(by_area_all[a]),
                 "pairs_surviving": len(by_area_surv.get(a, [])),
                 "distinct_targets_all": len(allt),
                 "distinct_targets_surviving": len(survt),
                 "target_yield_rate": round(len(survt) / len(allt), 4) if allt else None,
                 "by_phase_surviving": phase_table(by_area_surv.get(a, [])),
                 "by_phase_all": phase_table(by_area_all[a])}
    out["therapeutic_area"] = ta

    # ---------- step 7: controls ----------
    ctrl = {}
    for s in ("HMGCR", "PCSK9", "IL6R"):
        c = clump.get(s)
        ctrl[s] = {
            "in_GRCh37_GTF": s in coords,
            "ensembl_id": coords.get(s, {}).get("ensembl_id"),
            "in_eQTLGen_tested": (coords.get(s, {}).get("ensembl_id") in tested_ensg
                                  if s in coords else False),
            "n_hits_p_lt_5e-8_in_window": len(hits.get(s, [])),
            "n_independent_instruments": int(c["n_independent"]) if c else None,
            "best_p": float(c["best_p"]) if c else None,
            "best_snp": c["best_snp"] if c else None,
            "n_in_ld_panel": int(c["n_in_ld_panel"]) if c else None,
        }
    out["controls"] = ctrl

    json.dump(out, open(D + r"\census.json", "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
