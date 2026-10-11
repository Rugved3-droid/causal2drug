#!/usr/bin/env python
"""Deliverable 1b: reconcile the clump output onto the PRIMARY (ENSG) mapping join.

The main clump ran on the symbol-join hit set. The primary mapping method is the
stable-identifier join (Minikel symbol -> GRCh37 ENSG -> eQTLGen `Gene`). The two differ
for a handful of targets:
  - targets the symbol join missed entirely (eQTLGen uses Ensembl v71 symbols)
  - targets whose symbol is ambiguous in eQTLGen (maps to >1 ENSG) but which resolve to
    exactly one ENSG in the GRCh37 GTF

Those targets are clumped here under the ENSG join and merged in, so the census rests on
one consistent mapping method throughout.
"""
import collections
import csv
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d1b_clump as C  # noqa: E402

D = r"D:\causal2drug\deliverables\d1b"


def load_ensg_hits(coords):
    hits = collections.defaultdict(list)
    with open(D + r"\eqtlgen_hits_by_ensg.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = coords.get(r["MinikelSymbol"])
            if not c or r["SNPChr"] != c["chrom"]:
                continue
            p = int(r["SNPPos"])
            if int(c["start"]) - C.FLANK <= p <= int(c["end"]) + C.FLANK:
                hits[r["MinikelSymbol"]].append(r)
    return hits


def main():
    keep = C.eur_samples()
    coords = {r["symbol"]: r for r in csv.DictReader(
        open(D + r"\gene_coords_grch37.tsv", encoding="utf-8"), delimiter="\t")}
    ensg_hits = load_ensg_hits(coords)
    _, sym_hits = C.load_inputs()

    existing = {r["symbol"]: r for r in csv.DictReader(
        open(D + r"\clumped_instruments.tsv", encoding="utf-8"), delimiter="\t")}

    todo = []
    for s in ensg_hits:
        if s not in existing:
            todo.append((s, "absent_from_symbol_join"))
        elif len(ensg_hits[s]) != len(sym_hits.get(s, [])):
            todo.append((s, "hit_count_differs_under_ensg_join"))
    stale = [s for s in existing if s not in ensg_hits]

    print("targets under ENSG join (window)      :", len(ensg_hits))
    print("already clumped (symbol join)         :", len(existing))
    print("to (re)clump under ENSG join          :", len(todo))
    for s, why in sorted(todo):
        print("   %-10s %s" % (s, why))
    print("stale rows to drop (not in ENSG join) :", len(stale), stale)

    by_chrom = collections.defaultdict(list)
    for s, _ in todo:
        by_chrom[coords[s]["chrom"]].append(s)
    redone = {}
    for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
        ctx = C.ChromCtx(chrom, keep)
        for s in by_chrom[chrom]:
            r = C.clump_one(ctx, s, coords[s], ensg_hits[s])
            redone[s] = r
            print("   chr%-3s %-10s hits=%-5d in_panel=%-5d indep=%-3d best_p=%.3g"
                  % (chrom, s, r["n_hits"], r["n_in_ld_panel"], r["n_independent"],
                     r["best_p"]), flush=True)

    final = {}
    for s, r in existing.items():
        if s in ensg_hits and s not in redone:
            final[s] = r
    for s, r in redone.items():
        final[s] = {k: ("" if v is None else v) for k, v in r.items()}

    cols = ["symbol", "n_hits", "n_in_ld_panel", "n_independent", "best_p", "best_snp",
            "independent_snps", "note"]
    with open(D + r"\clumped_instruments_final.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for s in sorted(final):
            w.writerow({c: final[s].get(c, "") for c in cols})
    print("final rows:", len(final))
    assert set(final) == set(ensg_hits), "final clump set != ENSG-join hit set"
    print("final set == ENSG-join hit set: True")
    print("wrote", D + r"\clumped_instruments_final.tsv")


if __name__ == "__main__":
    main()
