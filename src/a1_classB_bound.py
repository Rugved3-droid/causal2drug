#!/usr/bin/env python
"""ANALYSIS 1 -- bound CLASS B rather than assigning it.

CLASS B = target is well instrumented in blood, but the drug-relevant mechanism operates
at a level total mRNA cannot represent (splicing, shedding, protein activity).
Only IL6R can be assigned directly. This bounds the class instead.

SOURCE AND WHY: GWAS Catalog associations (ontology-annotated full release). It carries,
per association, the mapped gene and the variant's VEP functional consequence in the
`CONTEXT` column -- already computed and versioned by the Catalog, so no VEP calls and no
consequence-calling of my own. ClinVar was rejected because it is dominated by rare
pathogenic variants rather than the common association signals that drug-target genetics
rests on; Open Targets genetics was rejected because its L2G scores are model outputs,
which would put a trained model inside the denominator of a bound.

RULE (mechanical): an instrumented CV target is flagged PROTEIN-ALTERING-SIGNAL if the
GWAS Catalog reports >=1 genome-wide association whose MAPPED_GENE includes that gene as
an exact token AND whose CONTEXT is a protein-altering consequence.

This is an UPPER BOUND on Class B, explicitly NOT an assignment: a coding association in
the gene shows the mechanism COULD act at protein level, not that the drug's does.
"""
import collections
import csv
import json
import re
import zipfile

RAW = r"D:\causal2drug\data\raw"
CV = r"D:\causal2drug\deliverables\cv"
OUT = r"D:\causal2drug\deliverables\cv"

CODING = {"missense_variant", "stop_gained", "frameshift_variant", "start_lost",
          "stop_lost", "inframe_deletion", "inframe_insertion",
          "protein_altering_variant", "splice_donor_variant",
          "splice_acceptor_variant", "coding_sequence_variant"}
CVTRAIT = ("cardio", "coronary", "heart", "myocard", "atrial", "arrhythm", "cholesterol",
           "lipid", "lipoprotein", "triglycerid", "blood pressure", "hypertens",
           "stroke", "atheroscler", "aneurysm", "thrombo", "fibrillation",
           "cardiac", "vascular", "artery", "aortic", "ischemi", "ischaemi")


def main():
    atlas = list(csv.DictReader(open(CV + r"\cardiovascular_atlas.csv", encoding="utf-8")))
    instr = {r["gene"]: r for r in atlas if r["state"] == "INSTRUMENTED"}
    print("instrumented CV targets: %d" % len(instr))

    z = zipfile.ZipFile(RAW + r"\gwas-catalog-associations-full.zip")
    name = z.namelist()[0]
    with z.open(name) as fh:
        hdr = fh.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
    i_g, i_c, i_s = (hdr.index("MAPPED_GENE"), hdr.index("CONTEXT"), hdr.index("SNPS"))
    i_t = hdr.index("MAPPED_TRAIT")
    i_p = hdr.index("P-VALUE")

    hits = collections.defaultdict(list)
    n_rows = 0
    with z.open(name) as fh:
        fh.readline()
        for line in fh:
            n_rows += 1
            f = line.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(f) <= max(i_g, i_c, i_t, i_p):
                continue
            ctx = f[i_c]
            if not ctx or not (set(re.split(r"[ x,]+", ctx)) & CODING):
                continue
            genes = set(re.split(r"[,;] *| - ", f[i_g]))
            for g in genes:
                g = g.strip()
                if g in instr:
                    try:
                        pv = float(f[i_p])
                    except ValueError:
                        pv = None
                    hits[g].append((f[i_s], ctx, f[i_t], pv))
    print("association rows scanned: %d" % n_rows)

    flagged = sorted(hits)
    print("\n=== UPPER BOUND ON CLASS B ===")
    print("  instrumented CV targets with >=1 protein-altering GWAS association:")
    print("    %d / %d = %.1f%%" % (len(flagged), len(instr),
                                    100 * len(flagged) / len(instr)))

    def phase_subset(pred):
        sub = {g: r for g, r in instr.items() if pred(r)}
        fl = [g for g in sub if g in hits]
        return len(fl), len(sub)

    a, b = phase_subset(lambda r: r["max_phase_num"] == "5")
    print("  ...restricted to LAUNCHED-drug targets   : %d / %d = %.1f%%"
          % (a, b, 100 * a / b if b else 0))
    c, d = phase_subset(lambda r: r["max_phase_num"] in ("4", "5"))
    print("  ...restricted to Phase III or Launched   : %d / %d = %.1f%%"
          % (c, d, 100 * c / d if d else 0))

    # CV-trait-restricted version
    cvhits = {g: [h for h in v if any(k in h[2].lower() for k in CVTRAIT)]
              for g, v in hits.items()}
    cvhits = {g: v for g, v in cvhits.items() if v}
    print("\n  tightened to CV-RELEVANT traits only:")
    print("    %d / %d = %.1f%% of instrumented CV targets"
          % (len(cvhits), len(instr), 100 * len(cvhits) / len(instr)))
    lc = [g for g in cvhits if instr[g]["max_phase_num"] == "5"]
    print("    ...of launched-drug targets: %d / %d" % (len(lc), b))

    print("\n=== WORKED EXAMPLES ===")
    for g in ("IL6R", "TYK2"):
        if g in hits:
            best = sorted(hits[g], key=lambda h: (h[3] if h[3] is not None else 1))[:3]
            print("  %-6s %d protein-altering associations; strongest:" % (g, len(hits[g])))
            for s, ctx, t, pv in best:
                print("        %-14s %-20s p=%-10s %s"
                      % (s, ctx[:20], ("%.1g" % pv) if pv else "NA", t[:44]))
        else:
            print("  %-6s no protein-altering association found" % g)

    print("\n=== FLAGGED LAUNCHED / PHASE III TARGETS (the bound's core) ===")
    print("  %-10s %-11s %-14s %-22s %s" % ("gene", "phase", "variant", "consequence", "trait"))
    rows_out = []
    for g in sorted(flagged, key=lambda x: (instr[x]["max_phase_num"] != "5", x)):
        r = instr[g]
        best = sorted(hits[g], key=lambda h: (h[3] if h[3] is not None else 1))[0]
        rows_out.append(dict(gene=g, max_clinical_phase=r["max_clinical_phase"],
                             n_protein_altering_assoc=len(hits[g]),
                             top_variant=best[0], consequence=best[1],
                             trait=best[2], p_value=best[3],
                             cv_relevant_trait=g in cvhits))
        if r["max_phase_num"] in ("4", "5"):
            print("  %-10s %-11s %-14s %-22s %s"
                  % (g, r["max_clinical_phase"], best[0][:14], best[1][:22], best[2][:40]))

    with open(OUT + r"\classB_upper_bound.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows_out)
    json.dump(dict(n_instrumented=len(instr), n_flagged=len(flagged),
                   pct=100 * len(flagged) / len(instr),
                   launched_flagged=a, launched_total=b,
                   ph3plus_flagged=c, ph3plus_total=d,
                   cv_trait_flagged=len(cvhits)),
              open(OUT + r"\classB_bound.json", "w"), indent=2)
    print("\nwrote classB_upper_bound.csv (%d rows)" % len(rows_out))


if __name__ == "__main__":
    main()
