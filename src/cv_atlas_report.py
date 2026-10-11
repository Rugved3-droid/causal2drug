#!/usr/bin/env python
"""CV atlas: fix the parquet write, and compute items 4, 5 and 6."""
import collections
import csv
import json
import os

OUT = r"D:\causal2drug\deliverables\cv"
D1F = r"D:\causal2drug\deliverables\d1f"
D1J = r"D:\causal2drug\deliverables\d1j"
FROZEN = r"D:\causal2drug\data\derived\pp_dedup_frozen.tsv"
RAW = r"D:\causal2drug\data\raw"

NUM = {"n_cv_pairs": "int", "max_phase_num": "int",
       "n_independent_instruments": "int", "best_eqtl_p": "float",
       "gtex_blood_tpm": "float", "gtex_heart_LV_tpm": "float",
       "gtex_liver_tpm": "float", "gtex_artery_aorta_tpm": "float"}
BOOL = ["cv_primary", "cv_extended", "instrumented_blood", "instrumented_heart_LV",
        "instrumented_heart_atrial", "instrumented_artery_aorta",
        "instrumented_liver", "instrumented_cardiac_union"]


def main():
    rows = list(csv.DictReader(open(OUT + r"\cardiovascular_atlas.csv", encoding="utf-8")))
    for r in rows:
        for k, t in NUM.items():
            v = r.get(k, "")
            if v == "" or v is None:
                r[k] = None
            else:
                r[k] = int(float(v)) if t == "int" else float(v)
        for k in BOOL:
            r[k] = str(r.get(k)).lower() == "true"
    import pyarrow as pa, pyarrow.parquet as pq
    pq.write_table(pa.Table.from_pylist(rows), OUT + r"\cardiovascular_atlas.parquet")
    print("parquet written: %d rows, %d cols" % (len(rows), len(rows[0])))

    # ---------------- item 4: named examples ----------------
    print("\n=== 4. CLASS A / CLASS B CV TARGETS, by clinical phase reached ===")
    ca = [r for r in rows if r["untestability_class"] == "CLASS_A"]
    cb = [r for r in rows if r["untestability_class"] == "CLASS_B"]
    print("CLASS A total %d ; CLASS B total %d" % (len(ca), len(cb)))
    launched = [r for r in ca if r["max_phase_num"] == 5]
    ph3 = [r for r in ca if r["max_phase_num"] == 4]
    print("  CLASS A reaching Launched : %d" % len(launched))
    print("  CLASS A reaching Phase III: %d" % len(ph3))
    hdr = "  %-10s %-13s %-11s %-9s %-8s %-9s %s"
    print("\nCLASS A, LAUNCHED drugs (the atlas core):")
    print(hdr % ("gene", "state", "phase", "bloodTPM", "heartTPM", "liverTPM", "required modality"))
    for r in sorted(launched, key=lambda x: -(x["gtex_heart_LV_tpm"] or 0)):
        print(hdr % (r["gene"], r["state"], r["max_clinical_phase"],
                     ("%.2f" % r["gtex_blood_tpm"]) if r["gtex_blood_tpm"] is not None else "-",
                     ("%.1f" % r["gtex_heart_LV_tpm"]) if r["gtex_heart_LV_tpm"] is not None else "-",
                     ("%.1f" % r["gtex_liver_tpm"]) if r["gtex_liver_tpm"] is not None else "-",
                     r["required_modality"][:44]))
    print("\nCLASS A, PHASE III:")
    for r in sorted(ph3, key=lambda x: x["gene"])[:25]:
        print(hdr % (r["gene"], r["state"], r["max_clinical_phase"],
                     ("%.2f" % r["gtex_blood_tpm"]) if r["gtex_blood_tpm"] is not None else "-",
                     ("%.1f" % r["gtex_heart_LV_tpm"]) if r["gtex_heart_LV_tpm"] is not None else "-",
                     ("%.1f" % r["gtex_liver_tpm"]) if r["gtex_liver_tpm"] is not None else "-",
                     r["required_modality"][:44]))

    # ---------------- item 5: modality summary ----------------
    print("\n=== 5. REQUIRED MODALITY for untestable CV targets (N=%d) ===" % len(ca))
    m = collections.Counter(r["required_modality"].split(" (")[0] for r in ca)
    for k, v in m.most_common():
        print("  %-46s %4d  (%.1f%%)" % (k, v, 100 * v / len(ca)))

    # ---------------- item 6: method availability, CV-restricted ----------------
    print("\n=== 6. METHOD AVAILABILITY, CV-restricted ===")
    areas = {r["topl"]: r["area"] for r in csv.DictReader(
        open(RAW + r"\minikel_areas.tsv", encoding="utf-8"), delimiter="\t")}
    topl = {r["indication_mesh_id"]: r["topl"] for r in csv.DictReader(
        open(RAW + r"\minikel_indic_topl_match.tsv", encoding="utf-8"), delimiter="\t")}
    TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
           "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")

    def is_cv(mesh_id, term):
        return (areas.get(topl.get(mesh_id)) == "cardiovascular"
                or any(k in term.lower() for k in TOK))

    pairs = list(csv.DictReader(open(D1F + r"\partc_pairs.tsv", encoding="utf-8"),
                               delimiter="\t"))
    cvpairs = [p for p in pairs if is_cv(p["mesh_id"], p["mesh_term"])]
    res = list(csv.DictReader(open(D1F + r"\partc_results.tsv", encoding="utf-8"),
                             delimiter="\t"))
    got = {(r["gene"], r["mesh_id"], r["gwas_acc"]) for r in res}
    lost = [p for p in cvpairs if (p["gene"], p["mesh_id"], p["gwas_acc"]) not in got]
    print("  CV pairs qualifying for the Part C mechanism analysis : %d" % len(cvpairs))
    print("  ...LOST to unindexed plain-gzip GWAS accessions        : %d (%.1f%%)"
          % (len(lost), 100 * len(lost) / len(cvpairs) if cvpairs else 0))

    r15 = list(csv.DictReader(open(D1J + r"\r15_results.tsv", encoding="utf-8"),
                              delimiter="\t"))
    cvr15 = [r for r in r15 if r["gene"] in {x["gene"] for x in rows}]
    sig = [r for r in cvr15 if r["gwas_genome_wide_sig"] == "True"]
    usable = [r for r in sig if r["susie_PP4"]]
    print("  CV target loci screened for colocalization             : %d" % len(cvr15))
    print("  ...with a genome-wide significant outcome signal       : %d (%.1f%%)"
          % (len(sig), 100 * len(sig) / len(cvr15) if cvr15 else 0))
    print("  ...yielding a usable credible-set colocalization       : %d (%.1f%%)"
          % (len(usable), 100 * len(usable) / len(cvr15) if cvr15 else 0))
    print("  => colocalization UNAVAILABLE for %.1f%% of CV target loci"
          % (100 * (1 - len(usable) / len(cvr15)) if cvr15 else 0))

    json.dump(dict(class_a=len(ca), class_b=len(cb),
                   class_a_launched=len(launched), class_a_phase3=len(ph3),
                   modality=dict(m),
                   cv_pairs_partc=len(cvpairs), cv_pairs_lost_unindexed=len(lost),
                   cv_loci_screened=len(cvr15), cv_loci_gwas_sig=len(sig),
                   cv_loci_coloc_usable=len(usable)),
              open(OUT + r"\cv_report.json", "w"), indent=2)


if __name__ == "__main__":
    main()
