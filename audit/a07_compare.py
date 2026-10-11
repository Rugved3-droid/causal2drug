"""AUDIT step 3: compare my outputs with the original, target by target.

Original files read here (first opened in Step 3):
  manuscript_prep/analysis/out/per_target_primary.csv     per-target eQTLGen status
  manuscript_prep/manuscript/recovery_map.csv             recovery map
  manuscript_prep/manuscript/rescue_controls/rescue_controls_results.csv
Run:  python audit/a07_compare.py
Writes audit/out/compare_status.csv, compare_recovery.csv
"""
import csv, collections, json, os

ROOT = r"D:\causal2drug"; OUT = os.path.join(ROOT, "audit", "out"); MP = os.path.join(ROOT, "manuscript_prep")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8")))
T = lambda v: str(v).strip().lower() == "true"

# ---- (a)(b) status ---------------------------------------------------------
mine = {r["symbol"]: r for r in rd(os.path.join(OUT, "audit_per_target.csv"))}
orig = {r["gene"]: r for r in rd(os.path.join(MP, "analysis", "out", "per_target_primary.csv"))}
print(f"[E01] target set: mine {len(mine)}, original {len(orig)}, only mine {sorted(set(mine)-set(orig))}, only original {sorted(set(orig)-set(mine))}")
SUB = {"autosomal gene not tested": "autosomal", "X-linked: chromosome not analysed by eQTLGen": "x_linked",
       "mitochondrial: genome not analysed by eQTLGen": "mitochondrial", "": ""}
PH = {"Preclinical": "1", "Phase I": "2", "Phase II": "3", "Phase III": "4", "Launched": "5"}
diff = collections.defaultdict(list)
with open(os.path.join(OUT, "compare_status.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["gene", "field", "mine", "original"])
    for g in sorted(set(mine) & set(orig)):
        m, o = mine[g], orig[g]
        pairs = [("status", m["status"], o["status"]), ("subtype", m["subtype"], SUB.get(o["no_eqtl_data_subtype"], o["no_eqtl_data_subtype"])),
                 ("ensembl_id", m["ensembl_id"], o["ensembl_id"]), ("chrom", m["chrom"], o["chromosome"]),
                 ("phase", m["phase"], PH.get(o["phase_PRIMARY"], o["phase_PRIMARY"])),
                 ("n_qualifying", str(m["n_qualifying"] or 0), str(o["n_cis_eqtl_p5e8_within_100kb"] or 0)),
                 ("lead_snp", m["lead_snp"], o["best_snp_within_100kb"])]
        for f, a, b in pairs:
            if a != b:
                diff[f].append(g); w.writerow([g, f, a, b])
for f in ("status", "subtype", "ensembl_id", "chrom", "phase", "n_qualifying", "lead_snp"):
    print(f"[E02] {f}: {len(set(mine)&set(orig))-len(diff[f])} match, {len(diff[f])} differ {diff[f][:12]}")

# ---- (c) recovery map ------------------------------------------------------
mr = {r["symbol"]: r for r in rd(os.path.join(OUT, "audit_recovery_map.csv"))}
orr = {r["gene"]: r for r in rd(os.path.join(MP, "manuscript", "recovery_map.csv"))}
print(f"[E03] recovery scope: mine {len(mr)}, original {len(orr)}, only mine {sorted(set(mr)-set(orr))}, only original {sorted(set(orr)-set(mr))}")
sset = lambda s: set(x.strip() for x in s.replace(";", "|").split("|") if x.strip())
rdiff = collections.defaultdict(list)
with open(os.path.join(OUT, "compare_recovery.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["gene", "field", "mine", "original"])
    for g in sorted(set(mr) & set(orr)):
        m, o = mr[g], orr[g]
        o_dec_inst = T(o["decode_instrument"])
        pairs = [("category", m["category"], o["category"]), ("tier", m["tier"], o["reporting_tier"]),
                 ("gtex_id_used", m["gtex_id_used"], o["gtex_ensembl_id"] if int(o["gtex_n_tissues_tested"] or 0) else ""),
                 ("gtex_n_tissues_tested", m["gtex_n_tissues_tested"], o["gtex_n_tissues_tested"]),
                 ("gtex_nonblood_tissues", sorted(sset(m["gtex_nonblood_tissues"])), sorted(sset(o["gtex_tissues_with_instrument_excl_blood"]))),
                 ("gtex_cv_tissues", sorted(sset(m["gtex_cv_tissues"])), sorted(sset(o["gtex_cv_relevant_tissues_with_instrument"]))),
                 ("gtex_whole_blood", T(m["gtex_whole_blood"]), T(o["gtex_whole_blood_instrument"])),
                 ("gtex_best_cv_tissue", m["gtex_best_cv_tissue"], o["gtex_best_cv_tissue"]),
                 ("gtex_best_cv_variant", m["gtex_best_cv_variant"], o["gtex_best_cv_variant_b38"]),
                 ("ukb_assayed", T(m["ukb_assayed"]), T(o["ukbppp_assayed"])), ("ukb_instrument", T(m["ukb_instrument"]), T(o["ukbppp_instrument"])),
                 ("ukb_sheet", m["ukb_sheet"], o["ukbppp_note"] if T(o["ukbppp_instrument"]) else ""),
                 ("decode_assayed", T(m["decode_listed"]), T(o["decode_assayed"])), ("decode_instrument", T(m["decode_instrument"]), o_dec_inst),
                 ("fenland_assayed", T(m["fenland_assayed"]), T(o["fenland_assayed"])), ("fenland_instrument", T(m["fenland_instrument"]), T(o["fenland_instrument"])),
                 ("x_flag", T(m["x_flag"]), T(o["x_linked"]))]
        import re
        opos = lambda v: (re.search(r":(\d+)\)", v).group(1) if re.search(r":(\d+)\)", v) else "")
        if T(m["ukb_instrument"]) and T(o["ukbppp_instrument"]):
            pairs.append(("ukb_lead_pos38", str(m["ukb_pos38"]), opos(o["ukbppp_lead_variant"])))
        if T(m["decode_instrument"]) and o_dec_inst:
            pairs.append(("decode_lead", m["decode_lead"], o["decode_lead_variant"].split(" (")[0]))
        if T(m["fenland_instrument"]) and T(o["fenland_instrument"]):
            pairs.append(("fenland_lead_pos37", re.search(r":(\d+)", m["fenland_lead"]).group(1), opos(o["fenland_lead_variant"])))
        for f, a, b in pairs:
            if a != b:
                rdiff[f].append(g); w.writerow([g, f, a, b])
n = len(set(mr) & set(orr))
for f in ["category", "tier", "gtex_id_used", "gtex_n_tissues_tested", "gtex_nonblood_tissues", "gtex_cv_tissues", "gtex_whole_blood",
          "gtex_best_cv_tissue", "gtex_best_cv_variant", "ukb_assayed", "ukb_instrument", "ukb_sheet", "decode_assayed",
          "decode_instrument", "fenland_assayed", "fenland_instrument", "x_flag", "ukb_lead_pos38", "decode_lead", "fenland_lead_pos37"]:
    print(f"[E04] {f}: {n-len(rdiff[f])} match, {len(rdiff[f])} differ {rdiff[f][:25]}")
C = collections.Counter
print("[E05] category mine", dict(sorted(C(r["category"] for r in mr.values()).items())), "| original", dict(sorted(C(r["category"] for r in orr.values()).items())))
print("[E05] tier mine", dict(sorted(C(r["tier"] for r in mr.values()).items())), "| original", dict(sorted(C(r["reporting_tier"] for r in orr.values()).items())))

# ---- (d) controls ----------------------------------------------------------
ctl = json.load(open(os.path.join(OUT, "audit_controls.json")))
oc = {r["control"]: r for r in rd(os.path.join(MP, "manuscript", "rescue_controls", "rescue_controls_results.csv"))}
o = oc["PCSK9_LDL"]; m = ctl["PCSK9"]
print(f"[E06] PCSK9->LDL-C mine  beta {m['mr_beta']:.4f} se {m['mr_se']:.4f} F {m['F']:.1f} exposure {m['exposure_sheet']} beta {m['exposure_beta']} se {m['exposure_se']}; outcome beta {m['outcome_beta_aligned']} se {m['outcome_se']}")
print(f"[E06] PCSK9->LDL-C orig  beta {float(o['mr_beta']):.4f} se {float(o['mr_se']):.4f} F {float(o['F']):.1f} exposure {o['resource']} beta {o['beta_exposure']} se {o['se_exposure']}; outcome beta {o['beta_outcome_harmonized']} se {o['se_outcome']}; study {o['outcome_study']}")
alt = m["outcome_beta_aligned"] / float(o["beta_exposure"])
print(f"[E06] reconciliation: my outcome beta / the original's ST9 exposure beta = {alt:.4f} (SE {abs(m['outcome_se']/float(o['beta_exposure'])):.4f}) -> the whole difference is the sheet")
m = ctl["HMGCR"]
print(f"[E07] HMGCR->LDL-C mine  beta {m['mr_beta']:.4f} se {m['mr_se']:.4f} p {m['mr_p']:.3g} F {m['F']:.1f} lead {m['exposure_variant']}  | KEY_NUMBERS: beta +0.389 (SE 0.012), p 4.2e-229, F 224")
