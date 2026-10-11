#!/usr/bin/env python
"""READ-ONLY verification of every key number for manuscript_prep/NUMBERS.md.

    python manuscript_prep/verify/verify_numbers.py > manuscript_prep/verify/verify_output.txt

Opens existing project files for reading only. Writes nothing outside
manuscript_prep/verify/. Does not import or execute any script in src/ (several of
those overwrite their own outputs when run).

Each block is tagged [Vnn] so NUMBERS.md can cite it.
"""
import collections
import csv
import hashlib
import json
import math
import os
import re
import sys

from scipy import stats

ROOT = r"D:\causal2drug"
RAW = ROOT + r"\data\raw"
DEL = ROOT + r"\deliverables"
HERE = os.path.dirname(os.path.abspath(__file__))
FROZEN = ROOT + r"\data\derived\pp_dedup_frozen.tsv"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
DROP = {"DAO-D022124", "HRH1-D000086382", "IFNAR2-D012141", "MC1R-D010787"}
PHASE = {0: "(none)", 1: "Preclinical", 2: "Phase I", 3: "Phase II", 4: "Phase III",
         5: "Launched"}

sys.stdout.reconfigure(encoding="utf-8")


def tsv(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def hdr(tag, title):
    print("\n" + "=" * 78)
    print("[%s] %s" % (tag, title))
    print("=" * 78)


def pct(a, b):
    return "%.1f%%" % (100.0 * a / b) if b else "n/a"


def or_ci(a, b, c, d):
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a + .5, b + .5, c + .5, d + .5)
    lor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return orr, math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se), p


# ---------------------------------------------------------------- V01 flow
hdr("V01", "FLOW: raw pp.tsv -> dedup by ti_uid -> drop 4 conflicting programs")
lines = open(RAW + r"\pp.tsv", encoding="utf-8").read().rstrip("\n").split("\n")
h, data = lines[0].split("\t"), lines[1:]
iu = h.index("ti_uid")
by = collections.defaultdict(list)
for l in data:
    by[l.split("\t")[iu]].append(l)
kept = {k: v for k, v in by.items() if k not in DROP}
print("raw data rows                         : %d" % len(data))
print("distinct ti_uid                       : %d" % len(by))
print("after dropping 4 conflicting programs : %d" % len(kept))
print("rows removed in total                 : %d" % (len(data) - len(kept)))
for k in sorted(DROP):
    vals = sorted({l.split("\t")[h.index("ccat")] for l in by.get(k, [])})
    print("   dropped %-18s present=%s rows=%d ccat values=%s"
          % (k, k in by, len(by.get(k, [])), vals))
frozen = tsv(FROZEN)
print("pp_dedup_frozen.tsv data rows         : %d" % len(frozen))
print("pp_dedup_frozen.tsv sha256            : %s" % sha(FROZEN))
print("frozen ti_uid set == recomputed set   : %s"
      % ({r["ti_uid"] for r in frozen} == set(kept)))
print("distinct genes / MeSH indications     : %d / %d"
      % (len({r["gene"] for r in frozen}), len({r["indication_mesh_id"] for r in frozen})))
print("pp.tsv sha1                           : %s"
      % hashlib.sha1(open(RAW + r"\pp.tsv", "rb").read()).hexdigest())

# ---------------------------------------------------------------- V02 CV sets
hdr("V02", "CV SET: CV-PRIMARY (MeSH C14) and CV-EXTENDED (+ 8 lipid/athero tokens)")
areas = {r["topl"]: r["area"] for r in tsv(RAW + r"\minikel_areas.tsv")}
topl = {r["indication_mesh_id"]: r["topl"] for r in tsv(RAW + r"\minikel_indic_topl_match.tsv")}
miss = {r["indication_mesh_id"] for r in frozen if r["indication_mesh_id"] not in topl}
print("frozen indications with no top-level mapping: %d" % len(miss))
for r in frozen:
    r["cvp"] = areas.get(topl.get(r["indication_mesh_id"])) == "cardiovascular"
    r["cve"] = r["cvp"] or any(k in r["indication_mesh_term"].lower() for k in TOK)
cvp = [r for r in frozen if r["cvp"]]
cve = [r for r in frozen if r["cve"]]
tp, te = {r["gene"] for r in cvp}, {r["gene"] for r in cve}
print("CV-PRIMARY  pairs %d  targets %d  indications %d"
      % (len(cvp), len(tp), len({r["indication_mesh_id"] for r in cvp})))
print("CV-EXTENDED pairs %d  targets %d  indications %d  (+%d targets)"
      % (len(cve), len(te), len({r["indication_mesh_id"] for r in cve}), len(te - tp)))
ext_only = sorted({r["indication_mesh_term"] for r in cve if not r["cvp"]})
print("indications added by the token rule (%d): %s" % (len(ext_only), "; ".join(ext_only)))
print("PCSK9 in CV-PRIMARY: %s ; in CV-EXTENDED: %s" % ("PCSK9" in tp, "PCSK9" in te))

# ---------------------------------------------------------------- V03 states
hdr("V03", "STATE / CLASS per CV target, recomputed from d1e/classA_annotation.tsv")
ann = {r["symbol"]: r for r in tsv(DEL + r"\d1e\classA_annotation.tsv")}
gw = collections.Counter(r["state"] for r in ann.values())
print("genome-wide (all %d Minikel targets): %s" % (len(ann), dict(gw)))
print("   + targets in frozen set not in annotation: %d"
      % len({r["gene"] for r in frozen} - set(ann)))
maxph = collections.defaultdict(int)
for r in cve:
    if r["ccatnum"].isdigit():
        maxph[r["gene"]] = max(maxph[r["gene"]], int(r["ccatnum"]))
maxph_p = collections.defaultdict(int)
for r in cvp:
    if r["ccatnum"].isdigit():
        maxph_p[r["gene"]] = max(maxph_p[r["gene"]], int(r["ccatnum"]))


def cls(sym):
    s = ann.get(sym, {}).get("state", "UNMAPPED")
    if s in ("ABSENT", "TESTED_NULL"):
        return s, "CLASS_A"
    if sym == "IL6R":
        return s, "CLASS_B"
    return s, ""


rec = {g: dict(gene=g, state=cls(g)[0], klass=cls(g)[1], ph=maxph.get(g, 0)) for g in te}
st = collections.Counter(v["state"] for v in rec.values())
N = len(rec)
for k in ("ABSENT", "TESTED_NULL", "INSTRUMENTED", "UNMAPPED"):
    print("  %-13s %4d  %s of %d" % (k, st[k], pct(st[k], N), N))
nA = st["ABSENT"] + st["TESTED_NULL"]
print("  CLASS_A       %4d  %s" % (nA, pct(nA, N)))
print("  CLASS_B       %4d  %s  (assigned by name: IL6R only)"
      % (sum(1 for v in rec.values() if v["klass"] == "CLASS_B"),
         pct(sum(1 for v in rec.values() if v["klass"] == "CLASS_B"), N)))
print("  CLASS_A as %% of MAPPED targets (%d): %s" % (N - st["UNMAPPED"],
                                                    pct(nA, N - st["UNMAPPED"])))
print("  CLASS_A + CLASS_B: %d  %s" % (nA + 1, pct(nA + 1, N)))
print("  UNMAPPED genes: %s" % ", ".join(sorted(g for g, v in rec.items()
                                                 if v["state"] == "UNMAPPED")))

# cross-check against the shipped atlas, row by row
atlas = list(csv.DictReader(open(DEL + r"\cv\cardiovascular_atlas.csv", encoding="utf-8")))
print("\nshipped atlas rows: %d, columns: %d" % (len(atlas), len(atlas[0])))
bad = 0
for a in atlas:
    v = rec.get(a["gene"])
    if (v is None or v["state"] != a["state"] or v["klass"] != a["untestability_class"]
            or str(v["ph"]) != a["max_phase_num"]):
        bad += 1
        print("   MISMATCH", a["gene"], a["state"], a["untestability_class"],
              a["max_phase_num"], "vs", v)
print("atlas rows disagreeing with recomputation (state/class/max phase): %d" % bad)
print("atlas genes == recomputed CV-EXTENDED genes: %s" % ({a["gene"] for a in atlas} == te))
for sub in ("cv", "supplementary"):
    for f in ("cardiovascular_atlas.csv", "classB_upper_bound.csv"):
        print("  sha256 %-14s %-26s %s" % (sub, f, sha(DEL + "\\" + sub + "\\" + f)))

# ---------------------------------------------------------------- V04 class x phase
hdr("V04", "CLASS x HIGHEST PHASE REACHED (target level; phase = max ccatnum over the "
           "target's CV-EXTENDED pairs)")
cats = ["ABSENT", "TESTED_NULL", "CLASS_A", "INSTR_other", "CLASS_B", "UNMAPPED", "ALL"]


def cat_of(v):
    if v["state"] == "UNMAPPED":
        return "UNMAPPED"
    if v["klass"] == "CLASS_B":
        return "CLASS_B"
    if v["state"] == "INSTRUMENTED":
        return "INSTR_other"
    return v["state"]


tab = collections.defaultdict(collections.Counter)
for v in rec.values():
    c = cat_of(v)
    tab[v["ph"]][c] += 1
    tab[v["ph"]]["ALL"] += 1
    if c in ("ABSENT", "TESTED_NULL"):
        tab[v["ph"]]["CLASS_A"] += 1
print("%-12s" % "phase" + "".join("%12s" % c for c in cats) + "   classA_%_of_phase")
out_rows = []
for ph in sorted(tab):
    t = tab[ph]
    print("%-12s" % PHASE[ph] + "".join("%12d" % t[c] for c in cats)
          + "   %s" % pct(t["CLASS_A"], t["ALL"]))
    out_rows.append(dict(phase=PHASE[ph], **{c: t[c] for c in cats},
                         classA_pct_of_phase=round(100 * t["CLASS_A"] / t["ALL"], 2)))
tot = collections.Counter()
for ph in tab:
    tot.update(tab[ph])
print("%-12s" % "TOTAL" + "".join("%12d" % tot[c] for c in cats)
      + "   %s" % pct(tot["CLASS_A"], tot["ALL"]))
with open(HERE + r"\class_by_phase.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), lineterminator="\n")
    w.writeheader(); w.writerows(out_rows)
print("raw ccat labels among CV-EXTENDED pairs: %s"
      % dict(collections.Counter((r["ccatnum"], r["ccat"]) for r in cve)))
print("CV-EXTENDED pairs with non-numeric ccatnum: %d"
      % sum(1 for r in cve if not r["ccatnum"].isdigit()))
# trend
a = [[tab[ph]["CLASS_A"], tab[ph]["INSTR_other"] + tab[ph]["CLASS_B"]] for ph in sorted(tab)]
chi2, pchi, dof, _ = stats.chi2_contingency(a)
print("descriptive only: chi-square ClassA-vs-instrumented across phases: chi2=%.2f dof=%d "
      "p=%.3g" % (chi2, dof, pchi))

# ---------------------------------------------------------------- V05 launched hook
hdr("V05", "HOOK: launched / Phase III CV targets by class")
L = [v for v in rec.values() if v["ph"] == 5]
P3 = [v for v in rec.values() if v["ph"] == 4]
for name, S in (("LAUNCHED", L), ("PHASE III", P3), ("PHASE III or LAUNCHED", L + P3)):
    c = collections.Counter(cat_of(v) for v in S)
    na = c["ABSENT"] + c["TESTED_NULL"]
    mapped = len(S) - c["UNMAPPED"]
    print("%s: N=%d  %s" % (name, len(S), dict(c)))
    print("   CLASS_A = %d  -> %s of all %d ; %s of the %d mapped"
          % (na, pct(na, len(S)), len(S), pct(na, mapped), mapped))
    print("   ABSENT only = %d -> %s of all" % (c["ABSENT"], pct(c["ABSENT"], len(S))))
    if name == "LAUNCHED":
        print("   CLASS_A beyond PCSK9 = %d -> %s of all %d launched"
              % (na - 1, pct(na - 1, len(S)), len(S)))
        print("   ABSENT beyond PCSK9  = %d" % (c["ABSENT"] - 1))
        print("   CLASS_A + CLASS_B(IL6R launched? %s)" % any(v["gene"] == "IL6R" for v in S))
    print("   unmapped genes: %s" % sorted(v["gene"] for v in S if v["state"] == "UNMAPPED"))
print("IL6R max CV phase: %s" % PHASE[rec["IL6R"]["ph"]])
# CV-PRIMARY robustness of the hook
Lp = [g for g in tp if maxph_p.get(g, 0) == 5]
cA = [g for g in Lp if cls(g)[1] == "CLASS_A"]
print("CV-PRIMARY: launched targets %d, CLASS_A %d (%s); not in CV-PRIMARY launched "
      "Class A but in EXTENDED: %s"
      % (len(Lp), len(cA), pct(len(cA), len(Lp)),
         sorted({v["gene"] for v in L if v["klass"] == "CLASS_A"} - set(cA))))

# detail rows for HEADLINE.md
t1 = {r["gene"]: r for r in csv.DictReader(
    open(DEL + r"\cv\table1_classA_launched_phase3.csv", encoding="utf-8"))}
arow = {a["gene"]: a for a in atlas}
det = []
for v in sorted(L + P3, key=lambda x: (-x["ph"], x["gene"])):
    if v["klass"] not in ("CLASS_A", "CLASS_B"):
        continue
    g = v["gene"]
    top = [r for r in cve if r["gene"] == g and r["ccatnum"] == str(v["ph"])]
    det.append(dict(
        gene=g, phase=PHASE[v["ph"]], klass=v["klass"], state=v["state"],
        n_cv_indications_at_max_phase=len(top),
        cv_indications_at_max_phase="; ".join(sorted(r["indication_mesh_term"] for r in top)),
        in_cv_primary=g in tp,
        year_launch="; ".join(sorted({r.get("year_launch", "") for r in top} - {""})),
        hcat_acat="; ".join(sorted({"%s|%s" % (r["hcat"], r["acat"]) for r in top})),
        table1_example_agent=t1.get(g, {}).get("example_agent", "NOT IN table1"),
        table1_indication=t1.get(g, {}).get("indication", ""),
        blood_tpm=arow[g]["gtex_blood_tpm"], heart_LV_tpm=arow[g]["gtex_heart_LV_tpm"],
        liver_tpm=arow[g]["gtex_liver_tpm"], required_modality=arow[g]["required_modality"]))
with open(HERE + r"\headline_targets.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(det[0].keys()), lineterminator="\n")
    w.writeheader(); w.writerows(det)
for d in det:
    print("  %-8s %-9s %-7s %-12s agent=%-20s yr=%-10s h|a=%s\n        CV indications at "
          "max phase (%d): %s"
          % (d["gene"], d["phase"], d["klass"], d["state"], d["table1_example_agent"],
             d["year_launch"], d["hcat_acat"], d["n_cv_indications_at_max_phase"],
             d["cv_indications_at_max_phase"]))
print("table1 genes == recomputed Class A launched+PhIII genes: %s"
      % (set(t1) == {d["gene"] for d in det if d["klass"] == "CLASS_A"}))
print("table1 rows whose single listed indication is NOT the only CV indication at max "
      "phase: %d of %d" % (sum(1 for d in det if d["n_cv_indications_at_max_phase"] > 1
                               and d["gene"] in t1), len(t1)))
print("table1 rows with example_agent 'not in Open Targets': %s"
      % sorted(g for g, r in t1.items() if r["example_agent"] == "not in Open Targets"))

# ---------------------------------------------------------------- V06 class x indication
hdr("V06", "CLASS x CV INDICATION (MeSH term level - the only indication grouping that "
           "exists in the repo)")
ind = collections.defaultdict(lambda: collections.Counter())
for r in cve:
    v = rec[r["gene"]]
    k = r["indication_mesh_term"]
    ind[k]["targets"] += 1          # ti_uid-unique so 1 row == 1 target for this indication
    ind[k][cat_of(v)] += 1
    if r["ccatnum"] == "5":
        ind[k]["launched_pairs"] += 1
        if v["klass"] == "CLASS_A":
            ind[k]["launched_classA"] += 1
rows6 = []
for k, c in ind.items():
    na = c["ABSENT"] + c["TESTED_NULL"]
    rows6.append(dict(indication_mesh_term=k,
                      cv_primary=any(r["cvp"] for r in cve if r["indication_mesh_term"] == k),
                      n_targets=c["targets"], absent=c["ABSENT"], tested_null=c["TESTED_NULL"],
                      class_A=na, class_A_pct=round(100 * na / c["targets"], 1),
                      instrumented=c["INSTR_other"] + c["CLASS_B"], unmapped=c["UNMAPPED"],
                      launched_pairs=c["launched_pairs"],
                      launched_pairs_classA=c["launched_classA"]))
rows6.sort(key=lambda x: -x["n_targets"])
with open(HERE + r"\class_by_indication_mesh.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows6[0].keys()), lineterminator="\n")
    w.writeheader(); w.writerows(rows6)
print("distinct CV-EXTENDED MeSH indications: %d  (full table -> class_by_indication_mesh.csv)"
      % len(rows6))
print("%-44s %5s %5s %6s %7s %5s %5s" % ("indication (top 30 by targets)", "N", "A", "A%",
                                          "instr", "L", "L_A"))
for r in rows6[:30]:
    print("%-44s %5d %5d %6s %7d %5d %5d"
          % (r["indication_mesh_term"][:44], r["n_targets"], r["class_A"],
             "%.1f" % r["class_A_pct"], r["instrumented"], r["launched_pairs"],
             r["launched_pairs_classA"]))
print("pair-level: CV-EXTENDED pairs whose target is Class A: %d of %d (%s)"
      % (sum(1 for r in cve if rec[r["gene"]]["klass"] == "CLASS_A"), len(cve),
         pct(sum(1 for r in cve if rec[r["gene"]]["klass"] == "CLASS_A"), len(cve))))
print("pair-level: pairs whose target is INSTRUMENTED: %d"
      % sum(1 for r in cve if rec[r["gene"]]["state"] == "INSTRUMENTED"))
print("any sub-area column (lipids/HF/arrhythmia/...) in atlas header: %s"
      % [c for c in atlas[0] if "area" in c.lower() or "indication" in c.lower()])

# ---------------------------------------------------------------- V07 enrichment
hdr("V07", "CLASS A ENRICHMENT for tissue-restricted expression (from shipped atlas)")
inh = [a for a in atlas if a["tissue_restricted"] in ("YES", "NO")]
q = lambda s, t: sum(1 for a in inh if a["state"] == s and a["tissue_restricted"] == t)
t22 = (q("ABSENT", "YES"), q("ABSENT", "NO"), q("INSTRUMENTED", "YES"),
       q("INSTRUMENTED", "NO"))
print("2x2 [[ABSENT TR, ABSENT low],[INSTR TR, INSTR low]] = %s" % (t22,))
print("OR %.3f  95%% CI [%.3f, %.3f]  Fisher p %.4g" % or_ci(*t22))
print("CV targets with tissue_restricted in YES/NO: %d of %d" % (len(inh), len(atlas)))
for lab, key in (("blood", "instrumented_blood"), ("heart LV", "instrumented_heart_LV"),
                 ("heart atrial", "instrumented_heart_atrial"),
                 ("artery aorta", "instrumented_artery_aorta"),
                 ("liver", "instrumented_liver"), ("cardiac union", "instrumented_cardiac_union")):
    n = sum(1 for a in atlas if a[key] == "True")
    print("  yield %-14s %4d / %d = %s" % (lab, n, len(atlas), pct(n, len(atlas))))
m = collections.Counter(a["required_modality"].split(" (")[0] for a in atlas
                        if a["untestability_class"] == "CLASS_A")
print("required modality, Class A: %s" % dict(m.most_common()))
print("Class A launched targets with heart LV TPM >= 1 -> labelled 'cardiac eQTL' while "
      "liver TPM is >10x higher:")
for a in atlas:
    if (a["untestability_class"] == "CLASS_A" and a["required_modality"] == "cardiac eQTL"
            and a["gtex_liver_tpm"] and a["gtex_heart_LV_tpm"]
            and float(a["gtex_liver_tpm"]) > 10 * max(float(a["gtex_heart_LV_tpm"]), 1e-9)
            and a["max_phase_num"] in ("4", "5")):
        print("    %-7s heartLV %.1f  liver %.1f  phase %s"
              % (a["gene"], float(a["gtex_heart_LV_tpm"]), float(a["gtex_liver_tpm"]),
                 a["max_clinical_phase"]))

# ---------------------------------------------------------------- V08 controls
hdr("V08", "POSITIVE CONTROLS")
print("-- Deliverable 1, GTEx v10 liver (d1/egene_status.json):")
for k, v in json.load(open(DEL + r"\d1\egene_status.json")).items():
    print("   %-6s top nominal cis p = %-12g eGene(FDR0.05) = %s ; resource eGenes %d / %d"
          % (v["symbol"], v["top_nominal_pvalue"], v["is_egene_fdr005"],
             v["n_egenes_fdr005"], v["n_genes_tested"]))
print("-- Deliverable 1c, eQTLGen (d1c/d1c_summary.json):")
for d in json.load(open(DEL + r"\d1c\d1c_summary.json")):
    mr = d["mr"]
    print("   %s VERDICT=%s  hits=%d  n_indep=%d (%s)" % (
        d["pair"], d["VERDICT"], d["n_eqtlgen_hits_p_lt_5e-8"],
        d["clumping"]["n_independent"], ",".join(d["clumping"]["independent_snps"])))
    print("      clump lead %s p_exp=%g F=%.1f" % (d["lead_snp"]["rsid"],
                                                 d["lead_snp"]["p_exp"], d["lead_snp"]["F"]))
    print("      MR lead %s lead_p_exposure=%g F=%.1f  beta=%.5f se=%.6f p=%.3g  dir=%s "
          "expected=%s  overlap=%d" % (mr["lead_snp"], mr["lead_p_exposure"], mr["F"],
                                       mr["beta"], mr["se"], mr["p"], mr["direction"],
                                       mr["expected_direction"],
                                       d["n_overlapping_after_harmonization"]))
    print("      pass flags: %s" % d["pass"])
for pair in ("HMGCR_LDL", "IL6R_CRP"):
    rows = tsv(DEL + r"\d1c\%s_harmonized_merge.tsv" % pair)
    lead = min(rows, key=lambda r: float(r["p_exp"]))
    b = float(lead["beta_out_harm"]) / float(lead["beta_exp"])
    se = abs(float(lead["se_out"]) / float(lead["beta_exp"]))
    z = b / se
    print("   independent Wald recompute %s: rows=%d lead(min p_exp in merge)=%s p_exp=%s "
          "F=%.1f -> beta=%.5f se=%.6f z=%.2f p=%.3g"
          % (pair, len(rows), lead["rsid"], lead["p_exp"], float(lead["F"]), b, se, z,
             2 * stats.norm.sf(abs(z))))
    print("      harmonization: %s" % dict(collections.Counter(r["harmonization"] for r in rows)))
    sig = [r for r in rows if float(r["p_exp"]) < 1e-10]
    pos = sum(1 for r in sig if float(r["beta_out_harm"]) / float(r["beta_exp"]) > 0)
    print("      variants p_exp<1e-10: %d, positive Wald %d, negative %d"
          % (len(sig), pos, len(sig) - pos))
print("-- atlas rows for the controls:")
for g in ("HMGCR", "IL6R", "PCSK9"):
    a = arow.get(g)
    print("   %-6s state=%s class=%s n_indep=%s best_eqtl_p=%s maxphase=%s blood_tpm=%s "
          "liver_tpm=%s" % (g, a["state"], a["untestability_class"],
                            a["n_independent_instruments"], a["best_eqtl_p"],
                            a["max_clinical_phase"], a["gtex_blood_tpm"], a["gtex_liver_tpm"]))
uni = open(DEL + r"\d1b\eqtlgen_tested_universe.tsv", encoding="utf-8").read()
print("-- PCSK9 in d1b/eqtlgen_tested_universe.tsv: ENSG00000169174 %s ; symbol token PCSK9 %s"
      % ("ENSG00000169174" in uni, bool(re.search(r"\bPCSK9\b", uni))))
print("   PCSK4/5/6/7 present: %s"
      % {s: bool(re.search(r"\b%s\b" % s, uni)) for s in ("PCSK4", "PCSK5", "PCSK6", "PCSK7")})
print("   (full raw-file scan: see pcsk9_scan.txt)")

# ---------------------------------------------------------------- V09 Class C
hdr("V09", "CLASS C refutation test (d1h/part3_headtohead.tsv)")
h2h = tsv(DEL + r"\d1h\part3_headtohead.tsv")
pr = [r for r in h2h if r["PP4_eqtlgen"] and r["PP4_gtexblood"]]
print("loci attempted %d ; paired %d ; unpaired notes: %s"
      % (len(h2h), len(pr), dict(collections.Counter(
          r["note"].split(":")[0] for r in h2h if not (r["PP4_eqtlgen"] and r["PP4_gtexblood"])))))
e = [float(r["PP4_eqtlgen"]) for r in pr]
g_ = [float(r["PP4_gtexblood"]) for r in pr]
hi = sum(1 for x, y in zip(e, g_) if x > y)
lo = sum(1 for x, y in zip(e, g_) if x < y)
import statistics as S
d = [math.log10(x) - math.log10(y) for x, y in zip(e, g_)]
print("PP4 higher in eQTLGen %d ; lower %d ; ties %d" % (hi, lo, len(pr) - hi - lo))
print("median PP4 GTEx blood %.5f ; eQTLGen %.5f" % (S.median(g_), S.median(e)))
print("median log10 difference (eQTLGen - GTEx): %+.4f" % S.median(d))
w = stats.wilcoxon(d)
print("Wilcoxon signed-rank on log10(PP4) diff: W=%.1f p=%.4f (two-sided)"
      % (w.statistic, w.pvalue))
print("exact binomial two-sided (%d of %d): p=%.4f"
      % (hi, hi + lo, stats.binomtest(hi, hi + lo, 0.5).pvalue))
ez = [abs(float(r["eqtlgen_lead_z"])) for r in pr]
gz = [abs(float(r["gtex_lead_z"])) for r in pr]
print("median |lead z| eQTLGen %.1f ; GTEx %.1f" % (S.median(ez), S.median(gz)))
rho = stats.spearmanr(ez, d)
print("Spearman(log10 diff, eQTLGen |lead z|): rho=%.3f p=%.5f" % (rho.statistic, rho.pvalue))
print("max PP4 in either arm: eQTLGen %.4f ; GTEx %.4f ; loci with PP4>0.5 in either arm: %d"
      % (max(e), max(g_), sum(1 for x, y in zip(e, g_) if max(x, y) > 0.5)))
print("distinct outcome GWAS among paired loci: %d" % len({r["gwas_acc"] for r in pr}))
print("paired loci whose gene is in the 630 CV atlas: %d of %d"
      % (sum(1 for r in pr if r["gene"] in te), len(pr)))

# ---------------------------------------------------------------- V10 method availability
hdr("V10", "METHOD AVAILABILITY (CV-restricted)")
r15 = tsv(DEL + r"\d1j\r15_results.tsv")
cvr = [r for r in r15 if r["gene"] in te]
sig = [r for r in cvr if r["gwas_genome_wide_sig"] == "True"]
use = [r for r in sig if r["susie_PP4"]]
print("R15 loci screened (all) %d ; genome-wide sig on outcome side %d"
      % (len(r15), sum(1 for r in r15 if r["gwas_genome_wide_sig"] == "True")))
print("CV-atlas genes among R15 loci %d ; outcome-significant %d (%s) ; with susie PP4 %d (%s)"
      % (len(cvr), len(sig), ", ".join(r["gene"] + "/" + r["mesh_term"] for r in sig),
         len(use), ", ".join(r["gene"] for r in use)))
print("   of those %d CV-gene loci, the locus's OWN indication is a CV-EXTENDED indication: %d"
      % (len(cvr), sum(1 for r in cvr if any(
          k in r["mesh_term"].lower() for k in TOK) or any(
          x["cvp"] for x in frozen if x["indication_mesh_term"] == r["mesh_term"]))))
pc = tsv(DEL + r"\d1f\partc_pairs.tsv")
res = tsv(DEL + r"\d1f\partc_results.tsv")
skp = tsv(DEL + r"\d1f\partc_skipped_accessions.tsv")
got = {(r["gene"], r["mesh_id"], r["gwas_acc"]) for r in res}
print("Part C: qualifying pairs %d ; result rows %d ; skipped accessions %d ; pairs lost "
      "(sum n_pairs_lost) %d ; lost (pairs not in results) %d -> %s"
      % (len(pc), len(res), len(skp), sum(int(r["n_pairs_lost"]) for r in skp),
         sum(1 for p in pc if (p["gene"], p["mesh_id"], p["gwas_acc"]) not in got),
         pct(sum(1 for p in pc if (p["gene"], p["mesh_id"], p["gwas_acc"]) not in got),
             len(pc))))
print("   skip reasons recorded: %s" % dict(collections.Counter(r["reason"] for r in skp)))
is_cv = lambda mid, term: (areas.get(topl.get(mid)) == "cardiovascular"
                           or any(k in term.lower() for k in TOK))
cvpc = [p for p in pc if is_cv(p["mesh_id"], p["mesh_term"])]
lost = [p for p in cvpc if (p["gene"], p["mesh_id"], p["gwas_acc"]) not in got]
print("CV Part C pairs %d (targets %d) ; lost %d (%s) ; accessions: %s"
      % (len(cvpc), len({p["gene"] for p in cvpc}), len(lost), pct(len(lost), len(cvpc)),
         dict(collections.Counter(p["gwas_acc"] for p in lost))))
# broad denominator, logic copied read-only from src/a3_cascade.py
harm = set()
for l in open(RAW + r"\harmonised_list.txt"):
    mm = re.search(r"/(GCST\d+)/harmonised/", l)
    if mm:
        harm.add(mm.group(1))
ANC = ["African", "East Asian", "South Asian", "Hispanic", "Latin", "Native",
       "Asian unspecified", "Oceanian", "Greater Middle Eastern", "Aboriginal",
       "Sub-Saharan", "NR"]
best = {}
with open(RAW + r"\gwas-catalog-download-studies-v1.0.3.1.txt", encoding="utf-8",
          errors="replace") as fh:
    rd = csv.reader(fh, delimiter="\t")
    hh = next(rd)
    it, ia, isz = hh.index("MAPPED_TRAIT"), hh.index("STUDY ACCESSION"), hh.index(
        "INITIAL SAMPLE SIZE")
    for f in rd:
        if len(f) <= it or f[ia] not in harm:
            continue
        s = f[isz]
        if not (("European" in s) and not any(a_ in s for a_ in ANC)):
            continue
        n = max([int(x.replace(",", "")) for x in re.findall(r"[\d,]{4,}", s)] or [0])
        for t in f[it].split(","):
            t = t.strip().lower()
            if t and (t not in best or n > best[t][1]):
                best[t] = (f[ia], n)
skipped = {r["gwas_acc"] for r in skp}
probed_ok = {r["gwas_acc"] for r in res}
for name, pairs in (("CV-EXTENDED", cve), ("CV-PRIMARY", cvp)):
    mp = [r for r in pairs if r["indication_mesh_term"].strip().lower() in best]
    acc = [best[r["indication_mesh_term"].strip().lower()][0] for r in mp]
    nl = sum(1 for a_ in acc if a_ in skipped)
    nok = sum(1 for a_ in acc if a_ in probed_ok)
    print("%s: pairs mappable to a EUR-only harmonised GWAS %d ; accession in skip list %d "
          "(%s) ; accession produced Part C results %d ; accession NEVER PROBED in Part C %d"
          % (name, len(mp), nl, pct(nl, len(mp)), nok, len(mp) - nl - nok))

# ---------------------------------------------------------------- V11 Class B bound
hdr("V11", "CLASS B upper bound (cv/classB_upper_bound.csv)")
cb = list(csv.DictReader(open(DEL + r"\cv\classB_upper_bound.csv", encoding="utf-8")))
ins = [a for a in atlas if a["state"] == "INSTRUMENTED"]
print("flagged rows %d of %d instrumented (%s) ; cv_relevant_trait True %d (%s)"
      % (len(cb), len(ins), pct(len(cb), len(ins)),
         sum(1 for r in cb if r["cv_relevant_trait"] == "True"),
         pct(sum(1 for r in cb if r["cv_relevant_trait"] == "True"), len(ins))))
insL = [a for a in ins if a["max_phase_num"] == "5"]
fl = {r["gene"]: r for r in cb}
print("instrumented launched %d ; flagged any-trait %d ; flagged CV-trait %d"
      % (len(insL), sum(1 for a in insL if a["gene"] in fl),
         sum(1 for a in insL if a["gene"] in fl and fl[a["gene"]]["cv_relevant_trait"] == "True")))
print("IL6R in bound file: %s" % ("IL6R" in fl))
print("\nDONE")
