#!/usr/bin/env python
"""Corrected classification for the manuscript: primary CV set + sensitivity sets.

    python manuscript_prep/analysis/classify.py

READ-ONLY on the project. Writes only under manuscript_prep/analysis/out/.
Network: HGNC REST for unmapped symbols (cached in verify/cache/hgnc/).

DECISIONS IMPLEMENTED (manuscript_prep/DECISIONS_LOG.md, 2026-10-10)
  CV sets
    PRIMARY  every indication Minikel tags 'cardiovascular' in indic.tsv `areas`,
             plus the 8 lipid/atherosclerosis tokens                       (93 terms)
    SENS1    original atlas set: top-level branch C14 + lipid tokens        (44 terms)
    SENS2    SENS1 + C23-filed indications tagged cardiovascular            (62 terms)
    C14ONLY  top-level branch C14 only (the old 'CV-PRIMARY')               (37 terms)
  Status of a target (window W = gene body +- 100 kb is PRIMARY; 500 kb, 1 Mb sensitivity)
    INSTRUMENTED               >=1 eQTLGen cis-eQTL, p < 5e-8, inside W
    NO_EQTL_DATA               gene not among the 19,250 genes eQTLGen tested
    NO_QUALIFYING_INSTRUMENT   tested, no p < 5e-8 variant inside W
    UNRESOLVED                 symbol has no unique GRCh37.87 gene record after HGNC rule
    'uninstrumentable' = NO_EQTL_DATA (lower bound)  ..  + NO_QUALIFYING_INSTRUMENT (upper)
    NO_EQTL_DATA has two named sub-types: 'autosomal gene not tested' and
    'X-linked: chromosome not analysed by eQTLGen'. The AUTOSOMAL-ONLY lower bound is
    autosomal no-eQTL-data targets / targets not on X, Y or MT, and is reported next to
    the overall lower bound everywhere.
  Drug verification (targets at Phase III or launched; ChEMBL 37)
    PERMISSIVE  an approved drug with any link to the gene has an approved indication
                under MeSH C14 or a lipid term
    STRICT      the link is direct (single protein, specific complex, or mRNA map) AND the
                approved indication sits only in cardiovascular MeSH branches
                (rule text in verify/build_chembl_links.py)
    Drug verification was completed on 2026-10-11 (drug_links/REVIEW_DECISION.md).
  Example eligibility   a target may be named as an example only if it is strict-verified
                        and holds the same phase in SENS2 (i.e. not through one of the 31
                        indications added outside C14 / signs-symptoms).
  Gene records   shipped mapping + HGNC Ensembl-ID rule for every unmapped symbol
  Launched       phase(target, S) = max ccatnum over the target's pairs in S; launched = 5
  Class B        dropped; IL6R is INSTRUMENTED
"""
import collections
import csv
import gzip
import json
import math
import os
import re
import statistics
import sys
import urllib.request
import zipfile

from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
RAW = ROOT + r"\data\raw"
D1B = ROOT + r"\deliverables\d1b"
D1F = ROOT + r"\deliverables\d1f"
HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
OUT = os.path.join(HERE, "out")
HCACHE = os.path.join(MP, "verify", "cache", "hgnc")
os.makedirs(OUT, exist_ok=True)
os.makedirs(HCACHE, exist_ok=True)

P_INSTR = 5e-8
WINDOWS = {"100kb": 100_000, "500kb": 500_000, "1Mb": 1_000_000}
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
PH = {0: "", 1: "Preclinical", 2: "Phase I", 3: "Phase II", 4: "Phase III", 5: "Launched"}
CARDIAC = [("QTD000256", "heart_LV"), ("QTD000251", "heart_atrial"),
           ("QTD000131", "artery_aorta"), ("QTD000266", "liver")]
RESTRICTED = {"Tissue enriched", "Group enriched", "Tissue enhanced"}
NOT_RESTRICTED = {"Low tissue specificity"}
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# ---- Table 3 indication areas (approved 2026-10-10; every placement is logged in
# ---- DECISIONS_LOG.md). Anything not listed here is 'Other'; the script asserts that
# ---- every PRIMARY term is either listed or in OTHER_TERMS, so nothing defaults silently.
LIP, IHD, HF, ARR, HTN, PHT, THR, CBV, OTH, UNS = (
    "Lipids / atherosclerosis", "Ischaemic heart disease", "Heart failure / cardiomyopathy",
    "Arrhythmia", "Hypertension", "Pulmonary hypertension", "Thrombosis",
    "Cerebrovascular", "Other", "Unspecified CV")
AREA_ORDER = [LIP, IHD, HF, ARR, HTN, PHT, THR, CBV, OTH, UNS]
AREA = {}
for a_, terms_ in {
    LIP: ["Atherosclerosis", "Hyperlipidemias", "Hypercholesterolemia", "Hypertriglyceridemia",
          "Hyperlipoproteinemia Type II", "Dyslipidemias", "Hyperlipoproteinemia Type I",
          "Hypoalphalipoproteinemias"],
    IHD: ["Myocardial Infarction", "Acute Coronary Syndrome", "Angina Pectoris",
          "Angina, Unstable", "Coronary Restenosis", "Coronary Thrombosis"],
    HF: ["Heart Failure", "Cardiomyopathies", "Cardiomyopathy, Hypertrophic",
         "Endomyocardial Fibrosis", "Edema, Cardiac", "Myocarditis", "Cardiomyopathy, Dilated",
         "Diabetic Cardiomyopathies"],
    ARR: ["Atrial Fibrillation", "Arrhythmias, Cardiac", "Ventricular Fibrillation",
          "Tachycardia, Ventricular", "Tachycardia, Supraventricular", "Long QT Syndrome",
          "Tachycardia", "Bradycardia"],
    HTN: ["Hypertension"],
    PHT: ["Hypertension, Pulmonary"],
    THR: ["Thrombosis", "Venous Thrombosis", "Pulmonary Embolism", "Thrombophlebitis"],
    CBV: ["Stroke", "Brain Ischemia", "Cerebral Hemorrhage", "Subarachnoid Hemorrhage",
          "Cerebral Infarction", "Intracranial Thrombosis", "Vasospasm, Intracranial",
          "CADASIL", "Dementia, Vascular", "Cerebral Amyloid Angiopathy"],
    UNS: ["Cardiovascular Diseases"],
    OTH: ["Peripheral Vascular Diseases", "Intermittent Claudication", "Reperfusion Injury",
          "Vasculitis", "Hypotension", "Retinal Vein Occlusion", "Retinal Artery Occlusion",
          "Optic Neuropathy, Ischemic", "Retinal Telangiectasis", "Raynaud Disease",
          "Vascular Diseases", "Venous Insufficiency", "Aneurysm", "Aortic Valve Stenosis",
          "Capillary Leak Syndrome", "Thromboangiitis Obliterans", "Erythromelalgia",
          "Endocarditis", "Pericarditis",
          # filed outside C14/C23 by Minikel
          "Diabetic Retinopathy", "Diabetic Foot", "Behcet Syndrome",
          "Telangiectasia, Hereditary Hemorrhagic", "Angioedemas, Hereditary",
          "Ehlers-Danlos Syndrome", "Granulomatosis with Polyangiitis", "Varicose Ulcer",
          "Churg-Strauss Syndrome", "Takayasu Arteritis", "Tetralogy of Fallot",
          "Fabry Disease", "Microscopic Polyangiitis", "Turner Syndrome", "MELAS Syndrome",
          "Mucocutaneous Lymph Node Syndrome", "Noonan Syndrome", "Ductus Arteriosus, Patent",
          "von Hippel-Lindau Disease", "Angioedema", "Glycogen Storage Disease Type IIb",
          "Sarcoglycanopathies", "Marfan Syndrome", "Ataxia Telangiectasia",
          "Heart Defects, Congenital", "Arteriovenous Malformations", "Alagille Syndrome"],
}.items():
    for t_ in terms_:
        assert t_ not in AREA, t_
        AREA[t_] = a_


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def pct(a, b):
    return round(100.0 * a / b, 1) if b else None


def or_ci(a, b, c, d):
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a + .5, b + .5, c + .5, d + .5)
    lor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return dict(table=[[a, b], [c, d]], odds_ratio=orr, ci_lo=math.exp(lor - 1.96 * se),
                ci_hi=math.exp(lor + 1.96 * se), p_fisher=p)


def hdr(tag, title):
    print("\n" + "=" * 86 + "\n[%s] %s\n" % (tag, title) + "=" * 86)


K = {}          # key numbers -> out/key_numbers.json

# ============================================================== K01 flow and sets
hdr("K01", "FLOW AND CV SETS")
ind = {r["indication_mesh_id"]: r for r in load(RAW + r"\minikel_indic.tsv")}
topl = {r["indication_mesh_id"]: r["topl"] for r in load(RAW + r"\minikel_indic_topl_match.tsv")}
raw_lines = open(RAW + r"\pp.tsv", encoding="utf-8").read().rstrip("\n").split("\n")
frozen = load(ROOT + r"\data\derived\pp_dedup_frozen.tsv")
iu = raw_lines[0].split("\t").index("ti_uid")
n_raw = len(raw_lines) - 1
n_uid = len({l.split("\t")[iu] for l in raw_lines[1:]})


def tagged(i):
    return "cardiovascular" in [a.strip() for a in ind[i]["areas"].split(",")]


def lipid(i):
    return any(k in ind[i]["indication_mesh_term"].lower() for k in TOK)


SETDEF = collections.OrderedDict([
    ("PRIMARY", lambda i: tagged(i) or lipid(i)),
    ("SENS1", lambda i: topl[i] == "C14" or lipid(i)),
    ("SENS2", lambda i: topl[i] == "C14" or lipid(i) or (topl[i] == "C23" and tagged(i))),
    ("C14ONLY", lambda i: topl[i] == "C14"),
])
SETLABEL = {"PRIMARY": "Primary: all indications Minikel tags cardiovascular + lipid tokens",
            "SENS1": "Sensitivity 1: original atlas set (MeSH C14 branch + lipid tokens)",
            "SENS2": "Sensitivity 2: sensitivity 1 + signs/symptoms terms tagged cardiovascular",
            "C14ONLY": "Sensitivity 3: MeSH C14 branch only"}
used = {r["indication_mesh_id"] for r in frozen}
setids = {k: {i for i in used if fn(i)} for k, fn in SETDEF.items()}
# every C14-top-level indication must also carry Minikel's cardiovascular tag
assert all(tagged(i) for i in setids["C14ONLY"]), "C14 term without cardiovascular tag"
assert setids["SENS1"] <= setids["SENS2"] <= setids["PRIMARY"]
pairs = {k: [r for r in frozen if r["indication_mesh_id"] in setids[k]] for k in SETDEF}
tset = {k: {r["gene"] for r in pairs[k]} for k in SETDEF}
phase = {k: collections.defaultdict(int) for k in SETDEF}
for k in SETDEF:
    for r in pairs[k]:
        phase[k][r["gene"]] = max(phase[k][r["gene"]], int(r["ccatnum"]))
K["flow"] = dict(raw_rows=n_raw, distinct_ti_uid=n_uid, frozen_pairs=len(frozen),
                 targets=len({r["gene"] for r in frozen}), indications=len(used))
print("raw rows %d -> distinct ti_uid %d -> frozen pairs %d ; targets %d ; indications %d"
      % (n_raw, n_uid, len(frozen), K["flow"]["targets"], len(used)))
K["sets"] = {}
for k in SETDEF:
    K["sets"][k] = dict(label=SETLABEL[k], terms=len(setids[k]), pairs=len(pairs[k]),
                        targets=len(tset[k]))
    print("%-8s terms %3d  pairs %5d  targets %4d   %s"
          % (k, len(setids[k]), len(pairs[k]), len(tset[k]), SETLABEL[k]))
missing_area = sorted(ind[i]["indication_mesh_term"] for i in setids["PRIMARY"]
                      if ind[i]["indication_mesh_term"] not in AREA)
assert not missing_area, "terms without an area: %s" % missing_area
assert set(AREA) == {ind[i]["indication_mesh_term"] for i in setids["PRIMARY"]}

# terms flagged as in the set by rule but not primarily cardiac (none is removed)
NOT_PRIMARILY_CARDIAC = {
    "Diabetic Retinopathy", "Diabetic Foot", "Angioedemas, Hereditary", "Angioedema",
    "Ehlers-Danlos Syndrome", "Fabry Disease", "Turner Syndrome", "MELAS Syndrome",
    "Noonan Syndrome", "Marfan Syndrome", "von Hippel-Lindau Disease",
    "Ataxia Telangiectasia", "Alagille Syndrome", "Sarcoglycanopathies",
    "Glycogen Storage Disease Type IIb"}
n_p = collections.Counter(r["indication_mesh_id"] for r in frozen)
n_4 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "4")
n_5 = collections.Counter(r["indication_mesh_id"] for r in frozen if r["ccatnum"] == "5")
terms = []
for i in setids["PRIMARY"]:
    t = ind[i]["indication_mesh_term"]
    rule = ("R1 Minikel top-level branch C14" if topl[i] == "C14" else
            "R3 tagged cardiovascular, filed under C23 signs/symptoms" if (
                tagged(i) and topl[i] == "C23") else
            "R4 tagged cardiovascular, filed under %s" % topl[i] if tagged(i) else
            "R2 lipid/atherosclerosis token")
    terms.append(dict(mesh_id=i, mesh_term=t, minikel_areas=ind[i]["areas"],
                      minikel_top_level=topl[i], rule=rule, indication_area=AREA[t],
                      pairs=n_p[i], phase3_pairs=n_4[i], launched_pairs=n_5[i],
                      in_sensitivity_1=i in setids["SENS1"], in_sensitivity_2=i in setids["SENS2"],
                      flag=("not primarily cardiac" if t in NOT_PRIMARILY_CARDIAC else "")))
terms.sort(key=lambda r: (AREA_ORDER.index(r["indication_area"]), -r["pairs"]))
with open(os.path.join(OUT, "cv_terms_primary.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(terms[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(terms)
fl = [r for r in terms if r["flag"]]
K["sets"]["PRIMARY"]["flagged_not_primarily_cardiac"] = dict(
    terms=len(fl), pairs=sum(r["pairs"] for r in fl),
    targets_only_via_flagged=len(tset["PRIMARY"] - {
        r["gene"] for r in pairs["PRIMARY"]
        if r["indication_mesh_term"] not in NOT_PRIMARILY_CARDIAC}))
print("flagged 'not primarily cardiac' (kept): %d terms, %d pairs ; targets in the set ONLY "
      "through a flagged term: %d" % (len(fl), sum(r["pairs"] for r in fl),
                                     K["sets"]["PRIMARY"]["flagged_not_primarily_cardiac"][
                                         "targets_only_via_flagged"]))
print("rule counts: %s" % dict(collections.Counter(r["rule"][:2] for r in terms)))

# ============================================================== K02 gene records
hdr("K02", "GENE RECORDS (shipped mapping + HGNC Ensembl-ID rule)")
coords = {r["symbol"]: dict(ensembl_id=r["ensembl_id"], chrom=r["chrom"], start=int(r["start"]),
                            end=int(r["end"]), note="shipped mapping")
          for r in load(D1B + r"\gene_coords_grch37.tsv")}
alltargets = sorted({r["gene"] for r in frozen})
unm = [s for s in alltargets if s not in coords]
gtf, g_ = {}, lambda k, s: re.search(k + ' "([^"]+)"', s).group(1)
with gzip.open(RAW + r"\Homo_sapiens.GRCh37.87.gtf.gz", "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t", 9)
        if f[2] == "gene":
            gtf[g_("gene_id", f[8])] = (g_("gene_name", f[8]), f[0], int(f[3]), int(f[4]))


def hgnc(sym):
    f = os.path.join(HCACHE, sym + ".json")
    if not os.path.exists(f):
        req = urllib.request.Request("https://rest.genenames.org/fetch/symbol/" + sym,
                                     headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            open(f, "wb").write(r.read())
    docs = json.load(open(f, encoding="utf-8"))["response"]["docs"]
    return docs[0] if docs else None


unresolved = {}
for s in unm:
    h = hgnc(s)
    eid = (h or {}).get("ensembl_gene_id", "")
    if h is None:
        unresolved[s] = "symbol not an approved HGNC symbol"
    elif not eid:
        unresolved[s] = "HGNC lists no Ensembl gene ID"
    elif eid not in gtf or gtf[eid][1] not in CHROMS:
        unresolved[s] = "HGNC Ensembl ID %s is not a primary-chromosome gene in GRCh37.87" % eid
    else:
        n_, c_, a_, b_ = gtf[eid]
        coords[s] = dict(ensembl_id=eid, chrom=c_, start=a_, end=b_,
                         note="HGNC rule: GRCh37.87 record named %s" % n_)
K["gene_records"] = dict(unmapped_as_shipped=len(unm), resolved=len(unm) - len(unresolved),
                         unresolved=len(unresolved),
                         unresolved_in_primary=sorted(s for s in unresolved
                                                      if s in tset["PRIMARY"]))
print("unmapped as shipped (all 2,517 targets): %d ; resolved by HGNC rule: %d ; unresolved: %d"
      % (len(unm), len(unm) - len(unresolved), len(unresolved)))
for k in SETDEF:
    u = sorted(s for s in unresolved if s in tset[k])
    print("  %-8s unmapped as shipped %2d -> unresolved %d %s"
          % (k, sum(1 for s in unm if s in tset[k]), len(u), u))
with open(os.path.join(OUT, "gene_record_resolution.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["gene", "in_primary_set", "result", "ensembl_id", "detail"])
    for s in unm:
        w.writerow([s, s in tset["PRIMARY"], "UNRESOLVED" if s in unresolved else "RESOLVED",
                    coords.get(s, {}).get("ensembl_id", ""),
                    unresolved.get(s) or coords[s]["note"]])

# ============================================================== K03 eQTLGen lookup
hdr("K03", "eQTLGen LOOKUP (one pass over the FDR<0.05 file; cached)")
tested = {r["Gene"] for r in load(D1B + r"\eqtlgen_tested_universe.tsv")}
e2s = collections.defaultdict(list)
for s, c in coords.items():
    if s in set(alltargets):
        e2s[c["ensembl_id"]].append(s)
cache = os.path.join(OUT, "eqtlgen_target_hits.tsv")
hit = {}
if os.path.exists(cache):
    for r in load(cache):
        hit[r["gene"]] = dict(any_fdr=r["any_fdr05_row"] == "1", n_sig=int(r["n_p5e8_rows"]),
                              min_dist=(int(r["min_dist_bp"]) if r["min_dist_bp"] else None),
                              n100=int(r["n_within_100kb"]), best_p_100=r["best_p_within_100kb"],
                              best_snp_100=r["best_snp_within_100kb"], best_p_any=r["best_p_any"])
    print("read cache %s (%d genes); delete it to rescan" % (cache, len(hit)))
else:
    for s in e2s.values():
        for x in s:
            hit[x] = dict(any_fdr=False, n_sig=0, min_dist=None, n100=0, best_p_100="",
                          best_snp_100="", best_p_any="")
    SIG = (RAW + r"\2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-"
           r"BonferroniAdded.txt.gz")
    with gzip.open(SIG, "rt", encoding="utf-8") as fh:
        hd = fh.readline().rstrip("\n").split("\t")
        ix = {c: hd.index(c) for c in ("Pvalue", "SNP", "SNPChr", "SNPPos", "Gene")}
        for line in fh:
            f = line.split("\t", 9)
            syms = e2s.get(f[ix["Gene"]])
            if not syms:
                continue
            p = float(f[ix["Pvalue"]])
            for s in syms:
                h, c = hit[s], coords[s]
                h["any_fdr"] = True
                if p >= P_INSTR or f[ix["SNPChr"]] != c["chrom"]:
                    continue
                pos = int(f[ix["SNPPos"]])
                d = max(c["start"] - pos, pos - c["end"], 0)
                h["n_sig"] += 1
                if h["min_dist"] is None or d < h["min_dist"]:
                    h["min_dist"] = d
                if h["best_p_any"] == "" or p < float(h["best_p_any"]):
                    h["best_p_any"] = f[ix["Pvalue"]]
                if d <= 100_000:
                    h["n100"] += 1
                    if h["best_p_100"] == "" or p < float(h["best_p_100"]):
                        h["best_p_100"], h["best_snp_100"] = f[ix["Pvalue"]], f[ix["SNP"]]
    with open(cache, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "ensembl_id", "any_fdr05_row", "n_p5e8_rows", "min_dist_bp",
                    "n_within_100kb", "best_p_within_100kb", "best_snp_within_100kb",
                    "best_p_any"])
        for s in sorted(hit):
            h = hit[s]
            w.writerow([s, coords[s]["ensembl_id"], int(h["any_fdr"]), h["n_sig"],
                        "" if h["min_dist"] is None else h["min_dist"], h["n100"],
                        h["best_p_100"], h["best_snp_100"], h["best_p_any"]])
    print("scanned eQTLGen FDR<0.05 file; wrote %s" % cache)


def status(s, w=100_000):
    if s not in coords:
        return "UNRESOLVED"
    if coords[s]["ensembl_id"] not in tested:
        return "NO_EQTL_DATA"
    d = hit[s]["min_dist"]
    return "INSTRUMENTED" if (d is not None and d <= w) else "NO_QUALIFYING_INSTRUMENT"


def subtype(s):
    if status(s) != "NO_QUALIFYING_INSTRUMENT":
        return ""
    if hit[s]["n_sig"]:
        return "T1 significant cis-eQTL only outside +-100 kb"
    return ("T2 cis-eQTL at FDR<0.05, none at p<5e-8" if hit[s]["any_fdr"]
            else "T3 no cis-eQTL at FDR<0.05")


# reproduce the shipped states for every target that was already mapped
ann = {r["symbol"]: r for r in load(ROOT + r"\deliverables\d1e\classA_annotation.tsv")}
OLD = {"INSTRUMENTED": "INSTRUMENTED", "ABSENT": "NO_EQTL_DATA",
       "TESTED_NULL": "NO_QUALIFYING_INSTRUMENT"}
bad = [s for s, r in ann.items() if r["state"] != "UNMAPPED" and OLD[r["state"]] != status(s)]
print("targets mapped as shipped: %d ; status disagreeing with shipped state: %d %s"
      % (sum(1 for r in ann.values() if r["state"] != "UNMAPPED"), len(bad), bad[:5]))
assert not bad
K["eqtlgen_tested_genes"] = len(tested)

# ============================================================== annotations
hpa = {}
with zipfile.ZipFile(RAW + r"\proteinatlas.tsv.zip").open("proteinatlas.tsv") as fh:
    for r in csv.DictReader((l.decode("utf-8") for l in fh), delimiter="\t"):
        if r.get("Ensembl", "").strip():
            hpa[r["Ensembl"].strip()] = r
tpm = {}
with gzip.open(RAW + r"\gtex_v10_median_tpm.gct.gz", "rt", encoding="utf-8") as fh:
    fh.readline(); fh.readline()
    hd = fh.readline().rstrip("\n").split("\t")
    want = {t: hd.index(t) for t in ("Whole_Blood", "Heart_Left_Ventricle",
                                     "Heart_Atrial_Appendage", "Liver", "Artery_Aorta",
                                     "Artery_Coronary")}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        try:
            v = [float(x) for x in f[2:]]
        except ValueError:
            continue
        e = f[0].split(".")[0]
        tpm[e] = {k: float(f[i]) for k, i in want.items()}
        tpm[e]["max_tissue"] = hd[2 + v.index(max(v))]
        tpm[e]["max_tpm"] = max(v)
card = {}
for qtd, lab in CARDIAC:
    s_ = set()
    with gzip.open(RAW + "\\%s.permuted.tsv.gz" % qtd, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                if float(r["pvalue"]) < P_INSTR:
                    s_.add(r["molecular_trait_id"])
            except ValueError:
                pass
    card[lab] = s_
chembl = {r["gene"]: r for r in csv.DictReader(
    open(os.path.join(MP, "drug_links", "chembl_target_summary.csv"), encoding="utf-8-sig"))}


def annot(s):
    e = coords.get(s, {}).get("ensembl_id", "")
    h, t = hpa.get(e), tpm.get(e, {})
    spec = (h.get("RNA tissue specificity") or "").strip() if h else ""
    secl = (h.get("Secretome location") or "").strip() if h else ""
    pcls = (h.get("Protein class") or "") if h else ""
    return dict(
        hpa_rna_tissue_specificity=spec,
        tissue_restricted=("YES" if spec in RESTRICTED else "NO" if spec in NOT_RESTRICTED
                           else "NA"),
        secreted=(("YES" if (secl or "secreted" in pcls.lower()) else "NO") if h else "NA"),
        gtex_blood_tpm=t.get("Whole_Blood", ""), gtex_heart_LV_tpm=t.get("Heart_Left_Ventricle", ""),
        gtex_heart_atrial_tpm=t.get("Heart_Atrial_Appendage", ""), gtex_liver_tpm=t.get("Liver", ""),
        gtex_artery_aorta_tpm=t.get("Artery_Aorta", ""),
        gtex_max_tissue=t.get("max_tissue", ""), gtex_max_tpm=t.get("max_tpm", ""),
        **{"instrument_%s_upper_bound" % lab: (e in card[lab]) for _, lab in CARDIAC},
        instrument_cardiac_union_upper_bound=any(e in card[lab] for _, lab in CARDIAC))


# the annotation logic must reproduce the shipped file for already-mapped targets
ab = sum(1 for s, r in ann.items() if r["state"] != "UNMAPPED"
         and (annot(s)["tissue_restricted"], annot(s)["secreted"]) !=
         (r["tissue_restricted"], r["secreted"]))
print("HPA annotation disagreements with shipped classA_annotation.tsv: %d" % ab)
assert ab == 0

# ============================================================== per-target tables
rows_all = []
for s in alltargets:
    c = coords.get(s, {})
    h = hit.get(s, {})
    r = dict(gene=s, ensembl_id=c.get("ensembl_id", ""),
             gene_record=("unresolved: " + unresolved[s]) if s in unresolved else c["note"],
             chromosome=c.get("chrom", ""),
             status=status(s), no_instrument_subtype=subtype(s),
             no_eqtl_data_subtype=("" if status(s) != "NO_EQTL_DATA" else
                                   "X-linked: chromosome not analysed by eQTLGen"
                                   if c["chrom"] == "X" else
                                   "mitochondrial: genome not analysed by eQTLGen"
                                   if c["chrom"] == "MT" else
                                   "Y-linked: chromosome not analysed by eQTLGen"
                                   if c["chrom"] == "Y" else
                                   "autosomal gene not tested"),
             n_cis_eqtl_p5e8_within_100kb=h.get("n100", ""),
             best_p_within_100kb=h.get("best_p_100", ""),
             best_snp_within_100kb=h.get("best_snp_100", ""),
             n_cis_eqtl_p5e8_any_distance=h.get("n_sig", ""),
             nearest_p5e8_variant_bp_from_gene_body=("" if h.get("min_dist") is None
                                                     else h["min_dist"]),
             best_p_any_distance=h.get("best_p_any", ""),
             status_window_500kb=status(s, 500_000), status_window_1Mb=status(s, 1_000_000))
    for k in SETDEF:
        r["in_" + k] = s in tset[k]
        r["phase_" + k] = PH[phase[k].get(s, 0)]
    r.update(annot(s))
    rows_all.append(r)
byg = {r["gene"]: r for r in rows_all}

prim = []
for r in rows_all:
    if not r["in_PRIMARY"]:
        continue
    s = r["gene"]
    mine = [p for p in pairs["PRIMARY"] if p["gene"] == s]
    ar = collections.defaultdict(int)
    for p in mine:
        a = AREA[p["indication_mesh_term"]]
        ar[a] = max(ar[a], int(p["ccatnum"]))
    top = max(int(p["ccatnum"]) for p in mine)
    ch = chembl.get(s, {})
    q = dict(r)
    q.update(n_cv_pairs=len(mine),
             cv_indications_at_highest_phase="; ".join(sorted(
                 p["indication_mesh_term"] for p in mine if int(p["ccatnum"]) == top)),
             indication_areas="; ".join(a for a in AREA_ORDER if a in ar),
             **{"area_phase: " + a: PH[ar.get(a, 0)] for a in AREA_ORDER},
             chembl_checked=bool(ch),
             chembl_verified_strict=(ch.get("verified_strict") == "True" if ch else ""),
             chembl_verified_permissive=(ch.get("verified_permissive") == "True" if ch else ""),
             chembl_strict_drugs=ch.get("strict_drugs", ""),
             chembl_strict_fail_reason=ch.get("strict_fail_reason", ""),
             example_eligible=bool(ch) and ch.get("verified_strict") == "True"
             and top >= 4 and phase["SENS2"].get(s, 0) == top,
             chembl_result=ch.get("chembl_supports_minikel_label", ""),
             chembl_approved_cv_drugs=ch.get("approved_drugs_with_approved_cv_indication", ""),
             chembl_phase3_cv_drugs=ch.get("drugs_with_phase3_cv_indication_only", ""))
    prim.append(q)
for name, data in (("per_target_all_2517.csv", rows_all), ("per_target_primary.csv", prim)):
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(data)
    print("wrote %s (%d rows, %d columns)" % (name, len(data), len(data[0])))

# ============================================================== K04 bounds
hdr("K04", "UNINSTRUMENTABLE TARGETS: LOWER AND UPPER BOUND, BY SET")
ST = ["INSTRUMENTED", "NO_EQTL_DATA", "NO_QUALIFYING_INSTRUMENT", "UNRESOLVED"]


NONAUTO = ("X", "Y", "MT")


def nonauto(s):
    return coords.get(s, {}).get("chrom") in NONAUTO


def tally(genes, w=100_000):
    c = collections.Counter(status(s, w) for s in genes)
    n = len(genes)
    lo, up = c["NO_EQTL_DATA"], c["NO_EQTL_DATA"] + c["NO_QUALIFYING_INSTRUMENT"]
    nx = sum(1 for s in genes if nonauto(s))
    xl = sum(1 for s in genes if nonauto(s) and status(s, w) == "NO_EQTL_DATA")
    bych = collections.Counter(coords[s]["chrom"] for s in genes
                               if nonauto(s) and status(s, w) == "NO_EQTL_DATA")
    # no_eqtl_data_x_linked is kept as the name of the combined X + Y + MT count (the name
    # predates the 2026-10-11 relabel); the per-chromosome counts sit beside it.
    return dict(n=n, instrumented=c["INSTRUMENTED"], no_eqtl_data=lo,
                no_eqtl_data_x_linked=xl, no_eqtl_data_autosomal=lo - xl,
                no_eqtl_data_chrX=bych["X"], no_eqtl_data_chrMT=bych["MT"],
                no_eqtl_data_chrY=bych["Y"],
                n_not_on_x_y_mt=n - nx, autosomal_lower_bound=lo - xl,
                autosomal_lower_pct=pct(lo - xl, n - nx),
                no_qualifying_instrument=c["NO_QUALIFYING_INSTRUMENT"],
                unresolved=c["UNRESOLVED"], lower_bound=lo, lower_pct=pct(lo, n),
                upper_bound=up, upper_pct=pct(up, n), instrumented_pct=pct(c["INSTRUMENTED"], n),
                no_qualifying_instrument_pct=pct(c["NO_QUALIFYING_INSTRUMENT"], n))


def show(label, t):
    print("%-34s N=%4d | instrumented %4d (%5.1f%%) | no eQTL data %4d (%5.1f%%) [LOWER] | "
          "no qualifying instrument %3d (%4.1f%%) | upper bound %4d (%5.1f%%) | unresolved %d"
          " || X-linked %d + mitochondrial %d + Y %d ; AUTOSOMAL-ONLY lower %d / %d (%s%%)"
          % (label, t["n"], t["instrumented"], t["instrumented_pct"], t["no_eqtl_data"],
             t["lower_pct"], t["no_qualifying_instrument"], t["no_qualifying_instrument_pct"],
             t["upper_bound"], t["upper_pct"], t["unresolved"], t["no_eqtl_data_chrX"],
             t["no_eqtl_data_chrMT"], t["no_eqtl_data_chrY"],
             t["autosomal_lower_bound"], t["n_not_on_x_y_mt"], t["autosomal_lower_pct"]))


K["overall"] = {}
for k in SETDEF:
    K["overall"][k] = tally(tset[k])
    show(k, K["overall"][k])
K["overall"]["ALL_2517"] = tally(alltargets)
show("ALL 2,517 Minikel targets", K["overall"]["ALL_2517"])
sub = collections.Counter(subtype(s) for s in tset["PRIMARY"] if subtype(s))
K["no_instrument_subtypes_primary"] = dict(sub)
print("PRIMARY, sub-types of 'no qualifying instrument': %s" % dict(sub))
K["no_instrument_subtypes_all"] = dict(collections.Counter(
    subtype(s) for s in alltargets if subtype(s)))

# ============================================================== K05 by phase
hdr("K05", "BY HIGHEST PHASE (phase taken over the same set that defines membership)")
K["by_phase"] = {}
for k in SETDEF:
    K["by_phase"][k] = {}
    print("-- %s" % k)
    for ph in (1, 2, 3, 4, 5):
        g = [s for s in tset[k] if phase[k][s] == ph]
        K["by_phase"][k][PH[ph]] = tally(g)
        show("   " + PH[ph], K["by_phase"][k][PH[ph]])
t5 = [[K["by_phase"]["PRIMARY"][PH[p]]["no_eqtl_data"],
       K["by_phase"]["PRIMARY"][PH[p]]["n"] - K["by_phase"]["PRIMARY"][PH[p]]["no_eqtl_data"]]
      for p in (1, 2, 3, 4, 5)]
chi = stats.chi2_contingency(t5)
K["by_phase_chi2_primary_no_eqtl_data"] = dict(chi2=chi[0], dof=chi[2], p=chi[1],
                                               note="descriptive, not pre-registered")
print("descriptive chi-square, no-eQTL-data vs rest across 5 phases (PRIMARY): "
      "chi2=%.2f dof=%d p=%.3f" % (chi[0], chi[2], chi[1]))

# ============================================================== K06 launched + ChEMBL
hdr("K06", "LAUNCHED TARGETS: STRICT, APPROVED (FDA/EMA RECORD) AND PERMISSIVE VERIFICATION")


# "approved" (added 2026-10-11): the Methods definition. Strict-verified in ChEMBL AND at least
# one strict-rule drug has an approval record at the FDA (Drugs@FDA) or the EMA (medicines
# table). The record check is manuscript/revision_sources/check_approved_35.py; it covers the
# strict-verified launched and Phase III targets of the primary set, which contain those of
# every other set.
_AP = os.path.join(MP, "manuscript", "revision_sources", "approved_35_by_target.csv")
APPROVED = {r["target"] for r in csv.DictReader(open(_AP, encoding="utf-8"))
            if r["approved_fda_or_ema"] == "True"}
APPROVED_CHECKED = {r["target"] for r in csv.DictReader(open(_AP, encoding="utf-8"))}


def ver(s, kind):
    if kind == "approved":
        if chembl.get(s, {}).get("verified_strict") != "True":
            return False
        assert s in APPROVED_CHECKED, "strict-verified target not in the FDA/EMA check: " + s
        return s in APPROVED
    return chembl.get(s, {}).get("verified_" + kind) == "True"


K["launched"] = {}
for k in SETDEF:
    L = sorted(s for s in tset[k] if phase[k][s] == 5)
    K["launched"][k] = dict(all_launched=tally(L))
    print("-- %s" % k)
    show("   Minikel launched", K["launched"][k]["all_launched"])
    for kind in ("strict", "approved", "permissive"):
        V = [s for s in L if ver(s, kind)]
        K["launched"][k][kind] = dict(
            tally=tally(V),
            no_eqtl_data=sorted(s for s in V if status(s) == "NO_EQTL_DATA"),
            no_eqtl_data_x_linked=sorted(s for s in V if status(s) == "NO_EQTL_DATA"
                                         and nonauto(s)),
            no_qualifying_instrument=sorted(s for s in V
                                            if status(s) == "NO_QUALIFYING_INSTRUMENT"),
            beyond_PCSK9=dict(
                lower=sum(1 for s in V if status(s) == "NO_EQTL_DATA" and s != "PCSK9"),
                upper=sum(1 for s in V if status(s) in ("NO_EQTL_DATA",
                                                        "NO_QUALIFYING_INSTRUMENT")
                          and s != "PCSK9")))
        show("   ...verified, %-10s" % kind.upper(), K["launched"][k][kind]["tally"])
        print("      no eQTL data: %s" % ", ".join(K["launched"][k][kind]["no_eqtl_data"]))
        print("      no qualifying instrument: %s ; beyond PCSK9 lower %d upper %d"
              % (", ".join(K["launched"][k][kind]["no_qualifying_instrument"]),
                 K["launched"][k][kind]["beyond_PCSK9"]["lower"],
                 K["launched"][k][kind]["beyond_PCSK9"]["upper"]))
    po = sorted(s for s in L if ver(s, "permissive") and not ver(s, "strict"))
    K["launched"][k]["permissive_only"] = {s: (status(s), chembl[s]["strict_fail_reason"])
                                           for s in po}
    K["launched"][k]["not_verified_either"] = {
        s: (status(s), chembl.get(s, {}).get("chembl_supports_minikel_label", ""))
        for s in L if not ver(s, "permissive")}
    print("      permissive only (in MANUAL_REVIEW.csv): %s" % ", ".join(
        "%s[%s]" % (s, status(s)[:6]) for s in po))
    print("      not verified by either rule: %s" % ", ".join(
        sorted(K["launched"][k]["not_verified_either"])))
P3 = sorted(s for s in tset["PRIMARY"] if phase["PRIMARY"][s] == 4)
K["phase3_with_chembl_approved_cv_drug_primary"] = dict(
    strict={s: status(s) for s in P3 if ver(s, "strict")},
    approved={s: status(s) for s in P3 if ver(s, "approved")},
    permissive={s: status(s) for s in P3 if ver(s, "permissive")})
print("PRIMARY, Minikel Phase III with an approved CV drug in ChEMBL (labels NOT edited):")
print("   strict: %s" % K["phase3_with_chembl_approved_cv_drug_primary"]["strict"])
print("   permissive: %s" % K["phase3_with_chembl_approved_cv_drug_primary"]["permissive"])
L = [s for s in tset["PRIMARY"] if phase["PRIMARY"][s] == 5]
K["chembl_verified_any_label_primary"] = {
    kind: tally([s for s in L + P3 if ver(s, kind)])
    for kind in ("strict", "approved", "permissive")}
K["strict_not_approved_primary"] = sorted(s for s in L + P3 if ver(s, "strict")
                                          and not ver(s, "approved"))
print("PRIMARY strict-verified without an FDA or EMA approval record: %s"
      % K["strict_not_approved_primary"])
for kind in ("strict", "approved", "permissive"):
    show("PRIMARY Phase III or launched, %s" % kind, K["chembl_verified_any_label_primary"][kind])
EL = sorted(q["gene"] for q in prim if q["example_eligible"])
K["example_eligible_primary"] = {
    "launched_no_eqtl_data": [g for g in EL if phase["PRIMARY"][g] == 5
                              and status(g) == "NO_EQTL_DATA"],
    "launched_no_qualifying_instrument": [g for g in EL if phase["PRIMARY"][g] == 5
                                          and status(g) == "NO_QUALIFYING_INSTRUMENT"],
    "launched_instrumented": [g for g in EL if phase["PRIMARY"][g] == 5
                              and status(g) == "INSTRUMENTED"],
    "phase3_uninstrumentable": [g for g in EL if phase["PRIMARY"][g] == 4
                                and status(g) != "INSTRUMENTED"]}
print("EXAMPLE-ELIGIBLE (strict-verified, same phase in SENS2): %s"
      % K["example_eligible_primary"])

mm = []
for s in sorted(tset["PRIMARY"]):
    ph = phase["PRIMARY"][s]
    if ph < 4:
        continue
    res = chembl[s]["chembl_supports_minikel_label"]
    kind = ("Launched in Minikel; not strict-verified in ChEMBL" if ph == 5
            and not ver(s, "strict") else
            "Phase III in Minikel; ChEMBL lists an approved drug with an approved CV indication"
            if ph == 4 and ver(s, "permissive") else "")
    if kind:
        mm.append(dict(gene=s, minikel_phase_primary_set=PH[ph], mismatch=kind,
                       verified_strict=ver(s, "strict"),
                       verified_permissive=ver(s, "permissive"),
                       strict_fail_reason=chembl[s]["strict_fail_reason"],
                       chembl_result=res, status=status(s),
                       minikel_cv_indications_at_that_phase=byg[s] and next(
                           q["cv_indications_at_highest_phase"] for q in prim
                           if q["gene"] == s),
                       chembl_approved_cv_drugs=chembl[s][
                           "approved_drugs_with_approved_cv_indication"],
                       chembl_phase3_cv_drugs=chembl[s]["drugs_with_phase3_cv_indication_only"],
                       n_chembl_drugs_with_mechanism=chembl[s]["n_drugs_with_mechanism"]))
with open(os.path.join(OUT, "label_chembl_mismatches.csv"), "w", newline="",
          encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(mm[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(mm)
K["label_mismatches_primary"] = dict(collections.Counter(m["mismatch"] for m in mm))
print("label / ChEMBL mismatches (PRIMARY): %s" % K["label_mismatches_primary"])

# ============================================================== K07 window sensitivity
hdr("K07", "CIS-WINDOW SENSITIVITY (+-100 kb primary; +-500 kb; +-1 Mb)")
K["window"] = {}
for k in ("PRIMARY", "SENS1", "SENS2"):
    K["window"][k] = {}
    L = [s for s in tset[k] if phase[k][s] == 5]
    for wn, w in WINDOWS.items():
        K["window"][k][wn] = dict(overall=tally(tset[k], w), launched=tally(L, w))
        show("%s +-%s all" % (k, wn), K["window"][k][wn]["overall"])
        show("%s +-%s launched" % (k, wn), K["window"][k][wn]["launched"])
moved = sorted((s, hit[s]["min_dist"]) for s in tset["PRIMARY"]
               if status(s) == "NO_QUALIFYING_INSTRUMENT" and hit[s]["n_sig"])
K["window"]["primary_targets_gaining_instrument"] = [
    dict(gene=s, nearest_bp=d, at=("500kb" if d <= 500_000 else "1Mb"),
         phase=PH[phase["PRIMARY"][s]]) for s, d in moved]
print("PRIMARY targets that gain an instrument in a wider window: %s"
      % K["window"]["primary_targets_gaining_instrument"])
dists = [hit[s]["min_dist"] for s in alltargets if s in hit and hit[s]["min_dist"] is not None]
K["window"]["max_distance_observed_bp"] = max(dists)
print("largest gene-body distance of any p<5e-8 variant eQTLGen reports for a target: %d bp "
      "(eQTLGen tests SNPs within 1 Mb of the gene CENTRE, so the +-1 Mb gene-body window is "
      "not fully covered)" % max(dists))

# ============================================================== K08 enrichment, expression
hdr("K08", "WHO IS UNINSTRUMENTABLE: TISSUE RESTRICTION, SECRETION, BLOOD EXPRESSION")
K["enrichment"] = {}
for k, genes in (("PRIMARY", sorted(tset["PRIMARY"])), ("SENS1", sorted(tset["SENS1"])),
                 ("SENS2", sorted(tset["SENS2"])), ("ALL_2517", alltargets)):
    K["enrichment"][k] = {}
    for fld in ("tissue_restricted", "secreted"):
        for grp in ("NO_EQTL_DATA", "NO_QUALIFYING_INSTRUMENT"):
            q = lambda st, v: sum(1 for s in genes if byg[s]["status"] == st
                                  and byg[s][fld] == v)
            r = or_ci(q(grp, "YES"), q(grp, "NO"), q("INSTRUMENTED", "YES"),
                      q("INSTRUMENTED", "NO"))
            K["enrichment"][k]["%s_vs_instrumented__%s" % (grp, fld)] = r
            print("%-8s %-25s vs INSTRUMENTED x %-17s OR %6.2f [%5.2f, %6.2f] p=%.3g  %s"
                  % (k, grp, fld, r["odds_ratio"], r["ci_lo"], r["ci_hi"], r["p_fisher"],
                     r["table"]))
K["blood_tpm"] = {}
for st in ST[:3]:
    v = [float(byg[s]["gtex_blood_tpm"]) for s in tset["PRIMARY"]
         if byg[s]["status"] == st and byg[s]["gtex_blood_tpm"] != ""]
    K["blood_tpm"][st] = dict(n=len(v), median=statistics.median(v),
                              pct_below_1=pct(sum(1 for x in v if x < 1), len(v)))
    print("PRIMARY GTEx whole-blood TPM, %-25s n=%3d median %.3f ; %%<1 TPM: %.1f"
          % (st, len(v), statistics.median(v), K["blood_tpm"][st]["pct_below_1"]))
K["tissue_yield_primary"] = {"whole_blood_eQTLGen": dict(
    n=K["overall"]["PRIMARY"]["instrumented"], pct=K["overall"]["PRIMARY"]["instrumented_pct"])}
N = len(tset["PRIMARY"])
for lab in [l for _, l in CARDIAC] + ["cardiac_union"]:
    n_ = sum(1 for s in tset["PRIMARY"] if byg[s]["instrument_%s_upper_bound" % lab])
    K["tissue_yield_primary"][lab] = dict(n=n_, pct=pct(n_, N))
    print("PRIMARY instrument yield (upper bound) %-14s %3d / %d = %.1f%%"
          % (lab, n_, N, pct(n_, N)))
ne = [s for s in tset["PRIMARY"] if byg[s]["status"] == "NO_EQTL_DATA"]
K["no_eqtl_data_with_tissue_instrument_primary"] = dict(
    n=len(ne), any_of_4_tissues=sum(1 for s in ne
                                    if byg[s]["instrument_cardiac_union_upper_bound"]))
print("PRIMARY: of %d no-eQTL-data targets, %d have a p<5e-8 top cis-eQTL in >=1 of the 4 "
      "GTEx tissues (upper bound)" % (len(ne), K["no_eqtl_data_with_tissue_instrument_primary"][
          "any_of_4_tissues"]))

# ============================================================== K09 infrastructure
hdr("K09", "OUTCOME-GWAS FILE ATTRITION (logic of src/a3_cascade.py, re-applied per set)")
harm = set()
for l in open(RAW + r"\harmonised_list.txt"):
    m_ = re.search(r"/(GCST\d+)/harmonised/", l)
    if m_:
        harm.add(m_.group(1))
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
        sz = f[isz]
        if not (("European" in sz) and not any(a in sz for a in ANC)):
            continue
        n_ = max([int(x.replace(",", "")) for x in re.findall(r"[\d,]{4,}", sz)] or [0])
        for t in f[it].split(","):
            t = t.strip().lower()
            if t and (t not in best or n_ > best[t][1]):
                best[t] = (f[ia], n_)
skipped = {r["gwas_acc"] for r in load(D1F + r"\partc_skipped_accessions.tsv")}
probed_ok = {r["gwas_acc"] for r in load(D1F + r"\partc_results.tsv")}
K["gwas_file_attrition"] = {}
for k in SETDEF:
    mp = [r for r in pairs[k] if r["indication_mesh_term"].strip().lower() in best]
    acc = [best[r["indication_mesh_term"].strip().lower()][0] for r in mp]
    lost = sum(1 for a in acc if a in skipped)
    never = sum(1 for a in acc if a not in skipped and a not in probed_ok)
    K["gwas_file_attrition"][k] = dict(
        pairs=len(pairs[k]), pairs_with_eur_harmonised_gwas=len(mp), lost_unindexed=lost,
        lost_pct=pct(lost, len(mp)), index_status_never_probed=never,
        distinct_indications_with_gwas=len({r["indication_mesh_id"] for r in mp}))
    print("%-8s pairs %5d ; indication exact-matches a EUR-only harmonised GWAS: %4d pairs "
          "(%d indications) ; file not range-queryable: %4d (%.1f%%) ; never probed: %d"
          % (k, len(pairs[k]), len(mp), K["gwas_file_attrition"][k][
              "distinct_indications_with_gwas"], lost, pct(lost, len(mp)), never))

# ============================================================== K10 sex chromosomes
hdr("K10", "WHY 'NO eQTL DATA': AUTOSOMAL vs SEX-CHROMOSOME / MITOCHONDRIAL GENES")
gchrom = {e: v[1] for e, v in gtf.items()}
K["x_linked_blood_expressed_primary"] = sorted(
    s for s in tset["PRIMARY"] if nonauto(s) and byg[s]["gtex_blood_tpm"] != ""
    and float(byg[s]["gtex_blood_tpm"]) >= 1)
tc = collections.Counter(("autosomal" if gchrom.get(e, "?").isdigit() else gchrom.get(e, "not in GRCh37.87"))
                         for e in tested)
K["eqtlgen_tested_genes_by_chromosome"] = dict(tc)
print("eQTLGen tested genes by chromosome class: %s" % dict(tc))
K["no_eqtl_data_by_chromosome"] = {}
for k, genes in list((k, sorted(tset[k])) for k in SETDEF) + [("ALL_2517", alltargets)]:
    ne_ = [s for s in genes if byg[s]["status"] == "NO_EQTL_DATA"]
    xs = [s for s in ne_ if byg[s]["chromosome"] in NONAUTO]
    allx = [s for s in genes if byg[s]["chromosome"] in NONAUTO]
    xb = [s for s in xs if byg[s]["gtex_blood_tpm"] != "" and float(byg[s]["gtex_blood_tpm"]) >= 1]
    au = [s for s in ne_ if s not in xs]
    aut = [s for s in au if byg[s]["gtex_blood_tpm"] != ""]
    K["no_eqtl_data_by_chromosome"][k] = dict(
        no_eqtl_data=len(ne_), on_X_Y_MT=len(xs), autosomal=len(au),
        on_X=sum(1 for s in xs if byg[s]["chromosome"] == "X"),
        on_MT=sum(1 for s in xs if byg[s]["chromosome"] == "MT"),
        on_Y=sum(1 for s in xs if byg[s]["chromosome"] == "Y"),
        X_with_blood_tpm_ge1=sum(1 for s in xb if byg[s]["chromosome"] == "X"),
        MT_with_blood_tpm_ge1=sum(1 for s in xb if byg[s]["chromosome"] == "MT"),
        targets_on_X_Y_MT_in_set=len(allx), X_Y_MT_with_blood_tpm_ge1=len(xb),
        autosomal_with_tpm=len(aut),
        autosomal_blood_tpm_below_1=sum(1 for s in aut if float(byg[s]["gtex_blood_tpm"]) < 1),
        n=len(genes), autosomal_lower_bound_pct=pct(len(au), len(genes) - len(allx)))
    _c = K["no_eqtl_data_by_chromosome"][k]
    print("%-8s no eQTL data %3d = autosomal %3d + X/Y/MT %3d [X %d, MT %d, Y %d] (all %d X/Y/MT "
          "targets in the set; blood TPM >= 1: %d of the X-linked, %d of the mitochondrial) ; "
          "autosomal with blood TPM < 1: %d of %d ; "
          "autosomal-only lower bound %d / %d = %.1f%%"
          % (k, len(ne_), len(au), len(xs), _c["on_X"], _c["on_MT"], _c["on_Y"], len(allx),
             _c["X_with_blood_tpm_ge1"], _c["MT_with_blood_tpm_ge1"],
             K["no_eqtl_data_by_chromosome"][k]["autosomal_blood_tpm_below_1"], len(aut),
             len(au), len(genes) - len(allx), pct(len(au), len(genes) - len(allx))))
ga = [s for s in tset["PRIMARY"] if byg[s]["chromosome"] not in NONAUTO]
q_ = lambda st, v: sum(1 for s in ga if byg[s]["status"] == st and byg[s]["tissue_restricted"] == v)
K["enrichment"]["PRIMARY_autosomal_only"] = {
    "NO_EQTL_DATA_vs_instrumented__tissue_restricted": or_ci(
        q_("NO_EQTL_DATA", "YES"), q_("NO_EQTL_DATA", "NO"), q_("INSTRUMENTED", "YES"),
        q_("INSTRUMENTED", "NO"))}
e_ = K["enrichment"]["PRIMARY_autosomal_only"]["NO_EQTL_DATA_vs_instrumented__tissue_restricted"]
print("PRIMARY autosomal only: no eQTL data vs instrumented x tissue-restricted OR %.2f "
      "[%.2f, %.2f] p=%.3g %s" % (e_["odds_ratio"], e_["ci_lo"], e_["ci_hi"], e_["p_fisher"],
                                  e_["table"]))

json.dump(K, open(os.path.join(OUT, "key_numbers.json"), "w"), indent=1, default=str)
print("\nwrote key_numbers.json\nDONE")
