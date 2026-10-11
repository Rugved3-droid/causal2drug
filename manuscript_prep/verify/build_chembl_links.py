#!/usr/bin/env python
"""Item 5 -- target -> drug -> indication link from ChEMBL, for manual review.

    python manuscript_prep/verify/build_chembl_links.py

READ-ONLY on the project. Network: ChEMBL REST API (mechanism, molecule,
drug_indication, target) and NLM MeSH RDF (tree numbers). Every response is cached
under verify/cache/ so a re-run is offline and identical.

TARGET SCOPE. Every gene whose highest Minikel phase is Phase III or Launched in ANY of
three candidate CV sets, so no re-run is needed whichever set is approved:
   S0   the original 630-target atlas set (MeSH C14 top-level + 8 lipid tokens)
   S0A  S0 + C23-filed indications that Minikel also tags 'cardiovascular'
   S0AB S0A + indications filed elsewhere that Minikel also tags 'cardiovascular'
        (S0AB is the PRIMARY set, approved 2026-10-10)
Phase labels are read from the frozen pair table; no instrument state is touched.

LINK.  gene symbol -> ChEMBL human targets containing that gene (single protein,
complex, family; type recorded) -> mechanism rows -> parent molecule -> drug_indication.

CARDIOVASCULAR FLAG for a ChEMBL indication (MeSH descriptor):
   ind_mesh_tree_C14     any MeSH tree number starts with C14   (NLM MeSH RDF)
   ind_lipid_token       heading matches one of the atlas's 8 lipid/athero tokens
   ind_in_S0 / S0A / S0AB  the descriptor is one of the Minikel indications in that set
   ind_is_cv = ind_mesh_tree_C14 OR ind_lipid_token
FLAG_approved_not_cv = indication is approved in ChEMBL (max_phase_for_ind == 4) and
ind_is_cv is False.

Outputs (manuscript_prep/drug_links/):
   chembl_target_drug_indication.csv   one row per target x drug x indication, phase >= 3
   chembl_target_summary.csv           one row per target -- start the review here
   chembl_targets_without_mechanism.csv
"""
import collections
import csv
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"D:\causal2drug"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache", "chembl")
MCACHE = os.path.join(HERE, "cache", "mesh")
OUT = os.path.join(os.path.dirname(HERE), "drug_links")
for d in (CACHE, MCACHE, OUT):
    os.makedirs(d, exist_ok=True)
API = "https://www.ebi.ac.uk/chembl/api/data/"
TOK = ("hyperlipid", "hypercholesterol", "dyslipid", "hyperlipoprotein",
       "hypertriglycerid", "lipoprotein", "atheroscler", "arterioscler")
PH = {0: "", 1: "Preclinical", 2: "Phase I", 3: "Phase II", 4: "Phase III", 5: "Launched"}


def load(p, d="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=d))


def fetch(url, cache_dir):
    f = os.path.join(cache_dir, hashlib.sha1(url.encode()).hexdigest() + ".json")
    if os.path.exists(f):
        return json.load(open(f, encoding="utf-8"))
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read()
            json.loads(raw)
            open(f, "wb").write(raw)
            time.sleep(0.05)
            return json.loads(raw)
        except Exception as e:                       # noqa: BLE001
            if "404" in str(e):
                open(f, "w").write("null")
                return None
            time.sleep(2 + 3 * attempt)
            err = e
    raise RuntimeError("failed: %s (%r)" % (url, err))


def chembl(resource, **params):
    """All pages of a ChEMBL list resource."""
    params = dict(params, limit=1000, format="json")
    url = API + resource + "?" + urllib.parse.urlencode(params)
    out = []
    while url:
        d = fetch(url, CACHE)
        key = [k for k in d if k != "page_meta"][0]
        out += d[key]
        nxt = d["page_meta"].get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
    return out


def chunks(xs, n):
    xs = list(xs)
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


# ---------------------------------------------------------------- target scope
ind = {r["indication_mesh_id"]: r for r in load(ROOT + r"\data\raw\minikel_indic.tsv")}
topl = {r["indication_mesh_id"]: r["topl"]
        for r in load(ROOT + r"\data\raw\minikel_indic_topl_match.tsv")}
frozen = load(ROOT + r"\data\derived\pp_dedup_frozen.tsv")


def tagged_cv(i):
    return "cardiovascular" in [a.strip() for a in ind[i]["areas"].split(",")]


def in_s0(i):
    return topl[i] == "C14" or any(k in ind[i]["indication_mesh_term"].lower() for k in TOK)


SETS = {"S0": lambda i: in_s0(i),
        "S0A": lambda i: in_s0(i) or (topl[i] == "C23" and tagged_cv(i)),
        "S0AB": lambda i: in_s0(i) or tagged_cv(i)}
setids = {k: {i for i in ind if fn(i)} for k, fn in SETS.items()}
mx = {k: collections.defaultdict(int) for k in SETS}
for r in frozen:
    for k in SETS:
        if r["indication_mesh_id"] in setids[k] and r["ccatnum"].isdigit():
            mx[k][r["gene"]] = max(mx[k][r["gene"]], int(r["ccatnum"]))
genes = sorted(g for g in mx["S0AB"] if mx["S0AB"][g] >= 4)
print("targets at Phase III or Launched: S0 %d, S0A %d, S0AB %d"
      % tuple(sum(1 for g in mx[k] if mx[k][g] >= 4) for k in ("S0", "S0A", "S0AB")))
status = fetch(API + "status.json", CACHE)
print("ChEMBL: %s released %s" % (status["chembl_db_version"], status["chembl_release_date"]))

# ---------------------------------------------------------------- gene -> targets
gene_targets = {}
for g in genes:
    ts = chembl("target", **{
        "target_components__target_component_synonyms__component_synonym__iexact": g,
        "organism": "Homo sapiens"})
    keep = []
    for t in ts:
        ok = any(s.get("syn_type") == "GENE_SYMBOL" and s.get("component_synonym", "").upper()
                 == g.upper() for c in t.get("target_components", [])
                 for s in c.get("target_component_synonyms", []))
        if ok:
            keep.append(t)
    gene_targets[g] = keep
tinfo = {t["target_chembl_id"]: t for ts in gene_targets.values() for t in ts}
print("genes %d ; with >=1 ChEMBL human target %d ; distinct ChEMBL targets %d"
      % (len(genes), sum(1 for g in genes if gene_targets[g]), len(tinfo)))

# ---------------------------------------------------------------- mRNA-directed drugs
# ChEMBL attaches antisense / RNAi mechanisms to an mRNA target that carries no gene
# component, so the gene-symbol lookup above cannot find them (Phase 1 missed mipomersen
# and volanesorsen for this reason). The map below was made BY HAND, on 2026-10-10, from
# the complete list of ChEMBL mechanisms with action_type ANTISENSE INHIBITOR or
# RNAI INHIBITOR (62 rows; verify/antisense_mechanisms.tsv), reading the gene from
# `mechanism_of_action`. It is applied to every in-scope gene, not only APOB and APOC3.
MRNA_TARGETS = {
    "APOB": ["CHEMBL2364185"], "APOC3": ["CHEMBL4495602", "CHEMBL4834388"],
    "LPA": ["CHEMBL4662965"], "PCSK9": ["CHEMBL4630662"], "ANGPTL3": ["CHEMBL4622825"],
    "TTR": ["CHEMBL4296221"], "SERPINC1": ["CHEMBL4634874"], "VEGFA": ["CHEMBL4630576"],
    "ICAM1": ["CHEMBL4622822"], "IRS1": ["CHEMBL4622821"], "ADRB2": ["CHEMBL4630578"],
    "C5": ["CHEMBL4627245"], "BCL2": ["CHEMBL4622828"], "TP53": ["CHEMBL4622827"],
    "STAT3": ["CHEMBL4622829"], "TGFB2": ["CHEMBL4622826"], "PRKCA": ["CHEMBL4622824"],
    "FLT1": ["CHEMBL4662970"], "TMPRSS6": ["CHEMBL4834387"], "SOD1": ["CHEMBL3833461"],
    "HTT": ["CHEMBL4296641"], "CLU": ["CHEMBL4630579"], "HSPB1": ["CHEMBL4622823"],
    "ALAS1": ["CHEMBL4296741"], "LDHA": ["CHEMBL4627242"], "HAO1": ["CHEMBL4524134"],
    "MIR122": ["CHEMBL4483924"], "MIR21": ["CHEMBL1615322"],
}
anti = []
for a in ("ANTISENSE INHIBITOR", "RNAI INHIBITOR"):
    anti += chembl("mechanism", action_type=a)
with open(os.path.join(HERE, "antisense_mechanisms.tsv"), "w", encoding="utf-8",
          newline="") as fh:
    cols_ = ["parent_molecule", "target_chembl_id", "action_type", "max_phase",
             "mechanism_of_action"]
    fh.write(chr(9).join(cols_) + chr(10))
    for m in anti:
        fh.write(chr(9).join(str(x) for x in (
            m.get("parent_molecule_chembl_id"), m["target_chembl_id"], m["action_type"],
            m["max_phase"], m["mechanism_of_action"])) + chr(10))
mrna_added = {}
for g_, tids in MRNA_TARGETS.items():
    if g_ not in genes:
        continue
    for t in chembl("target", target_chembl_id__in=",".join(tids)):
        t = dict(t, _mrna=True)
        gene_targets[g_].append(t)
        tinfo[t["target_chembl_id"]] = t
        mrna_added.setdefault(g_, []).append(t["target_chembl_id"])
print("mRNA-target additions (manual map) for in-scope genes: %s" % mrna_added)

# ---------------------------------------------------------------- mechanisms
mech = []
for ch in chunks(sorted(tinfo), 40):
    mech += chembl("mechanism", target_chembl_id__in=",".join(ch))
mech_by_t = collections.defaultdict(list)
for m in mech:
    mech_by_t[m["target_chembl_id"]].append(m)
parents = sorted({m.get("parent_molecule_chembl_id") or m["molecule_chembl_id"] for m in mech})
print("mechanism rows %d ; distinct parent molecules %d" % (len(mech), len(parents)))

# ---------------------------------------------------------------- molecules
mol = {}
for ch in chunks(parents, 40):
    for m in chembl("molecule", molecule_chembl_id__in=",".join(ch),
                    only="molecule_chembl_id,pref_name,max_phase,first_approval,"
                         "withdrawn_flag,molecule_type"):
        mol[m["molecule_chembl_id"]] = m

# ---------------------------------------------------------------- indications
indic = collections.defaultdict(list)
for ch in chunks(parents, 40):
    for d in chembl("drug_indication", parent_molecule_chembl_id__in=",".join(ch)):
        indic[d["parent_molecule_chembl_id"]].append(d)
# collapse duplicates (ChEMBL stores one row per EFO mapping) to one row per MeSH id
for p in indic:
    best = {}
    for d in indic[p]:
        k = d.get("mesh_id") or ("EFO:" + str(d.get("efo_id")))
        ph = float(d.get("max_phase_for_ind") or 0)
        if k not in best or ph > best[k][0]:
            best[k] = (ph, d)
    indic[p] = [dict(d, _ph=ph) for ph, d in best.values()]
print("parent molecules with >=1 indication row: %d" % len(indic))


# ---------------------------------------------------------------- MeSH tree numbers
def mesh_trees(mid):
    if not mid:
        return []
    d = fetch("https://id.nlm.nih.gov/mesh/%s.json" % mid, MCACHE)
    if not d:
        return []
    tn = d.get("treeNumber", [])
    if isinstance(tn, str):
        tn = [tn]
    return sorted(x.rsplit("/", 1)[-1] for x in tn)


need = sorted({d.get("mesh_id") for p in indic for d in indic[p]
               if d["_ph"] >= 3 and d.get("mesh_id")})
print("distinct MeSH descriptors among phase>=3 indications: %d" % len(need))
trees = {m: mesh_trees(m) for m in need}
print("descriptors with no tree number returned: %d" % sum(1 for m in need if not trees[m]))

# ---------------------------------------------------------------- strict rule
# STRICT verification (added 2026-10-10 on the owner's instruction). A target is
# strict-verified when ONE drug satisfies all of:
#   S1  the drug is approved (ChEMBL max_phase = 4);
#   S2  the drug is linked to THIS gene directly: ChEMBL target type SINGLE PROTEIN or
#       PROTEIN COMPLEX containing the gene, or an mRNA target in the manual map.
#       Links through PROTEIN FAMILY, PROTEIN COMPLEX GROUP, SELECTIVITY GROUP,
#       CHIMERIC PROTEIN or PROTEIN-PROTEIN INTERACTION targets do not count;
#   S3  that drug has an approved indication (max_phase_for_ind = 4) that is strictly
#       cardiovascular: the MeSH descriptor is under C14 or matches a lipid token, AND
#       every disease branch it sits in is one of C14, C23 (signs/symptoms), C10
#       (nervous system: stroke), C08 (respiratory: pulmonary hypertension/embolism),
#       or, for lipid terms only, C18 and C16. A descriptor that also sits in any other
#       branch (C06, C11, C12, C16, C17, C19, C20, C04 ...) is not strictly cardiovascular.
# PERMISSIVE (the earlier rule) = S1 + any link + approved indication under C14 or lipid.
DIRECT_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX"}
CV_BRANCHES = {"C14", "C23", "C10", "C08"}


def strict_cv(tr, head):
    lip = any(k in head.lower() for k in TOK)
    br = {x.split(".")[0] for x in tr if x[:1] == "C"}
    if lip:
        return br <= (CV_BRANCHES | {"C18", "C16"})
    return any(x.startswith("C14") for x in tr) and br <= CV_BRANCHES


mk_launched = collections.defaultdict(set)   # gene -> MeSH ids labelled Launched, primary set
for r in frozen:
    if r["indication_mesh_id"] in setids["S0AB"] and r["ccatnum"] == "5":
        mk_launched[r["gene"]].add(r["indication_mesh_id"])

# ---------------------------------------------------------------- assemble
rows, nomech = [], []
summ = []
review = []
for g in genes:
    base = dict(gene=g, minikel_phase_S0=PH[mx["S0"].get(g, 0)],
                minikel_phase_S0A=PH[mx["S0A"].get(g, 0)],
                minikel_phase_S0AB=PH[mx["S0AB"].get(g, 0)])
    drugs = {}
    for t in gene_targets[g]:
        for m in mech_by_t.get(t["target_chembl_id"], []):
            p = m.get("parent_molecule_chembl_id") or m["molecule_chembl_id"]
            drugs.setdefault(p, []).append((t, m))
    if not drugs:
        nomech.append(dict(base, n_chembl_targets=len(gene_targets[g]),
                           chembl_targets=";".join("%s(%s)" % (t["target_chembl_id"],
                                                                t["target_type"])
                                                   for t in gene_targets[g])))
    appr, appr_cv, p3_cv, appr_noncv_only = [], [], [], []
    appr_strict, ev = [], []
    for p, tm in drugs.items():
        mi = mol.get(p, {})
        mphase = float(mi.get("max_phase") or 0)
        name = mi.get("pref_name") or p
        t, m = sorted(tm, key=lambda x: (not x[0].get("_mrna"),
                                         x[0]["target_type"] != "SINGLE PROTEIN",
                                         x[0]["target_chembl_id"]))[0]
        inds = indic.get(p, [])
        cv_appr = cv_p3 = False
        direct = any(x[0].get("_mrna") or x[0]["target_type"] in DIRECT_TYPES for x in tm)
        link_types = ";".join(sorted({("NUCLEIC-ACID (mRNA map)" if x[0].get("_mrna")
                                       else x[0]["target_type"]) for x in tm}))
        d_strict = False
        for d in inds:
            if d["_ph"] < 3:
                continue
            mid = d.get("mesh_id") or ""
            head = d.get("mesh_heading") or ""
            tr = trees.get(mid, [])
            c14 = any(x.startswith("C14") for x in tr)
            lip = any(k in head.lower() for k in TOK)
            iscv = c14 or lip
            cv_appr |= iscv and d["_ph"] >= 4
            cv_p3 |= iscv
            scv = strict_cv(tr, head)
            row_perm = mphase >= 4 and d["_ph"] >= 4 and iscv
            row_strict = row_perm and direct and scv
            d_strict |= row_strict
            pair = mid in mk_launched[g]
            if row_perm:
                ev.append(dict(drug=name, ind=head, mid=mid, tree=";".join(tr), direct=direct,
                               link_types=link_types, scv=scv, strict=row_strict, pair=pair,
                               action=m.get("action_type") or "",
                               moa=m.get("mechanism_of_action") or "",
                               first_approval=mi.get("first_approval") or "",
                               withdrawn=mi.get("withdrawn_flag")))
            rows.append(dict(
                base, target_chembl_id=t["target_chembl_id"], target_type=t["target_type"],
                link_method=("mRNA target, manual gene map" if t.get("_mrna")
                             else "protein target, gene symbol"),
                target_name=t["pref_name"], n_chembl_targets_for_drug_gene=len(tm),
                drug_chembl_id=p, drug_name=name, drug_max_phase=mphase,
                drug_first_approval=mi.get("first_approval") or "",
                drug_withdrawn=mi.get("withdrawn_flag"),
                molecule_type=mi.get("molecule_type") or "",
                action_type=m.get("action_type") or "",
                mechanism_of_action=m.get("mechanism_of_action") or "",
                indication_mesh_id=mid, indication_mesh_heading=head,
                indication_efo_term=d.get("efo_term") or "",
                max_phase_for_indication=d["_ph"],
                indication_approved=d["_ph"] >= 4,
                mesh_tree_numbers=";".join(tr), ind_mesh_tree_C14=c14, ind_lipid_token=lip,
                ind_in_S0=mid in setids["S0"], ind_in_S0A=mid in setids["S0A"],
                ind_in_S0AB=mid in setids["S0AB"], ind_is_cv=iscv,
                link_types_for_this_gene=link_types, link_direct=direct,
                ind_strict_cv=scv, row_permissive=row_perm, row_strict=row_strict,
                same_indication_minikel_launched=pair,
                FLAG_approved_not_cv=(d["_ph"] >= 4 and not iscv)))
        if d_strict:
            appr_strict.append(name)
        if mphase >= 4:
            appr.append(name)
            (appr_cv if cv_appr else appr_noncv_only).append(name)
        if cv_p3 and not cv_appr:
            p3_cv.append(name)
    top = max(mx["S0AB"].get(g, 0), 0)
    summ.append(dict(
        base, n_chembl_targets=len(gene_targets[g]), n_drugs_with_mechanism=len(drugs),
        n_approved_drugs=len(appr), n_approved_drugs_with_approved_cv_indication=len(appr_cv),
        approved_drugs_with_approved_cv_indication="; ".join(sorted(appr_cv)[:15]),
        n_approved_drugs_no_approved_cv_indication=len(appr_noncv_only),
        approved_drugs_no_approved_cv_indication="; ".join(sorted(appr_noncv_only)[:15]),
        drugs_with_phase3_cv_indication_only="; ".join(sorted(p3_cv)[:15]),
        chembl_supports_minikel_label=(
            "no ChEMBL mechanism for this gene" if not drugs else
            "YES approved drug with approved CV indication" if appr_cv else
            ("PARTIAL phase 3 CV indication only" if p3_cv else
             "NO approved or phase 3 CV indication found")),
        review_priority=("HIGH" if top == 5 and not appr_cv else
                         "MEDIUM" if top == 4 and not (appr_cv or p3_cv) else "low"),
        verified_permissive=bool(appr_cv), verified_strict=bool(appr_strict),
        n_strict_drugs=len(appr_strict),
        strict_drugs="; ".join(sorted(appr_strict)[:15]),
        same_indication_as_minikel_launched=any(e["pair"] and e["strict"] for e in ev),
        strict_fail_reason=(
            "" if appr_strict or not appr_cv else
            "no direct link: every approved CV drug is linked through a family or group target"
            if not any(e["direct"] for e in ev) else
            "no strictly cardiovascular approved indication on a directly linked drug"
            if not any(e["direct"] and e["scv"] for e in ev) else "other")))
    if appr_cv and not appr_strict:
        for e in sorted(ev, key=lambda e: (not e["direct"], not e["scv"], e["drug"], e["ind"])):
            review.append(dict(
                gene=g, minikel_phase_primary=PH[mx["S0AB"].get(g, 0)],
                minikel_phase_sens2=PH[mx["S0A"].get(g, 0)],
                minikel_launched_cv_indications="; ".join(sorted(
                    ind[i]["indication_mesh_term"] for i in mk_launched[g])),
                drug=e["drug"], drug_first_approval=e["first_approval"],
                drug_withdrawn=e["withdrawn"], action_type=e["action"],
                mechanism_of_action=e["moa"], chembl_link_types=e["link_types"],
                link_is_direct=e["direct"], approved_indication=e["ind"],
                indication_mesh_id=e["mid"], mesh_tree_numbers=e["tree"],
                indication_strictly_cv=e["scv"],
                why_row_fails_strict="; ".join(
                    ([] if e["direct"] else ["link is family/group level, not this gene"])
                    + ([] if e["scv"] else ["indication also sits outside cardiovascular "
                                            "MeSH branches"])),
                same_indication_as_minikel_launched=e["pair"],
                REVIEWER_accept_as_approved_cv_drug_for_this_target="",
                REVIEWER_note=""))


def write(name, data, cols=None):
    if not data:
        open(os.path.join(OUT, name), "w").write("")
        return
    cols = cols or list(data[0].keys())
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(data)


rows.sort(key=lambda r: (r["gene"], -r["max_phase_for_indication"], r["drug_name"],
                         r["indication_mesh_heading"]))
write("chembl_target_drug_indication.csv", rows)
write("chembl_target_summary.csv", summ)
write("chembl_targets_without_mechanism.csv", nomech)
write("MANUAL_REVIEW.csv", review)
print()
print("MANUAL_REVIEW.csv: %d evidence rows for %d targets where permissive = verified and "
      "strict = not verified" % (len(review), len({r["gene"] for r in review})))
for k in ("S0", "S0A", "S0AB"):
    L_ = [s_ for s_ in summ if s_["minikel_phase_" + k] == "Launched"]
    print("  %s launched %d: permissive %d, strict %d ; permissive-only: %s"
          % (k, len(L_), sum(s_["verified_permissive"] for s_ in L_),
             sum(s_["verified_strict"] for s_ in L_),
             ", ".join(sorted(s_["gene"] for s_ in L_ if s_["verified_permissive"]
                              and not s_["verified_strict"]))))
P3_ = [s_ for s_ in summ if s_["minikel_phase_S0AB"] == "Phase III"]
print("  S0AB Phase III %d: permissive %d (%s) ; strict %d (%s)"
      % (len(P3_), sum(s_["verified_permissive"] for s_ in P3_),
         ", ".join(sorted(s_["gene"] for s_ in P3_ if s_["verified_permissive"])),
         sum(s_["verified_strict"] for s_ in P3_),
         ", ".join(sorted(s_["gene"] for s_ in P3_ if s_["verified_strict"]))))
print("  strict-fail reasons: %s" % dict(collections.Counter(
    s_["strict_fail_reason"] for s_ in summ if s_["strict_fail_reason"])))
json.dump(dict(chembl=status, n_genes=len(genes), n_rows=len(rows),
               built_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
          open(os.path.join(OUT, "chembl_build_info.json"), "w"), indent=1)

print("\nrows (target x drug x indication, phase>=3): %d" % len(rows))
print("   approved-indication rows: %d ; of which FLAG_approved_not_cv: %d"
      % (sum(1 for r in rows if r["indication_approved"]),
         sum(1 for r in rows if r["FLAG_approved_not_cv"])))
print("genes without any ChEMBL mechanism row: %d -> %s"
      % (len(nomech), ", ".join(r["gene"] for r in nomech)))
for k in ("S0", "S0A", "S0AB"):
    L = [s for s in summ if s["minikel_phase_" + k] == "Launched"]
    c = collections.Counter(s["chembl_supports_minikel_label"] for s in L)
    print("Minikel Launched in %s: %d targets -> %s" % (k, len(L), dict(c)))
    print("    not supported: %s" % ", ".join(sorted(
        s["gene"] for s in L if not s["chembl_supports_minikel_label"].startswith("YES"))))
    P3 = [s for s in summ if s["minikel_phase_" + k] == "Phase III"]
    print("Minikel Phase III in %s: %d targets -> %s ; with an APPROVED CV drug in ChEMBL: %s"
          % (k, len(P3), dict(collections.Counter(s["chembl_supports_minikel_label"]
                                                  for s in P3)),
             ", ".join(sorted(s["gene"] for s in P3
                              if s["chembl_supports_minikel_label"].startswith("YES")))))
