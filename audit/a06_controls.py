"""AUDIT step 2d: PCSK9 -> LDL-C and HMGCR -> LDL-C control estimates, from the rule documents.

  outcome study     PREREGISTRATION R2 + amendment R2a (re-applied here to the raw GWAS Catalog files)
  harmonization     PREREGISTRATION R6 (own implementation)
  estimator         R12.2 / RESCUE_CONTROLS_PROTOCOL: lead-variant Wald ratio, SE = |SE_out / beta_exp|
  PCSK9 exposure    RESCUE_CONTROLS_PROTOCOL: UKB-PPP plasma cis-pQTL, the lead qualifying variant of MY recovery map
  HMGCR exposure    eQTLGen whole blood (R3a), lead = smallest p inside gene body +-100 kb that survives
                    liftover and harmonization. eQTLGen publishes Z and N only; beta and SE are derived with the
                    standard conversion beta = z / sqrt(2 f (1-f) (N + z^2)), SE = beta / z  (f from eQTLGen's
                    allele-frequency file). THE RULE DOCUMENTS DO NOT STATE THIS CONVERSION - see AUDIT.md.

Uses my own remote tabix reader (audit/tabix_region.py); nothing is imported from the project's src/.
Run:  python audit/a06_controls.py
Writes audit/out/audit_controls.json
"""
import csv, gzip, json, math, os, re, subprocess, sys, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabix_region as tx
from pyliftover import LiftOver

ROOT = r"D:\causal2drug"; RAW = os.path.join(ROOT, "data", "raw"); OUT = os.path.join(ROOT, "audit", "out")
FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
ANC = ["african", "asian", "hispanic", "latin", "admixed", "native", "oceani", "arab", "middle east", "aborigin",
       "other", "indian", "japanese", "chinese", "korean", "greater", "sub-saharan", "caribbean", "pacific", "nr"]


def pnorm2(z):
    return math.erfc(abs(z) / math.sqrt(2))


# ---- R2 + R2a --------------------------------------------------------------
def select_gwas(efo):
    harm = {}
    for line in open(os.path.join(RAW, "harmonised_list.txt")):
        m = re.search(r"(GCST\d+)\.h\.tsv\.gz$", line.strip())
        if m:
            harm[m.group(1)] = line.strip().lstrip("./")
    cands = []
    with open(os.path.join(RAW, "gwas-catalog-download-studies-v1.0.3.1.txt"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if efo not in (r["MAPPED_TRAIT_URI"] or ""):
                continue
            acc = r["STUDY ACCESSION"]
            if acc not in harm:
                continue
            ss = r["INITIAL SAMPLE SIZE"] or ""
            low = ss.lower()
            if "european" not in low or any(re.search(r"\b" + a, low) for a in ANC):
                continue                                     # R2a: European ancestry only
            ints = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", ss)]
            if ints:
                cands.append((-max(ints), int(acc[4:]), acc, ss))
    cands.sort()
    for negn, _, acc, ss in cands:                           # largest N first, ties by lowest accession
        url = FTP + harm[acc]
        try:
            req = urllib.request.Request(url + ".tbi", method="HEAD")
            ok = urllib.request.urlopen(req, timeout=60).status == 200
        except Exception:
            ok = False
        if ok:
            return acc, -negn, ss, url, len(cands)
    return None


def harmonize(g, ea_exp, oa_exp, f_exp):
    """R6. g: outcome row dict. Exposure effect allele ea_exp / other oa_exp. Returns (beta_out aligned, status)."""
    ge, go = g["effect_allele"].upper(), g["other_allele"].upper()
    if len(ge) != 1 or len(go) != 1 or len(ea_exp) != 1 or len(oa_exp) != 1:
        return None, "indel"
    b = float(g["beta"])
    pal = {ea_exp, oa_exp} in ({"A", "T"}, {"C", "G"})
    if pal:
        try:
            fg = float(g["effect_allele_frequency"]); fe = float(f_exp)
        except (TypeError, ValueError):
            return None, "ambiguous_unresolved"
        # put both frequencies on the exposure effect allele, by label
        if (ge, go) == (ea_exp, oa_exp):
            same = True
        elif (ge, go) == (oa_exp, ea_exp):
            same = False
        else:
            return None, "allele_mismatch"
        lo = lambda x: x < 0.42
        hi = lambda x: x > 0.58
        if (lo(fg) and lo(fe)) or (hi(fg) and hi(fe)):
            return (b if same else -b), "palindromic_resolved_" + ("aligned" if same else "flipped")
        if (lo(fg) and hi(fe)) or (hi(fg) and lo(fe)):
            return (-b if same else b), "palindromic_resolved_strand_" + ("flipped" if same else "aligned")
        return None, "ambiguous_unresolved"
    if (ge, go) == (ea_exp, oa_exp):
        return b, "aligned"
    if (ge, go) == (oa_exp, ea_exp):
        return -b, "flipped"
    ce, co = COMP.get(ge), COMP.get(go)
    if (ce, co) == (ea_exp, oa_exp):
        return b, "strand_flipped_aligned"
    if (ce, co) == (oa_exp, ea_exp):
        return -b, "strand_flipped_flipped"
    return None, "allele_mismatch"


def wald(b_out, se_out, b_exp, se_exp):
    beta = b_out / b_exp; se = abs(se_out / b_exp)
    return dict(mr_beta=beta, mr_se=se, ci_lo=beta - 1.96 * se, ci_hi=beta + 1.96 * se, mr_z=beta / se,
                mr_p=pnorm2(beta / se), F=(b_exp / se_exp) ** 2)


res = {}
sel = select_gwas("EFO_0004611")
acc, n, ss, url, ncand = sel
print(f"[D01] R2+R2a for EFO_0004611: {acc} (N={n:,}; '{ss}'); {ncand} European-only harmonised candidates")
hdr = tx.header(url)
ix = tx.read_index(url)
res["outcome"] = dict(accession=acc, n=n, sample=ss, url=url)


def outcome_rows(chrom, pos):
    lines, _ = tx.fetch(url, str(chrom), pos, pos, ix)
    return [dict(zip(hdr, l.split("\t"))) for l in lines]


# ---- PCSK9 (UKB-PPP pQTL from my recovery map) ------------------------------
rm = {r["symbol"]: r for r in csv.DictReader(open(os.path.join(OUT, "audit_recovery_map.csv"), encoding="utf-8"))}
r = rm["PCSK9"]
c, pos37, a0, a1 = r["ukb_lead"].split(":")[:4]
b_exp, se_exp = float(r["ukb_beta"]), float(r["ukb_se"])
# A1 frequency for the palindrome rule (not needed unless palindromic): read back from the sheet
g = outcome_rows(c, int(r["ukb_pos38"]))
print(f"[D02] PCSK9 exposure: UKB-PPP {r['ukb_sheet']} {r['ukb_lead']} hg38 pos {r['ukb_pos38']} "
      f"beta(A1={a1}) {b_exp} se {se_exp} log10p {r['ukb_log10p']} F {(b_exp/se_exp)**2:.1f}; outcome rows at position: {len(g)}")
out = None
for row in g:
    bo, st = harmonize(row, a1, a0, None)
    if bo is not None:
        out = dict(wald(bo, float(row["standard_error"]), b_exp, se_exp), harmonization=st, outcome_rsid=row.get("rsid"),
                   outcome_beta_aligned=bo, outcome_se=float(row["standard_error"]), outcome_p=row["p_value"],
                   outcome_ea=row["effect_allele"], outcome_oa=row["other_allele"], exposure_variant=r["ukb_lead"],
                   exposure_sheet=r["ukb_sheet"], exposure_beta=b_exp, exposure_se=se_exp)
        break
res["PCSK9"] = out or dict(status="not harmonizable", tried=[(x["effect_allele"], x["other_allele"]) for x in g])
if out:
    ok = out["F"] >= 10 and out["mr_p"] < 0.05 and out["mr_beta"] > 0
    out["verdict"] = "PASS" if ok else "FAIL"
    print(f"[D03] PCSK9 -> LDL-C  Wald beta {out['mr_beta']:+.4f} (SE {out['mr_se']:.4f}; 95% CI {out['ci_lo']:.4f} to {out['ci_hi']:.4f}); "
          f"p {out['mr_p']:.3g}; F {out['F']:.1f}; {out['harmonization']}; expected sign +; {out['verdict']}")

# ---- HMGCR (eQTLGen) -------------------------------------------------------
tg = {x["symbol"]: x for x in csv.DictReader(open(os.path.join(OUT, "audit_per_target.csv"), encoding="utf-8"))}["HMGCR"]
lo_, hi_ = int(tg["start"]) - 100_000, int(tg["end"]) + 100_000
cands = []
with open(os.path.join(OUT, "eqtlgen_target_rows.tsv"), encoding="utf-8") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = dict(zip(h, line.rstrip("\n").split("\t")))
        if f["Gene"] == tg["ensembl_id"] and f["SNPChr"] == tg["chrom"] and float(f["Pvalue"]) < 5e-8 \
                and lo_ <= int(f["SNPPos"]) <= hi_:
            cands.append(f)
cands.sort(key=lambda f: (float(f["Pvalue"]), -abs(float(f["Zscore"]))))
print(f"[D04] HMGCR: {len(cands)} eQTLGen variants at p<5e-8 in {tg['chrom']}:{lo_}-{hi_} (GRCh37); "
      f"top by p: {[(f['SNP'], f['Pvalue'], f['Zscore']) for f in cands[:3]]}")
need = {f["SNP"] for f in cands[:40]}
af = {}
with gzip.open(os.path.join(RAW, "eqtlgen_SNP_AF.txt.gz"), "rt") as fh:
    next(fh)
    for line in fh:
        if line[:line.index("\t")] in need:
            f = line.rstrip("\n").split("\t")
            af[f[0]] = (f[3], f[4], float(f[8]))            # AlleleA, AlleleB, freq of AlleleB
            if len(af) == len(need):
                break
lifto = LiftOver(os.path.join(RAW, "hg19ToHg38.over.chain.gz"))
out = None; log = []
for f in cands[:40]:
    lifted = lifto.convert_coordinate("chr" + f["SNPChr"], int(f["SNPPos"]) - 1)
    if not lifted:
        log.append((f["SNP"], "liftover_failed")); continue
    pos38 = lifted[0][1] + 1
    ea, oa = f["AssessedAllele"].upper(), f["OtherAllele"].upper()
    z, n_ = float(f["Zscore"]), int(f["NrSamples"])
    if f["SNP"] not in af:
        log.append((f["SNP"], "no_allele_frequency")); continue
    aA, aB, fB = af[f["SNP"]]
    fa = fB if aB == ea else 1 - fB if aA == ea else None
    if fa is None:
        log.append((f["SNP"], "af_allele_mismatch")); continue
    b_exp = z / math.sqrt(2 * fa * (1 - fa) * (n_ + z * z)); se_exp = b_exp / z
    got = None
    for row in outcome_rows(f["SNPChr"], pos38):
        bo, st = harmonize(row, ea, oa, fa)
        if bo is not None:
            got = (row, bo, st); break
        log.append((f["SNP"], st))
    if not got:
        log.append((f["SNP"], "absent_or_not_harmonizable")); continue
    row, bo, st = got
    out = dict(wald(bo, float(row["standard_error"]), b_exp, se_exp), harmonization=st, exposure_variant=f["SNP"],
               pos37=int(f["SNPPos"]), pos38=pos38, exposure_z=z, exposure_n=n_, exposure_p=f["Pvalue"],
               exposure_ea=ea, exposure_eaf=fa, exposure_beta=b_exp, exposure_se=se_exp,
               outcome_beta_aligned=bo, outcome_se=float(row["standard_error"]), outcome_p=row["p_value"],
               outcome_ea=row["effect_allele"], outcome_oa=row["other_allele"], skipped_before=log[:])
    # scale-free form, independent of the Z -> beta conversion
    out["outcome_z_aligned"] = bo / float(row["standard_error"])
    break
res["HMGCR"] = out or dict(status="no harmonizable instrument", log=log)
if out:
    ok = out["F"] >= 10 and out["mr_p"] < 0.05 and out["mr_beta"] > 0
    out["verdict"] = "PASS" if ok else "FAIL"
    print(f"[D05] HMGCR lead {out['exposure_variant']} (GRCh37 {out['pos37']} -> GRCh38 {out['pos38']}), z {out['exposure_z']}, "
          f"N {out['exposure_n']}, EAF {out['exposure_eaf']:.4f}, beta {out['exposure_beta']:.5f}, se {out['exposure_se']:.5f}; "
          f"skipped before it: {out['skipped_before']}")
    print(f"[D06] HMGCR -> LDL-C  Wald beta {out['mr_beta']:+.4f} (SE {out['mr_se']:.4f}; 95% CI {out['ci_lo']:.4f} to {out['ci_hi']:.4f}); "
          f"z {out['mr_z']:.2f}; p {out['mr_p']:.3g}; F {out['F']:.1f}; {out['harmonization']}; expected sign +; {out['verdict']}")
json.dump(res, open(os.path.join(OUT, "audit_controls.json"), "w"), indent=1, default=str)
