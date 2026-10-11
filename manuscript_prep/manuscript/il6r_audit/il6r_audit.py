#!/usr/bin/env python
"""Audit of the IL6R -> CRP positive control. Read-only with respect to the pipeline outputs.

    python manuscript_prep/manuscript/il6r_audit/il6r_audit.py > .../il6r_audit_output.txt

1. eQTLGen rows for IL6R (ENSG00000160712): the lead variant, the variant the MR used,
   rs2228145 (Asp358Ala) and its proxies.
2. LD in 1000 Genomes phase 3 EUR (the pipeline's own panel and client, GRCh37).
3. CRP associations for the same variants in the outcome file the control used.
4. Re-derivation of the harmonization of the MR variant by hand.
"""
import csv
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import d1b_clump as C               # noqa: E402
import tabix_http as T              # noqa: E402
from harmonize import harmonize_one  # noqa: E402
from pyliftover import LiftOver     # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
RAW = os.path.join(ROOT, "data", "raw")
ENSG = "ENSG00000160712"
GENE = ("1", 154377669, 154441926)                       # GRCh37.87, from the provenance file
FOCUS = ["rs10908839", "rs2229238", "rs2228145", "rs4129267", "rs4537545", "rs12126142",
         "rs58229059"]
URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90179001-GCST90180000/GCST90179146/harmonised/GCST90179146.h.tsv.gz")

prov = json.load(open(os.path.join(ROOT, "deliverables", "d1c", "IL6R_CRP_provenance.json")))
print("[A01] control as run: eQTLGen lead %s (assessed %s, z %+.2f, p %.3g); MR variant %s; "
      "MR beta %+.4f, p %.3g; expected %s, observed %s; verdict %s"
      % (prov["lead_snp"]["rsid"], prov["lead_snp"]["assessed"], prov["lead_snp"]["z"],
         prov["lead_snp"]["p_exp"], prov["mr"]["lead_snp"], prov["mr"]["beta"],
         prov["mr"]["p"], prov["mr"]["expected_direction"], prov["mr"]["direction"],
         prov["VERDICT"]))
print("[A01] harmonization counts in the control: %s" % prov["harmonization"])

# ---------------------------------------------------------------- 1. eQTLGen rows
cache = os.path.join(HERE, "eqtlgen_IL6R_rows.tsv")
if not os.path.exists(cache):
    import glob
    full = [f for f in glob.glob(os.path.join(RAW, "*cis-eQTLsFDR-ProbeLevel*")) if
            f.endswith(".gz")][0]
    with gzip.open(full, "rt", encoding="utf-8") as fh, \
            open(cache, "w", encoding="utf-8", newline="") as out:
        head = fh.readline()
        out.write(head)
        gi = head.rstrip("\n").split("\t").index("Gene")
        for line in fh:
            if ENSG in line and line.split("\t")[gi] == ENSG:
                out.write(line)
rows = list(csv.DictReader(open(cache, encoding="utf-8"), delimiter="\t"))
for r in rows:
    r["z"], r["p"], r["pos"] = float(r["Zscore"]), float(r["Pvalue"]), int(r["SNPPos"])
inwin = [r for r in rows if r["SNPChr"] == GENE[0]
         and GENE[1] - 100000 <= r["pos"] <= GENE[2] + 100000]
sig = sorted([r for r in inwin if r["p"] < 5e-8], key=lambda r: r["p"])
print("\n[A02] eQTLGen full file: %d rows for IL6R; %d within gene body +-100 kb; %d of those "
      "at p < 5e-8" % (len(rows), len(inwin), len(sig)))
print("[A02] five strongest within the window:")
for r in sig[:5]:
    print("        %-11s chr%s:%d assessed %s other %s z %+7.2f p %.3g N %s"
          % (r["SNP"], r["SNPChr"], r["pos"], r["AssessedAllele"], r["OtherAllele"], r["z"],
             r["p"], r["NrSamples"]))
byrs = {r["SNP"]: r for r in rows}
print("[A03] variants of interest in eQTLGen (IL6R expression):")
for rs in FOCUS:
    r = byrs.get(rs)
    if r is None:
        print("        %-11s NOT in the eQTLGen full file for IL6R" % rs)
    else:
        print("        %-11s chr%s:%d assessed %s other %s z %+7.2f p %.3g rank among window "
              "hits %s" % (rs, r["SNPChr"], r["pos"], r["AssessedAllele"], r["OtherAllele"],
                           r["z"], r["p"],
                           (sig.index(r) + 1) if r in sig else "not significant at 5e-8"))

# ---------------------------------------------------------------- 2. LD, 1000G EUR GRCh37
keep = C.eur_samples()
ctx = C.ChromCtx(GENE[0], keep)
lo, hi = GENE[1] - 100000, GENE[2] + 100000
vcf = {}
# matched on GRCh37 position and allele pair, as the pipeline's clumping does
WANT = {(byrs[rs]["pos"], frozenset((byrs[rs]["AssessedAllele"], byrs[rs]["OtherAllele"]))): rs
        for rs in FOCUS if rs in byrs}
nrec = 0
for raw in ctx.region(lo, hi):
    if not raw or raw[:1] == b"#":
        continue
    f = raw.split(b"\t")
    if len(f) != 9 + ctx.n_samples_total:
        continue
    nrec += 1
    try:
        pos_ = int(f[1])
    except ValueError:
        continue
    rsid = WANT.get((pos_, frozenset((f[3].decode(), f[4].decode()))))
    if rsid is None:
        continue
    blob = b"".join(g[:3] for g in f[9:])
    arr = np.frombuffer(blob, dtype="S1").reshape(-1, 3)
    d = (arr[:, 0] == b"1").astype(float) + (arr[:, 2] == b"1").astype(float)
    vcf[rsid] = dict(pos=int(f[1]), ref=f[3].decode(), alt=f[4].decode(), d=d[ctx.idx])
print("\n[A04] 1000 Genomes phase 3, %d EUR samples (pipeline panel), GRCh37; %d complete records "
      "read in the window:" % (ctx.n_eur, nrec))
for rs in FOCUS:
    v = vcf.get(rs)
    print("        %-11s %s" % (rs, "not found in the panel" if v is None else
                                "chr1:%d REF %s ALT %s, ALT frequency %.3f"
                                % (v["pos"], v["ref"], v["alt"], v["d"].mean() / 2)))
ref = vcf.get("rs2228145")
print("[A05] LD with rs2228145 (Asp358Ala; ALT = %s):" % (ref["alt"] if ref else "?"))
ld = {}
for rs in FOCUS:
    v = vcf.get(rs)
    if v is None or ref is None or rs == "rs2228145":
        continue
    r = float(np.corrcoef(v["d"], ref["d"])[0, 1])
    ld[rs] = r
    print("        %-11s r = %+.3f, r2 = %.3f; the %s allele of %s travels with the %s allele "
          "of rs2228145" % (rs, r, r * r, v["alt"], rs, ref["alt"] if r > 0 else ref["ref"]))

# ---------------------------------------------------------------- 3. CRP outcome rows
chain = LiftOver(os.path.join(RAW, "hg19ToHg38.over.chain.gz"))
lines, _ = T.header_lines(URL, nbytes=1 << 16)
gcols = lines[0].rstrip("\r").split("\t")
gi = {k: i for i, k in enumerate(gcols)}
p38 = {}
for rs in FOCUS:
    pos37 = byrs[rs]["pos"] if rs in byrs else (vcf[rs]["pos"] if rs in vcf else None)
    if pos37:
        res = chain.convert_coordinate("chr1", pos37 - 1)
        p38[rs] = res[0][1] + 1 if res else None
grows, _ = T.query(URL, "1", min(v for v in p38.values() if v) - 10,
                   max(v for v in p38.values() if v) + 10)
gmap = {}
for g in grows:
    gmap.setdefault(int(g[gi["base_pair_location"]]), []).append(g)
print("\n[A06] CRP outcome file GCST90179146 (GRCh38), columns: %s" % gcols[:8])
crp = {}
for rs in FOCUS:
    g = (gmap.get(p38.get(rs)) or [None])[0]
    if g is None:
        print("        %-11s GRCh38 pos %s: NOT in the outcome file" % (rs, p38.get(rs)))
        continue
    crp[rs] = dict(ea=g[gi["effect_allele"]], oa=g[gi["other_allele"]],
                   beta=float(g[gi["beta"]]), se=float(g[gi["standard_error"]]),
                   p=g[gi["p_value"]], eaf=g[gi["effect_allele_frequency"]])
    print("        %-11s chr1:%d effect allele %s other %s beta %+.4f SE %.4f p %s EAF %s"
          % (rs, p38[rs], crp[rs]["ea"], crp[rs]["oa"], crp[rs]["beta"], crp[rs]["se"],
             crp[rs]["p"], crp[rs]["eaf"]))

# ---------------------------------------------------------------- 4. one table, one allele
print("\n[A07] Both associations expressed for the SAME allele of each variant (the eQTLGen "
      "assessed allele), and the single-variant Wald ratio:")
print("        %-11s %-7s %-22s %-22s %-10s %s" % ("variant", "allele", "IL6R expression z",
                                                   "CRP beta (SE)", "Wald", "r2 with rs2228145"))
for rs in FOCUS:
    e, o = byrs.get(rs), crp.get(rs)
    if e is None or o is None:
        continue
    h = harmonize_one(e["OtherAllele"], e["AssessedAllele"], o["ea"], o["oa"], o["beta"])
    if h["sign"] is None:
        print("        %-11s harmonization: %s" % (rs, h["status"]))
        continue
    sgn = "same sign" if (e["z"] > 0) == (h["beta_out_harmonized"] > 0) else "opposite sign"
    print("        %-11s %-7s %+7.2f (p %.2g)%8s %+.4f (%.4f)%8s %-10s %s  [%s]"
          % (rs, e["AssessedAllele"], e["z"], e["p"], "", h["beta_out_harmonized"], o["se"], "",
             "positive" if sgn == "same sign" else "negative",
             ("%.3f" % (ld[rs] ** 2)) if rs in ld else ("1 (itself)" if rs == "rs2228145"
                                                        else "n/a"), h["status"]))

# ---------------------------------------------------------------- 5. the MR variant by hand
m = prov["mr"]
e, o = byrs[m["lead_snp"]], crp[m["lead_snp"]]
print("\n[A08] hand check of the MR variant %s: eQTLGen assessed allele %s (z %+.2f); outcome "
      "effect allele %s, beta %+.5f, EAF %s. The outcome effect allele is the eQTLGen %s "
      "allele, so the outcome beta for the assessed allele is %+.5f. Pipeline value: %+.5f. "
      "Pipeline status: %s."
      % (m["lead_snp"], e["AssessedAllele"], e["z"], o["ea"], o["beta"], o["eaf"],
         "assessed" if o["ea"] == e["AssessedAllele"] else "other",
         o["beta"] if o["ea"] == e["AssessedAllele"] else -o["beta"],
         m["beta_outcome_harmonized"], m["harmonization_of_lead"]))
mer = {r["rsid"]: r for r in csv.DictReader(open(os.path.join(
    ROOT, "deliverables", "d1c", "IL6R_CRP_harmonized_merge.tsv"), encoding="utf-8"),
    delimiter="\t")}
x = mer[m["lead_snp"]]
print("[A08] allele-frequency check: eQTLGen frequency of the assessed allele %s = %.3f; "
      "outcome frequency of its effect allele %s = %s, i.e. %.3f for the assessed allele."
      % (x["exp_assessed"], float(x["af_assessed_exp"]), x["gwas_effect"], x["af_effect_out"],
         float(x["af_effect_out"]) if x["gwas_effect"] == x["exp_assessed"]
         else 1 - float(x["af_effect_out"])))
print("[A09] why the eQTLGen lead %s was not the MR variant: %s"
      % (prov["lead_snp"]["rsid"],
         "it has no row at its lifted position in the outcome file" if
         prov["lead_snp"]["rsid"] not in mer else "it is in the merge file"))
same = sum(1 for r in mer.values() if (float(r["beta_exp"]) > 0) == (float(r["beta_out_harm"]) > 0))
print("[A10] of the %d harmonized variants in the control's merge file, %d give a positive "
      "single-variant ratio and %d a negative one" % (len(mer), same, len(mer) - same))
