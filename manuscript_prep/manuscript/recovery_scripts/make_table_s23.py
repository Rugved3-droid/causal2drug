#!/usr/bin/env python
"""Supplementary Table S23: for the 11 named targets, the primary site of drug action and
whether a candidate instrument exists in the matching GTEx v8 tissue.

    python manuscript_prep/manuscript/recovery_scripts/make_table_s23.py

Nothing is reclassified. The candidate definition is the recovery map's (p < 5e-8, gene body
+-100 kb); results are read from recovery_out/recovery_gtex_detail.csv. The site of action and
the tissue it is matched to are stated per target with their source; where the regulatory text
does not name an organ this is said, and the site given by the project owner is marked as such.
"""
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.dirname(HERE)
OUT = os.path.join(MS, "display")
D = {(r["gene"], r["tissue"]): r for r in csv.DictReader(
    open(os.path.join(MS, "recovery_out", "recovery_gtex_detail.csv"), encoding="utf-8"))}
N = {r["tissueSiteDetailId"]: int(r["eqtl_samples"]) for r in csv.DictReader(
    open(os.path.join(MS, "revision_sources", "gtex_v8_tissue_sample_sizes.csv"),
         encoding="utf-8"))}
T2 = {r["Target"].split(" ")[0].rstrip("†‡§"): r for r in csv.DictReader(
    open(os.path.join(OUT, "table2_named_targets.csv"), encoding="utf-8-sig"))}
ART = ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial"]
GUT = ["Small_Intestine_Terminal_Ileum", "Colon_Transverse", "Colon_Sigmoid"]
DM = "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid="
# target: (primary site of drug action, matching GTEx v8 tissues, source, note)
SPEC = [
 ("AGTR1", "Vasculature (angiotensin II type 1 receptor on vascular smooth muscle); adrenal "
           "(aldosterone release); kidney", ART + ["Adrenal_Gland", "Kidney_Cortex"],
  "FDA label, losartan, section 12.1: " + DM + "021cd76a-b093-4704-8410-5e7d01e20a54",
  "Kidney was specified by the authors; the label text names vasoconstriction and aldosterone "
  "secretion by the adrenal cortex"),
 ("ANGPTL3", "Liver (ANGPTL3 is 'expressed primarily in the liver'); the antibody binds the "
             "circulating protein", ["Liver"],
  "FDA label, evinacumab (Evkeeza), section 12.1: " + DM
  + "73412138-6d8f-4ea6-bb72-a740190470ff", ""),
 ("APOB", "Liver (antisense to apo B-100 mRNA; activity characterised in hepatocytes)",
  ["Liver"], "FDA label, mipomersen (Kynamro) 2013, section 12.1: https://www.accessdata.fda.gov/"
             "drugsatfda_docs/label/2013/203568s000lbl.pdf",
  "The label describes hepatoma cells and primary hepatocytes in vitro; it does not state the "
  "organ of action in patients"),
 ("APOC3", "Liver (antisense to apoC-III mRNA)", ["Liver"],
  "EMA product information, volanesorsen (Waylivra), section 5.1: https://www.ema.europa.eu/en/"
  "documents/product-information/waylivra-epar-product-information_en.pdf",
  "UNVERIFIED from the regulatory text: section 5.1 does not name the organ. Liver was "
  "specified by the authors"),
 ("AVPR2", "Kidney (V2-receptor antagonism increases urine water excretion)", ["Kidney_Cortex"],
  "FDA label, tolvaptan (Samsca), section 12.1: " + DM
  + "5526617c-c7b9-4556-886d-729bbabbc566",
  "GTEx v8 has eQTLs for kidney cortex only, not medulla"),
 ("F2", "Circulating thrombin (active site), in blood", [],
  "FDA label, argatroban, section 12.1: " + DM + "0cdf3f67-5bf3-4d58-b330-03febdc652bf",
  "The drug acts on a plasma protein; the matching measurement is a plasma pQTL, which exists "
  "(Table 2)"),
 ("PCSK9", "Liver (circulating PCSK9 binds LDL receptors on hepatocytes)", ["Liver"],
  "FDA label, alirocumab (Praluent), section 12.1: " + DM
  + "446f6b5c-0dd4-44ff-9bc2-c2b41f2806b4", "A plasma pQTL also exists (Table 2)"),
 ("REN", "Kidney (renin is secreted by the kidney and acts on circulating angiotensinogen)",
  ["Kidney_Cortex"], "FDA label, aliskiren (Tekturna), section 12.1: " + DM
  + "3dd61fa5-1620-4ebf-bbdb-ede29b92fce2", "A plasma pQTL also exists (Table 2)"),
 ("SLC5A1", "Intestine (SGLT1 inhibition reduces intestinal absorption of glucose and sodium)",
  GUT, "FDA label, sotagliflozin (Inpefa), section 12.1: " + DM
  + "1a46614e-05f6-421a-b6f4-d6f8760d643a",
  "The label states that the mechanism of cardiovascular benefit is not established"),
 ("ITGA2B", "Platelets (glycoprotein IIb/IIIa)", [],
  "FDA label, eptifibatide, section 12.1: " + DM + "00ad7bbe-cadd-4647-b74b-adfc9f2b5e94",
  "GTEx has no platelet tissue. Nearest, whole blood: gene tested, no candidate"),
 ("NPC1L1", "Small intestine (brush border)", ["Small_Intestine_Terminal_Ileum"],
  "FDA label, ezetimibe (Zetia), section 12.1: " + DM
  + "174e45e0-45bb-48b7-9e60-8f6cfc1d22bb", ""),
]


def tis(t):
    return t.replace("_", " ")


rows, n_match, n_cand, cand_names = [], 0, 0, []
for g, site, tissues, src, note in SPEC:
    if not tissues:
        res, any_c = "No matching GTEx tissue", None
    else:
        n_match += 1
        parts, any_c = [], False
        for t in tissues:
            r = D.get((g, t))
            if r is None:
                parts.append("%s (n=%d): gene not tested" % (tis(t), N[t]))
            elif r["instrument"] == "True":
                any_c = True
                parts.append("%s (n=%d): YES, lead p %s" % (tis(t), N[t], r["lead_p"]))
            else:
                parts.append("%s (n=%d): no" % (tis(t), N[t]))
        res = "; ".join(parts)
        if any_c:
            n_cand += 1
            cand_names.append(g)
    rows.append([g, T2[g]["Approved cardiovascular drug(s) in ChEMBL"], site, res,
                 "Yes" if any_c else "No" if any_c is False else "Not applicable", src, note])
title = ("Supplementary Table S23. Site of drug action and candidate instruments in the "
         "matching GTEx tissue, for the 11 named targets")
head = ["Target", "Drug(s)", "Primary site of drug action",
        "Candidate instrument in the matching GTEx v8 tissue (eQTL sample size)",
        "Candidate in a matching tissue", "Source for the site of action", "Note"]
foot = ["A candidate is a variant at p < 5e-8 within the gene body +/- 100 kb in that tissue "
        "(recovery map; rules unchanged). 'no' means the gene was tested in the tissue and no "
        "variant qualified. The site of action is taken from the mechanism-of-action section of "
        "the regulatory label or product information named in the source column, read on "
        "2026-10-11; where that text does not name an organ, the note says so. Matching a site "
        "to GTEx tissues is our own step: vasculature to the three arteries, intestine to "
        "terminal ileum and colon, kidney to kidney cortex. Sample sizes are small and unequal "
        "(Table S20), so 'no' may reflect power.",
        "Summary: %d of 11 targets have a matching GTEx tissue; %d of those have a candidate "
        "instrument in it (%s)." % (n_match, n_cand, ", ".join(cand_names))]
with open(os.path.join(OUT, "tableS23_site_of_action_matching_tissue.csv"), "w", newline="",
          encoding="utf-8-sig") as fh:
    w = csv.writer(fh)
    w.writerow(head)
    w.writerows(rows)

import re as _re                                                   # noqa: E402

_SUP = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def style(v):
    """House style for anything shown to a reader: one notation for P values and windows,
    no ISO dates. Applied to titles, headers, cells and footnotes of every table."""
    if not isinstance(v, str):
        return v
    v = _re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)e-(\d+)\b",
                lambda m: "%s\u00d710%s" % (m.group(1),
                                           ("-" + str(int(m.group(2)))).translate(_SUP)), v)
    v = v.replace("p < 5\u00d710", "P<5\u00d710").replace("p<5\u00d710", "P<5\u00d710")
    v = v.replace("+/- 100 kb", "\u00b1100 kb").replace("+/-100 kb", "\u00b1100 kb")
    v = v.replace("+-100 kb", "\u00b1100 kb").replace("+- 100 kb", "\u00b1100 kb")
    v = v.replace("F >= 10", "F\u226510").replace("MR p < 0.05", "MR P<0.05")
    v = v.replace("FDR < 0.05", "FDR<0.05").replace(", p = ", ", P=").replace("lead p ", "lead P ")
    v = v.replace("accessed 2026-10-11", "accessed October 2026")
    v = v.replace("Accessed 2026-10-11", "Accessed October 2026")
    v = v.replace("on 2026-10-11", "in October 2026").replace("2026-10-11", "October 2026")
    return v

rows = [[style(c) for c in r] for r in rows]
foot = [style(x) for x in foot]
md = ["**%s**" % title, "", "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
md += ["| " + " | ".join(str(c).replace("|", "/") for c in r) + " |" for r in rows]
md += [""] + foot
open(os.path.join(OUT, "tableS23_site_of_action_matching_tissue.md"), "w", encoding="utf-8",
     newline="\n").write("\n".join(md) + "\n")
print("[S23] named targets %d; with a matching GTEx tissue %d; with a candidate instrument in a "
      "matching tissue %d (%s)" % (len(rows), n_match, n_cand, ", ".join(cand_names)))
for r in rows:
    print("[S23] %-7s %s" % (r[0], r[3]))
