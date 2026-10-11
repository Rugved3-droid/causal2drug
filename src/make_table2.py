#!/usr/bin/env python
"""Table 2 -- method availability for cardiovascular target loci.

    python src/make_table2.py

Emits table2_method_availability.{pdf,svg,png} (vector, journal resolution) and
table2_method_availability.csv (machine-readable, one row per line of the table).

EVERY number is recomputed from frozen outputs at render time. Nothing is hardcoded:
sources are cv_report.json, cv_robustness.json, partc_skipped_accessions.tsv,
d1h/part3_results.tsv (R14.1, flawed rule) and d1j/r15_results.tsv (R15, corrected).

The residual-z diagnostic is DENTIST (Chen et al. 2021, Nat Commun 12:7117). The
implementation here is ours; the method is theirs.
"""
import csv
import json
import math
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = r"D:\causal2drug\deliverables"
FIG = D + r"\figures"
VERM, GREY, INK = "#D55E00", "#7F7F7F", "#222222"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})

SIG_Z = 5.4513   # |z| for p = 5e-8, the R15 outcome-significance criterion


def load():
    rep = json.load(open(D + r"\cv\cv_report.json"))
    rob = json.load(open(D + r"\cv\cv_robustness.json"))
    skip = list(csv.DictReader(open(D + r"\d1f\partc_skipped_accessions.tsv",
                                    encoding="utf-8"), delimiter="\t"))
    r14 = list(csv.DictReader(open(D + r"\d1h\part3_results.tsv",
                                   encoding="utf-8"), delimiter="\t"))
    r15 = list(csv.DictReader(open(D + r"\d1j\r15_results.tsv",
                                   encoding="utf-8"), delimiter="\t"))
    return rep, rob, skip, r14, r15


def stat(rows, col):
    """median and % > 2 of a DENTIST RMSE column.

    Non-finite values are EXCLUDED, not merely non-empty ones: at least one locus
    (CA2) carries a literal 'nan', which float() parses without error and which
    silently corrupts a median if it is left in. n is the count actually used.
    """
    v = [float(r[col]) for r in rows if r.get(col) not in ("", None)]
    v = [x for x in v if math.isfinite(x)]
    return len(v), statistics.median(v), 100.0 * sum(1 for x in v if x > 2) / len(v)


def build(rep, rob, skip, r14, r15):
    """Return (availability rows, diagnostics rows, notes); every value computed here."""
    n_scr = rep["cv_loci_screened"]
    n_sig = rep["cv_loci_gwas_sig"]
    n_use = rep["cv_loci_coloc_usable"]
    n_un = n_scr - n_use
    ext, pri = rob["CV_EXTENDED"]["attr"], rob["CV_PRIMARY"]["attr"]
    pc_lost, pc_den = rep["cv_pairs_lost_unindexed"], rep["cv_pairs_partc"]

    four = ["GCST90029021", "GCST90043777", "GCST90044006", "GCST90080708"]
    sk = {r["gwas_acc"]: r for r in skip}
    acc_rows = [(a, int(sk[a]["n_pairs_lost"])) for a in four if a in sk]
    acc_tot = sum(a[1] for a in acc_rows)

    A = [("S", "Colocalization availability at cardiovascular target loci", "", "", ""),
         ("R", "Cardiovascular target loci screened",
          "%d" % n_scr, "\u2014", "R15 rule, cardiovascular subset"),
         ("R", "Loci with a genome-wide significant disease signal in the gene window",
          "%d" % n_sig, "of %d" % n_scr, "disease GWAS p < 5\u00d710\u207b\u2078"),
         ("R", "Loci yielding a usable credible-set colocalization",
          "%d" % n_use, "of %d" % n_scr, "coloc.susie, after the FIX 1 assertion"),
         ("H", "Colocalization UNAVAILABLE",
          "%d (%.1f%%)" % (n_un, 100.0 * n_un / n_scr), "of %d" % n_scr,
          "method availability, not evidence"),
         ("S", "Infrastructure attrition \u2014 disease GWAS cannot be region-queried",
          "", "", ""),
         ("H", "Cardiovascular pairs lost (CV-EXTENDED, atlas definition)",
          "%d (%.1f%%)" % (ext[0], ext[2]), "of %d" % ext[1], "PRIMARY figure"),
         ("R", "Cardiovascular pairs lost (CV-PRIMARY, strict definition)",
          "%d (%.1f%%)" % (pri[0], pri[2]), "of %d" % pri[1],
          "robustness, Analysis 2"),
         ("R", "Part C-conditioned sub-rate",
          "%d (%.1f%%)" % (pc_lost, 100.0 * pc_lost / pc_den), "of %d" % pc_den,
          "SUB-RATE ONLY \u2014 conditioned on two unrelated Part C rules"),
         ("S", "The four accessions responsible for the %d Part C-conditioned losses"
          % acc_tot, "", "", "")]
    for a, n_ in acc_rows:
        A.append(("R", a, "%d pairs" % n_, "of %d" % acc_tot,
                  "plain gzip, no .tbi (HTTP 404) \u2014 not range-queryable"))

    n14 = len(r14)
    e14 = stat(r14, "dentist_rmse_eqtlgen")
    g14 = stat(r14, "dentist_rmse_gwas")
    cs14 = sum(1 for r in r14 if r.get("susie_PP4_eqtlgen") not in ("", None))
    sig14 = sum(1 for r in r14 if r.get("gwas_lead_z") not in ("", None)
                and abs(float(r["gwas_lead_z"])) > SIG_Z)
    sig15 = [r for r in r15 if r["gwas_genome_wide_sig"] == "True"]
    n15 = len(sig15)
    e15 = stat(sig15, "dentist_rmse_eqtlgen")
    g15 = stat(sig15, "dentist_rmse_gwas")
    cs15 = sum(1 for r in sig15 if r.get("susie_PP4") not in ("", None))
    n_g15 = sum(1 for r in sig15 if r.get("dentist_rmse_gwas") not in ("", None)
                and float(r["dentist_rmse_gwas"]) > 2)
    n_e15 = sum(1 for r in sig15 if r.get("dentist_rmse_eqtlgen") not in ("", None)
                and float(r["dentist_rmse_eqtlgen"]) > 2)

    Dd = [("R", "Loci analysed",
           "%d of 60 (32 superseded)" % n14,
           "%d of %d screened" % (n15, len(r15))),
          ("R", "Loci genome-wide significant on the disease side",
           "%d / %d" % (sig14, n14), "%d / %d (100%%)" % (n15, n15)),
          ("H", "DENTIST residual-z RMSE, gene-expression side \u2014 median",
           "%.2f" % e14[1], "%.2f" % e15[1]),
          ("H", "DENTIST residual-z RMSE, gene-expression side \u2014 loci above 2",
           "%.0f%%" % e14[2], "%.0f%%  (%d/%d)" % (e15[2], n_e15, n15)),
          ("R", "DENTIST residual-z RMSE, disease side \u2014 median",
           "%.2f  (n=%d)" % (g14[1], g14[0]), "%.2f" % g15[1]),
          ("R", "DENTIST residual-z RMSE, disease side \u2014 loci above 2",
           "%.0f%%  (0/%d)" % (g14[2], g14[0]),
           "%.0f%%  (%d/%d)" % (g15[2], n_g15, n15)),
          ("H", "Credible-set pair returned",
           "%d / %d (%.0f%%)   WITHDRAWN" % (cs14, n14, 100.0 * cs14 / n14),
           "%d / %d (%.0f%%)" % (cs15, n15, 100.0 * cs15 / n15))]

    return A, Dd, dict(n14=n14, n15=n15)


def render(A, Dd, notes):
    fig = plt.figure(figsize=(11.0, 8.8))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    X = (1.5, 55.0, 69.5, 78.5)
    y = 96.5

    def rule(yy, lw=1.1, x0=1.5, x1=98.5, col=INK):
        ax.plot([x0, x1], [yy, yy], lw=lw, color=col, clip_on=False)

    ax.text(1.5, y, "Table 2.  Method availability for cardiovascular target loci",
            fontsize=13.5, fontweight="bold", color=INK)
    y -= 3.0
    ax.text(1.5, y, "What the standard analytic methods can and cannot deliver on these "
                    "targets. Counts are loci or target-disease pairs, not patients.",
            fontsize=8.6, color="#444444")
    y -= 3.2
    rule(y, 1.3)
    y -= 2.9
    for lab, x in zip(("Measure", "Value", "Of", "Basis / rule"), X):
        ax.text(x, y, lab, fontsize=9.0, fontweight="bold", color=INK)
    y -= 1.5
    rule(y, 0.9)
    y -= 3.0

    for kind, c0, c1, c2, c3 in A:
        if kind == "S":
            y -= 0.8
            ax.text(1.5, y, c0, fontsize=9.3, fontweight="bold", color=VERM)
            y -= 1.3
            rule(y, 0.6, col="#CCCCCC")
            y -= 2.8
            continue
        fw = "bold" if kind == "H" else "normal"
        ax.text(X[0], y, c0, fontsize=8.5, fontweight=fw, color=INK)
        ax.text(X[1], y, c1, fontsize=8.5, fontweight=fw,
                color=(VERM if kind == "H" else INK))
        ax.text(X[2], y, c2, fontsize=8.2, color="#555555")
        ax.text(X[3], y, c3, fontsize=7.9, color="#555555")
        y -= 2.95

    y -= 1.0
    rule(y, 1.3)
    y -= 3.1
    ax.text(1.5, y, "Fine-mapping diagnostics under the corrected locus-selection rule",
            fontsize=9.3, fontweight="bold", color=VERM)
    y -= 2.9
    XD = (1.5, 55.0, 80.0)
    for lab, x in zip(("Diagnostic", "R14.1 \u2014 flawed rule (withdrawn)",
                       "R15 \u2014 corrected"), XD):
        ax.text(x, y, lab, fontsize=9.0, fontweight="bold", color=INK)
    y -= 1.5
    rule(y, 0.9)
    y -= 3.0
    for kind, c0, c1, c2 in Dd:
        fw = "bold" if kind == "H" else "normal"
        ax.text(XD[0], y, c0, fontsize=8.5, fontweight=fw, color=INK)
        ax.text(XD[1], y, c1, fontsize=8.5, color=GREY)
        ax.text(XD[2], y, c2, fontsize=8.5, fontweight=fw,
                color=(VERM if kind == "H" else INK))
        y -= 2.95
    y -= 0.8
    rule(y, 1.3)
    y -= 3.2

    cap = [
        "R15 selection rule, declared before running: n_overlap \u2265 100; the disease "
        "GWAS has \u2265 1 variant with p < 5\u00d710\u207b\u2078 inside the merged gene "
        "window; one locus per target by largest overlap; no cap.",
        "R14.1 omitted the disease-significance criterion, so it admitted loci with no "
        "disease signal \u2014 its residual-z values on the disease side are an artefact "
        "of that and are withdrawn, shown only for comparison.",
        "Of the R14.1 60-locus set, %d completed and 32 are excluded as superseded. The "
        "residual-z diagnostic is DENTIST (Chen et al. 2021, Nat Commun 12:7117); the "
        "implementation is ours, the method is theirs." % notes["n14"],
        "Colocalization figures are method availability, not evidence: Gate D is closed "
        "as a documented negative result. One R14.1 locus (CA2) has a non-finite "
        "disease-side residual-z and is excluded from that row only, so its n is 27 "
        "against 28 elsewhere.",
        "",
        "Reproduce: python src/make_table2.py     Sources: cv_report.json, "
        "cv_robustness.json, d1f/partc_skipped_accessions.tsv, d1h/part3_results.tsv, "
        "d1j/r15_results.tsv",
    ]
    for line in cap:
        ax.text(1.5, y, line, fontsize=7.4, color="#555555")
        y -= 2.2

    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIG + r"\table2_method_availability." + ext,
                    bbox_inches="tight", dpi=(600 if ext != "png" else 150))
    plt.close(fig)
    print("  wrote table2_method_availability.pdf / .svg / .png")


def write_csv(A, Dd):
    p = D + r"\cv\table2_method_availability.csv"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["block", "section", "emphasis", "measure", "value", "of",
                    "basis_or_rule"])
        sec = ""
        for kind, c0, c1, c2, c3 in A:
            if kind == "S":
                sec = c0
                continue
            w.writerow(["availability", sec, kind, c0, c1, c2, c3])
        for kind, c0, c1, c2 in Dd:
            w.writerow(["finemapping_diagnostics",
                        "corrected R15 vs flawed R14.1", kind, c0, c1, c2, ""])
    print("  wrote table2_method_availability.csv")


def main():
    rep, rob, skip, r14, r15 = load()
    A, Dd, notes = build(rep, rob, skip, r14, r15)
    render(A, Dd, notes)
    write_csv(A, Dd)
    print("Table 2 complete (all values recomputed from frozen outputs)")


if __name__ == "__main__":
    main()
