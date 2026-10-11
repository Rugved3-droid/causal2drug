#!/usr/bin/env python
"""Manuscript figures. Each figure is written as PDF, SVG (vector) and PNG (300 dpi).

    python manuscript_prep/analysis/classify.py     # first
    python manuscript_prep/figures/make_figures.py            # all four
    python manuscript_prep/figures/make_figures.py fig2       # one

  fig1   Figure 1   derivation of the target set and instrument status (flow diagram)
  fig2   Figure 2   instrument status by highest clinical phase
  fig3   Suppl. Figure S2   launched and Phase III targets with no blood instrument
         (was Figure 3 until 2026-10-11; Figure 3 is now the recovery map, built by
         manuscript/recovery_scripts/make_recovery_outputs.py)
  figS1  Suppl. Figure S1   Class C pre-registered test (negative result)

Reads manuscript_prep/analysis/out/ and, for figS1 only,
deliverables/d1h/part3_headtohead.tsv (see DEVIATION_CLASS_C_FIGURE.md).

RULES APPLIED (owner, 2026-10-10)
  * "No eQTL data" is shown with its sub-types: autosomal gene not tested (solid) and
    X-linked or mitochondrial, not analysed by eQTLGen (hatched). Relabelled 2026-10-11:
    the hatched group in the primary set is 27 X-linked and 2 mitochondrial genes.
  * Only EXAMPLE-ELIGIBLE targets are named: strict-verified approved CV drug in ChEMBL
    and the same phase in sensitivity set 2. No target that enters through a non-cardiac
    addition (GHR, GLA, PDE4C ...) is drawn in Figure 3.
  * Drug-verified rows: verification completed 2026-10-11 (drug_links/REVIEW_DECISION.md).

Colour: three categorical slots validated for colour-vision deficiency
(blue / orange / aqua; all-pairs CVD dE >= 9.2). Status is never carried by colour
alone: every segment and row is labelled, and the X-linked sub-type uses texture.
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
AOUT = os.path.join(MP, "analysis", "out")
K = json.load(open(os.path.join(AOUT, "key_numbers.json")))
PRIM = list(csv.DictReader(open(os.path.join(AOUT, "per_target_primary.csv"),
                                encoding="utf-8")))

C_INS, C_NOD, C_NOQ, C_UNR = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b2"
C_NOD_X = "#f6b79d"                      # same hue as C_NOD, lighter; always hatched
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2,
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "text.color": INK, "axes.labelcolor": INK,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.facecolor": SURF,
    "savefig.facecolor": SURF, "hatch.linewidth": 0.6, "hatch.color": "#a8502b"})


def save(fig, stem):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(HERE, "%s.%s" % (stem, ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s.{pdf,svg,png}" % stem)


def fmt(n):
    return "{:,}".format(n)


# ====================================================================== Figure 1
def fig1():
    f, s, o = K["flow"], K["sets"]["PRIMARY"], K["overall"]["PRIMARY"]
    fig, ax = plt.subplots(figsize=(7.2, 7.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, lines):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    fc=SURF, ec=INK2, lw=0.8))
        n = len(lines)
        for i, t in enumerate(lines):
            ax.text(x + w / 2, y + h - (i + 0.5) * h / n, t, ha="center", va="center",
                    fontsize=8.5 if i == 0 else 7.6,
                    fontweight="bold" if i == 0 else "normal", color=INK if i == 0 else INK2)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9, shrinkA=0, shrinkB=0))

    cx, w = 18, 46
    box(cx, 90, w, 8, ["%s rows" % fmt(f["raw_rows"]),
                       "Minikel et al. target-indication table"])
    box(cx, 76, w, 8, ["%s target-indication pairs" % fmt(f["distinct_ti_uid"]),
                       "one row per pair identifier"])
    box(cx, 62, w, 8, ["%s pairs" % fmt(f["frozen_pairs"]),
                       "%s targets, %s indications" % (fmt(f["targets"]),
                                                       fmt(f["indications"]))])
    box(cx, 46, w, 10, ["%s cardiovascular pairs" % fmt(s["pairs"]),
                        "%d indications tagged cardiovascular or lipid" % s["terms"],
                        "%d targets" % s["targets"]])
    for y1, y2 in ((90, 84), (76, 70), (62, 56)):
        arrow(cx + w / 2, y1, cx + w / 2, y2)
    ex = [("%d duplicate rows removed" % (f["raw_rows"] - f["distinct_ti_uid"]), 87),
          ("4 programmes removed:\nconflicting phase labels", 73),
          ("%s pairs in %d non-cardiovascular\nindications set aside"
           % (fmt(f["frozen_pairs"] - s["pairs"]), f["indications"] - s["terms"]), 59)]
    for t, y in ex:
        ax.plot([cx + w / 2, 70], [y, y], color=INK2, lw=0.7)
        ax.text(71.5, y, t, ha="left", va="center", fontsize=7.2, color=INK2)
    ax.text(cx + w / 2, 41.2, "Is there a whole-blood cis-eQTL instrument?\n"
            "eQTLGen, p < 5×10⁻⁸, gene body ± 100 kb", ha="center", va="top",
            fontsize=7.6, color=INK2)
    arrow(cx + w / 2, 46, cx + w / 2, 42)
    ys, h = 15, 14
    spec = [(2, "Instrumented", o["instrumented"], C_INS, "white", ""),
            (27, "No qualifying\ninstrument", o["no_qualifying_instrument"], C_NOQ, INK, ""),
            (52, "No eQTL data", o["no_eqtl_data"], C_NOD, "white",
             "%d autosomal, not tested\n%d X-linked, not analysed\n%d mitochondrial, not analysed"
             % (o["no_eqtl_data_autosomal"], o["no_eqtl_data_chrX"],
                o["no_eqtl_data_chrMT"])),
            (77, "Unresolved\ngene record", o["unresolved"], C_UNR, INK, "")]
    for x, lab, n, col, tc, sub in spec:
        ax.add_patch(FancyBboxPatch((x, ys), 21, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    fc=col, ec="none"))
        ax.text(x + 10.5, ys + h - 2.9, lab, ha="center", va="center", fontsize=8.2,
                fontweight="bold", color=tc, linespacing=1.05)
        ax.text(x + 10.5, ys + (8.0 if sub else 4.2), "%d  (%.1f%%)" % (n, 100.0 * n / o["n"]),
                ha="center", va="center", fontsize=8.2, color=tc)
        if sub:
            ax.text(x + 10.5, ys + 3.4, sub, ha="center", va="center", fontsize=5.6, color=tc,
                    linespacing=1.15)
        ax.plot([x + 10.5, x + 10.5], [ys + h, 33], color=INK2, lw=0.7)
    ax.plot([12.5, 87.5], [33, 33], color=INK2, lw=0.7)
    ax.plot([cx + w / 2, cx + w / 2], [33, 35], color=INK2, lw=0.7)
    ax.plot([27, 73], [11.5, 11.5], color=INK, lw=0.9)
    ax.plot([27, 27], [11.5, 13], color=INK, lw=0.9)
    ax.plot([73, 73], [11.5, 13], color=INK, lw=0.9)
    ax.text(50, 8.2, "Uninstrumentable: lower bound %d (%.1f%%, no eQTL data), "
            "upper bound %d (%.1f%%)" % (o["lower_bound"], o["lower_pct"], o["upper_bound"],
                                         o["upper_pct"]),
            ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(50, 4.6, "Autosomal-only lower bound: %d of %d targets not on X, Y or MT (%.1f%%)"
            % (o["autosomal_lower_bound"], o["n_not_on_x_y_mt"], o["autosomal_lower_pct"]),
            ha="center", va="center", fontsize=7.6, color=INK2)
    save(fig, "figure1_flow")


# ====================================================================== Figure 2
def fig2():
    rows = [(p, K["by_phase"]["PRIMARY"][p]) for p in
            ("Preclinical", "Phase I", "Phase II", "Phase III", "Launched")]
    rows += [("All targets", K["overall"]["PRIMARY"]),
             ("Launched, drug ever approved\nby FDA or EMA",
              K["launched"]["PRIMARY"]["approved"]["tally"]),
             ("Launched, approved CV drug,\npermissive ChEMBL rule",
              K["launched"]["PRIMARY"]["permissive"]["tally"])]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ypos = [8.5, 7.5, 6.5, 5.5, 4.5, 3.1, 1.7, 0.6]
    for (lab, t), y in zip(rows, ypos):
        n = float(t["n"])
        parts = [(t["no_eqtl_data_autosomal"], C_NOD, None),
                 (t["no_eqtl_data_x_linked"], C_NOD_X, "/////"),
                 (t["no_qualifying_instrument"], C_NOQ, None),
                 (t["instrumented"], C_INS, None), (t["unresolved"], C_UNR, None)]
        left = 0.0
        for cnt, col, hatch in parts:
            v = 100.0 * cnt / n
            if v > 0:
                ax.barh(y, v, left=left, height=0.62, color=col, edgecolor=SURF,
                        linewidth=1.4, hatch=hatch)
                if hatch:
                    ax.barh(y, v, left=left, height=0.62, color="none", edgecolor=SURF,
                            linewidth=1.4)
            left += v
        lo = 100.0 * t["no_eqtl_data"] / n
        nq = 100.0 * t["no_qualifying_instrument"] / n
        ins = 100.0 * t["instrumented"] / n
        ax.text(100.0 * t["no_eqtl_data_autosomal"] / n / 2, y,   # value = both sub-types
                "%d\n%.0f%%" % (t["no_eqtl_data"], lo), ha="center", va="center",
                fontsize=6.9, color="white", linespacing=0.95)
        if nq >= 5.9:
            ax.text(lo + nq / 2, y, "%d\n%.0f%%" % (t["no_qualifying_instrument"], nq),
                    ha="center", va="center", fontsize=6.9, color=INK, linespacing=0.95)
        elif nq >= 2.5:
            ax.text(lo + nq / 2, y, "%d" % t["no_qualifying_instrument"], ha="center",
                    va="center", fontsize=6.9, color=INK)
        ax.text(lo + nq + ins / 2, y, "%d\n%.0f%%" % (t["instrumented"], ins), ha="center",
                va="center", fontsize=6.9, color="white", linespacing=0.95)
        ax.plot([lo + nq, lo + nq], [y - 0.42, y + 0.42], color=INK, lw=1.1)
        ax.text(101.5, y, "n = %d" % t["n"], ha="left", va="center", fontsize=7.4, color=INK2)
    ax.set_yticks(ypos)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.8, color=INK)
    ax.axhline(3.8, color=GRID, lw=0.8)
    ax.axhline(2.4, color=GRID, lw=0.8)
    ax.set_xlim(0, 100); ax.set_ylim(0.05, 9.05)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("Share of cardiovascular drug targets (%)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    hnd = [Patch(fc=C_NOD), Patch(fc=C_NOD_X, hatch="/////", ec="#a8502b", lw=0),
           Patch(fc=C_NOQ), Patch(fc=C_INS), Patch(fc=C_UNR),
           plt.Line2D([0], [0], color=INK, lw=1.1)]
    ax.legend(hnd, ["No eQTL data: autosomal, not tested", "No eQTL data: X-linked or mitochondrial",
                    "No qualifying instrument", "Instrumented", "Unresolved gene record",
                    "Upper bound"],
              loc="lower center", bbox_to_anchor=(0.42, 1.0), ncol=3, frameon=False,
              fontsize=7.1, handlelength=1.3, columnspacing=1.4, handletextpad=0.5)
    save(fig, "figure2_status_by_phase")


# ====================================================================== Figure 3
def tname(t):
    return t.replace("_", " ").replace(" - ", " ").replace("  ", " ").lower()[:32]


def fig3():
    def pick(phase_label):
        # IL1B is not drawn (owner, 2026-10-11): its approved drug, rilonacept, binds both
        # IL-1 alpha and IL-1 beta and its approved indication is pericarditis, while the
        # source label is atherosclerosis. It stays in Table 2 with a footnote.
        sel = [r for r in PRIM if r["phase_PRIMARY"] == phase_label
               and r["example_eligible"] == "True"
               and r["gene"] not in ("IL1B", "PLG")   # PLG: no FDA or EMA approval record found
               and r["status"] in ("NO_EQTL_DATA", "NO_QUALIFYING_INSTRUMENT")]
        sel.sort(key=lambda r: (r["status"] != "NO_EQTL_DATA",
                                -(float(r["gtex_max_tpm"]) if r["gtex_max_tpm"] else -1),
                                r["gene"]))
        return sel

    def panel(ax, sel, title, nrows):
        lo = 0.008
        clip = lambda v: max(v, lo)
        for i, r in enumerate(sel):
            y = nrows - 1 - i
            col = C_NOD if r["status"] == "NO_EQTL_DATA" else C_NOQ
            xl = " (X)" if r["chromosome"] == "X" else ""
            ax.text(-0.02, y, r["gene"] + xl, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=7.2, fontweight="bold", color=INK)
            if r["gtex_blood_tpm"] == "":
                ax.text(lo * 1.3, y, "no expression record", ha="left", va="center",
                        fontsize=6.2, color=INK2, style="italic")
                continue
            b, m = float(r["gtex_blood_tpm"]), float(r["gtex_max_tpm"])
            ax.plot([clip(b), clip(m)], [y, y], color="#c9c8c2", lw=1.3, zorder=1,
                    solid_capstyle="round")
            ax.scatter([clip(m)], [y], s=34, facecolor=SURF, edgecolor=INK2, linewidth=0.9,
                       zorder=3)
            ax.scatter([clip(b)], [y], s=17, color=col, edgecolor=SURF, linewidth=0.5,
                       zorder=4, marker="o" if r["status"] == "NO_EQTL_DATA" else "D")
            ax.text(clip(m) * 1.6, y, tname(r["gtex_max_tissue"]), ha="left", va="center",
                    fontsize=6.2, color=INK2)
        ax.axvline(1, color=INK2, lw=0.6, ls=(0, (3, 3)), zorder=0)
        ax.set_xscale("log"); ax.set_xlim(lo * 0.8, 2.5e6)
        ax.set_ylim(-0.8, nrows - 0.2)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xticks([0.01, 1, 100, 10000])
        ax.set_xticklabels(["≤0.01", "1", "100", "10,000"])
        ax.minorticks_off()
        ax.set_xlabel("Median expression (TPM, GTEx v10)")
        ax.set_title(title % len(sel), fontsize=8.4, fontweight="bold", loc="left", pad=6)
        ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)

    sa, sb = pick("Launched"), pick("Phase III")
    nrows = max(len(sa), len(sb), 1)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 0.24 * nrows + 1.5),
                                 gridspec_kw=dict(wspace=0.5))
    panel(a1, sa, "a  Launched (n = %d)", nrows)
    panel(a2, sb, "b  Phase III in source data,\n    approved drug in ChEMBL (n = %d)", nrows)
    hnd = [plt.Line2D([0], [0], marker="o", ls="", color=C_NOD, markersize=5.5),
           plt.Line2D([0], [0], marker="D", ls="", color=C_NOQ, markersize=5),
           plt.Line2D([0], [0], marker="o", ls="", markerfacecolor=SURF, markeredgecolor=INK2,
                      markersize=5)]
    fig.legend(hnd, ["Whole blood, target has no eQTL data",
                     "Whole blood, target has no qualifying instrument",
                     "Tissue with the highest expression (named)"],
               loc="lower center", bbox_to_anchor=(0.5, 0.97), ncol=3, frameon=False,
               fontsize=7.2, handletextpad=0.3, columnspacing=1.4)
    save(fig, "figureS2_expression_launched_phase3")
    print("   fig3 names: launched %s ; phase III %s"
          % ([r["gene"] for r in sa], [r["gene"] for r in sb]))


# ====================================================================== Figure S1
def figS1():
    p = r"D:\causal2drug\deliverables\d1h\part3_headtohead.tsv"
    rows = [r for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")
            if r["PP4_eqtlgen"] and r["PP4_gtexblood"]]
    x = [float(r["PP4_gtexblood"]) for r in rows]
    y = [float(r["PP4_eqtlgen"]) for r in rows]
    hi = sum(1 for a, b in zip(x, y) if b > a)
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    lo_, hi_ = 3e-5, 1.0
    ax.plot([lo_, hi_], [lo_, hi_], color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.axvspan(0.8, hi_, color=GRID, zorder=0); ax.axhspan(0.8, hi_, color=GRID, zorder=0)
    ax.scatter(x, y, s=20, color=C_INS, edgecolor=SURF, linewidth=0.7, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_aspect("equal")
    ax.set_xlabel("PP4, GTEx whole blood exposure (n = 853)")
    ax.set_ylabel("PP4, eQTLGen exposure (n = 31,684)")
    ax.text(0.04, 0.95, "higher with eQTLGen: %d loci" % hi, transform=ax.transAxes,
            ha="left", va="top", fontsize=7.4)
    ax.text(0.96, 0.05, "lower with eQTLGen: %d loci" % (len(rows) - hi),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.4)
    ax.text(0.9, 0.93, "PP4 > 0.8", transform=ax.transAxes, ha="right", va="top",
            fontsize=6.6, color=INK2)
    ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
    save(fig, "figureS1_classC_negative_result")
    print("   figS1: %d paired loci, %d higher, %d lower, max PP4 x %.3f y %.3f"
          % (len(rows), hi, len(rows) - hi, max(x), max(y)))


if __name__ == "__main__":
    want = sys.argv[1:] or ["fig1", "fig2", "fig3", "figS1"]
    for w in want:
        {"fig1": fig1, "fig2": fig2, "fig3": fig3, "figS1": figS1}[w]()
