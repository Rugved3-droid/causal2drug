#!/usr/bin/env python
"""Publication figures for the cardiovascular untestability atlas.

One command regenerates all four from frozen outputs:
    python src/make_figures.py

Vector PDF + SVG (plus PNG previews), Okabe-Ito colourblind-safe palette, no chartjunk.
Axis labels and captions are written for a cardiology audience: "measurable in blood"
rather than "cis-eQTL instrument available".
"""
import collections
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

CV = r"D:\causal2drug\deliverables\cv"
FIG = r"D:\causal2drug\deliverables\figures"

BLUE, ORANGE, GREEN, VERM = "#0072B2", "#E69F00", "#009E73", "#D55E00"
PURPLE, GREY = "#CC79A7", "#7F7F7F"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})

NL = chr(10)


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(FIG, "%s.%s" % (name, ext)),
                    bbox_inches="tight", dpi=(600 if ext != "png" else 150))
    plt.close(fig)
    print("  wrote %s.pdf / .svg / .png" % name)


def caption(fig, text, y=-0.02):
    fig.text(0.02, y, text, fontsize=7.4, color="#555555", va="top", ha="left")


def load():
    return list(csv.DictReader(open(CV + r"\cardiovascular_atlas.csv", encoding="utf-8")))


# ---------------------------------------------------------------- Figure 1
def fig1(atlas):
    st = collections.Counter(r["state"] for r in atlas)
    lau = collections.Counter(r["state"] for r in atlas if r["max_phase_num"] == "5")
    n = len(atlas)
    ca = st["ABSENT"] + st["TESTED_NULL"]
    ca_lau = lau["ABSENT"] + lau["TESTED_NULL"]

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(4, 100)
    ax.axis("off")

    def box(x, y, w, h, col, rows):
        """rows: list of (relative_y, text, fontsize, style) drawn inside the box."""
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.5,rounding_size=1.5",
                                    facecolor=col, edgecolor="none"))
        for ry, txt, fs, sty in rows:
            ax.text(x + w / 2, y + h * ry, txt, ha="center", va="center", color="white",
                    fontsize=fs, fontweight=("bold" if sty == "b" else "normal"),
                    style=("italic" if sty == "i" else "normal"))

    box(22, 84, 56, 15, GREY,
        [(0.76, "630 cardiovascular drug targets", 11.5, "b"),
         (0.45, "every target with a cardiovascular indication", 8.8, "n"),
         (0.16, "%d have an approved medicine" % sum(lau.values()), 8.0, "i")])

    ax.add_patch(FancyArrowPatch((50, 84), (50, 76), arrowstyle="-|>",
                                 mutation_scale=14, lw=1.2, color="#444444"))
    ax.text(53, 80,
            "Can the target's activity be measured genetically" + NL
            + "from a routine blood sample?",
            fontsize=8.6, va="center", ha="left", color="#333333", linespacing=1.5)

    def state_box(x, col, title, key, desc):
        box(x, 44, 30, 28, col,
            [(0.83, title, 11, "b"),
             (0.62, "%d targets  (%.1f%%)" % (st[key], 100 * st[key] / n), 9.6, "n"),
             (0.33, desc, 8.8, "n"),
             (0.10, "%d have an approved medicine" % lau[key], 8.0, "i")])

    state_box(1, VERM, "ABSENT", "ABSENT", "not detectable in blood")
    state_box(35, ORANGE, "TESTED, NO SIGNAL", "TESTED_NULL", "no usable genetic marker")
    state_box(69, BLUE, "MEASURABLE", "INSTRUMENTED", "genetic marker available")
    for xt in (16, 50, 84):
        ax.add_patch(FancyArrowPatch((50, 76), (xt, 72.5), arrowstyle="-|>",
                                     mutation_scale=12, lw=1.0, color="#666666"))

    ax.add_patch(FancyBboxPatch((1, 6), 64, 27,
                                boxstyle="round,pad=0.5,rounding_size=1.5",
                                facecolor="#F7F2F0", edgecolor=VERM, lw=1.6))
    ax.text(33, 28.0, "CANNOT BE TESTED", ha="center", fontsize=12.5,
            fontweight="bold", color=VERM)
    ax.text(33, 21.5, "%d targets  (%.1f%% of all cardiovascular targets)"
            % (ca, 100 * ca / n), ha="center", fontsize=10.2, color=VERM)
    ax.text(33, 14.0, "including %d with an approved medicine, among them" % ca_lau,
            ha="center", fontsize=8.6, color="#333333")
    ax.text(33, 9.5, "PCSK9, NPC1L1 (ezetimibe), APOC3 and ANGPTL3",
            ha="center", fontsize=8.6, color="#333333")
    ax.text(84, 22.5, "%d targets" % st["INSTRUMENTED"], ha="center", fontsize=11.5,
            color=BLUE, fontweight="bold")
    ax.text(84, 17.5, "can be tested", ha="center", fontsize=11.5,
            color=BLUE, fontweight="bold")

    caption(fig,
            "Testable means a common genetic variant exists that measurably changes the "
            "target's activity in blood, so the target's effect on disease can be "
            "estimated" + NL
            + "from human genetics. A further %d targets (%.1f%%) could not be matched to "
            "a gene record and are excluded from the three groups above." % (
                st["UNMAPPED"], 100 * st["UNMAPPED"] / n) + NL
            + "Reproduce: python src/cv_atlas.py ; python src/make_figures.py")
    save(fig, "figure1_untestability_cascade")
    return st, lau


# ---------------------------------------------------------------- Figure 2
def fig2(atlas):
    n = len(atlas)
    keys = [("instrumented_blood", "Blood" + NL + "(convenient)", BLUE),
            ("instrumented_cardiac_union",
             "Any of the four" + NL + "heart / vessel /" + NL + "liver tissues", GREEN),
            ("instrumented_artery_aorta", "Aorta", ORANGE),
            ("instrumented_heart_LV", "Heart muscle" + NL + "(left ventricle)", VERM),
            ("instrumented_heart_atrial", "Heart" + NL + "(atrium)", VERM),
            ("instrumented_liver", "Liver", VERM)]
    vals = [100 * sum(1 for r in atlas if r[k] == "True") / n for k, _, _ in keys]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.bar(range(len(vals)), vals, color=[c for _, _, c in keys], width=.64)
    for i, v in enumerate(vals):
        ax.text(i, v + 1.6, "%.1f%%" % v, ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([l for _, l, _ in keys], fontsize=8.5)
    ax.set_ylabel("Cardiovascular targets that can be tested (%)", fontsize=9.2)
    ax.set_ylim(0, 78)
    ax.axhline(vals[0], color=BLUE, lw=.8, ls=":", alpha=.55)
    ax.axhline(vals[3], color=VERM, lw=.8, ls=":", alpha=.55)
    ratio = vals[0] / vals[3]
    ax.annotate("", xy=(3.5, vals[3]), xytext=(3.5, vals[0]),
                arrowprops=dict(arrowstyle="<->", color="#333333", lw=1.3))
    ax.text(3.68, (vals[0] + vals[3]) / 2,
            "%.1f-fold fewer" % ratio + NL + "targets testable",
            fontsize=10.5, fontweight="bold", color="#333333", va="center",
            linespacing=1.5)
    ax.set_title("The tissue that matters clinically is the one we cannot test in",
                 fontsize=11, pad=12, loc="left")
    caption(fig,
            "Blood is used because samples are abundant, not because it is where "
            "cardiovascular drugs act. Heart, vessel and liver samples are scarce, so far "
            "fewer" + NL
            + "targets reach the threshold for a usable genetic test. "
            "Reproduce: python src/cv_atlas.py ; python src/make_figures.py")
    save(fig, "figure2_tissue_yield")
    return vals


# ---------------------------------------------------------------- Figure 3
def fig3(atlas, yields):
    ca = [r for r in atlas if r["untestability_class"] == "CLASS_A"]
    need = collections.Counter(r["required_modality"].split(" (")[0] for r in ca)
    tot = len(ca)
    _blood, _union, aorta, hlv, _hat, liver = yields
    rows = [("Heart / vessel tissue", need["cardiac eQTL"], hlv, VERM),
            ("Liver tissue", need["liver eQTL"], liver, ORANGE),
            ("Vessel wall (aorta)", need["vascular eQTL"], aorta, GREEN),
            ("Circulating protein", need["pQTL"], None, PURPLE)]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    for i, (lab, nneed, yld, col) in enumerate(rows):
        pct = 100 * nneed / tot
        ax.barh(i - .19, pct, height=.34, color=col)
        ax.text(pct + 1.0, i - .19, "%.1f%% need this  (n=%d)" % (pct, nneed),
                va="center", fontsize=8.6, fontweight="bold")
        if yld is None:
            ax.barh(i + .19, 0.5, height=.34, color="white", edgecolor=col,
                    hatch="////", lw=1.0)
            ax.text(1.8, i + .19, "no accessible resource", va="center",
                    fontsize=8.6, color=GREY, style="italic")
        else:
            ax.barh(i + .19, yld, height=.34, color=col, alpha=.40)
            ax.text(yld + 1.0, i + .19, "only %.1f%% testable this way" % yld,
                    va="center", fontsize=8.6, color="#444444")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.4)
    ax.set_xlabel("Percentage of cardiovascular targets", fontsize=9.2)
    ax.set_xlim(0, 66)
    ax.invert_yaxis()
    ax.set_title("What untestable targets need, and what that measurement delivers",
                 fontsize=11, pad=12, loc="left")
    caption(fig,
            "Dark bar: share of the %d untestable cardiovascular targets that would "
            "require this measurement.   Pale bar: share of all cardiovascular targets "
            "that" % tot + NL
            + "measurement can actually test today. Nearly half of untestable targets need "
            "heart or vessel tissue, the measurement that can test only %.1f%% of all "
            "cardiovascular targets." % hlv + NL
            + "Circulating-protein data would cover a further %.1f%%, but no resource is "
            "accessible: UK Biobank-PPP requires credentials and deCODE could not be "
            "verified." % (100 * need["pQTL"] / tot) + NL
            + "Reproduce: python src/cv_atlas.py ; python src/make_figures.py")
    save(fig, "figure3_closed_loop")


# ---------------------------------------------------------------- Figure 4
def fig4():
    rob = json.load(open(CV + r"\cv_robustness.json"))
    gw = json.load(open(r"D:\causal2drug\deliverables\d1e\classA_tests.json"))[
        "absent_vs_instrumented_tissue_restricted"]

    def ntot(tab):
        return sum(tab[0]) + sum(tab[1])

    rows = [("Cardiovascular targets" + NL + "(atlas definition)",
             rob["CV_EXTENDED"]["or_"], rob["CV_EXTENDED"]["or_lo"],
             rob["CV_EXTENDED"]["or_hi"], ntot(rob["CV_EXTENDED"]["or_table"]), VERM),
            ("Cardiovascular targets" + NL + "(strict definition)",
             rob["CV_PRIMARY"]["or_"], rob["CV_PRIMARY"]["or_lo"],
             rob["CV_PRIMARY"]["or_hi"], ntot(rob["CV_PRIMARY"]["or_table"]), ORANGE),
            ("All drug targets", gw["odds_ratio"], gw["ci95"][0], gw["ci95"][1],
             ntot(gw["table"]), BLUE)]
    fig, ax = plt.subplots(figsize=(7.8, 2.9))
    for i, (lab, o, lo, hi, nn, col) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=col, lw=2.4, solid_capstyle="round")
        ax.plot([o], [i], "o", color=col, ms=8.5, zorder=3)
        ax.text(26, i, "%.1f   (95%% CI %.1f-%.1f)" % (o, lo, hi),
                va="center", fontsize=9)
    ax.axvline(1, color="#999999", lw=.9, ls="--")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(["%s%sn = %s targets compared" % (r[0], NL, "{:,}".format(r[4]))
                        for r in rows], fontsize=8.6)
    ax.set_xscale("log")
    ax.set_xlim(0.85, 75)
    ax.set_xticks([1, 2, 5, 10, 20])
    ax.set_xticklabels(["1", "2", "5", "10", "20"])
    ax.set_xlabel("Odds that an untestable target is expressed in only a few tissues"
                  + NL + "(odds ratio vs testable targets, log scale)", fontsize=9)
    ax.invert_yaxis()
    ax.set_title("Targets we cannot test are disproportionately tissue-restricted",
                 fontsize=11, pad=12, loc="left")
    caption(fig,
            "Values above 1 mean untestable targets are more likely to be expressed in "
            "only a few tissues. The attrition is therefore not random: it removes a "
            "specific kind" + NL
            + "of target, and does so slightly more strongly among cardiovascular targets "
            "than across all drug targets. Counts are the targets carrying the tissue-"
            "specificity" + NL
            + "annotation needed for the comparison, which is fewer than the full sets "
            "(630 cardiovascular, 594 strict)." + NL
            + "Reproduce: python src/a2_robustness.py ; python src/make_figures.py",
            y=-0.17)
    save(fig, "figure4_classA_enrichment_forest")


def main():
    atlas = load()
    print("generating figures from %d atlas rows" % len(atlas))
    fig1(atlas)
    yields = fig2(atlas)
    fig3(atlas, yields)
    fig4()
    print("all figures written to", FIG)


if __name__ == "__main__":
    main()
