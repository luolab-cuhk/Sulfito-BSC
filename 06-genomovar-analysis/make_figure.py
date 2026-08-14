#!/usr/bin/env python3
"""
Publication figure for the genomovar analysis (reviewer response).

Panel a  Distribution of within-clade pairwise ANI, restricted to the four
         focal clades (C1, C2, C9, C12).
Panel b  Number of genomovars recovered as a function of ANI threshold.
Panel c  Partitioning of within-clade polymorphism relative to genomovars.

Panel a is shown to define the range over which genomovars were delineated,
not as a test of the reported 99.2-99.8% gap: SI Text 1.3 notes that with at
most 20 genomes per clade, the shape of this distribution is not resolvable
at the density at which the gap was originally described.

Outputs: Fig_genomovar.pdf / .png / .tif
Run genomovar_analysis.py first (this script re-uses its functions).
"""
import csv
import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import genomovar_analysis as G

OUT = os.path.dirname(os.path.abspath(__file__))
GAP_LO, GAP_HI = G.GAP_LO, G.GAP_HI

plt.rcParams.update({
    "font.family": "sans-serif",
    # Helvetica on macOS is a .ttc collection whose metrics matplotlib cannot
    # read reliably (breaks layout); Arial is a plain .ttf and renders identically
    # for our purposes.
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Keep the set and ordering of clades shown in every panel in one place.
FOCAL_CLADES = ("C1", "C2", "C9", "C12")
CLADE_COLORS = {"C1": "#3B6FB6", "C12": "#D1495B", "C2": "#00857C", "C9": "#E8A33D"}
GAP_COLOR = "#BBBBBB"


def collect():
    """Recompute everything needed for plotting."""
    strains = G.load_strains()
    _, get = G.load_ani()
    all_clades = defaultdict(list)
    for s, c in strains.items():
        all_clades[c].append(s)
    for c in all_clades:
        all_clades[c].sort()
    # Keep all plotting inputs restricted to the four focal clades.
    clades = {c: all_clades[c] for c in FOCAL_CLADES if c in all_clades}
    focal = list(clades)

    within = {c: G.within_pairs(clades[c], get) for c in focal}

    thresholds = np.arange(99.0, 100.01, 0.05)
    ngv = {c: [len(set(G.cluster_at(clades[c], get, t, "average").values()))
               for t in thresholds] for c in focal}
    keep = set(s for c in focal for s in clades[c])
    seqs = G.load_alignment(keep)
    part = {}
    for c in focal:
        members = clades[c]
        arr, mask = G.unambiguous_mask(seqs, members)
        pos, sub = G.polymorphic_columns(arr, mask)
        rowc = {}
        for thr in (99.2, 99.5, 99.8, 99.9):
            gvv = G.cluster_at(members, get, thr, "average")
            r = G.genomovar_partition_test(members, gvv, sub, pos)
            # A genomovar containing one genome is internally fixed by
            # definition, so a partition made of singletons inflates
            # pct_between without saying anything about population structure.
            # Record the largest genomovar size so the panel can report it:
            # the size, not a verdict, is what tells the reader how much the
            # partition actually constrains.
            sizes = Counter(gvv.values())
            max_sz = max(sizes.values())
            rowc[thr] = (r["between"], r["within"], r["n_gv"], max_sz)
        part[c] = rowc
    return clades, focal, thresholds, ngv, part, within


def panel_a(ax, within, clades):
    ax.axvspan(GAP_LO, GAP_HI, color=GAP_COLOR, alpha=0.45, lw=0, zorder=0)
    bins = np.arange(99.0, 100.0001, 0.02)
    bottom = np.zeros(len(bins) - 1)
    n_pairs = n_gap = 0
    for c in FOCAL_CLADES:
        if c not in within:
            continue
        vals = within[c]
        n_pairs += len(vals)
        n_gap += sum(1 for v in vals if GAP_LO <= v <= GAP_HI)
        counts, _ = np.histogram(vals, bins=bins)
        ax.bar(bins[:-1], counts, width=np.diff(bins), align="edge",
               bottom=bottom, color=CLADE_COLORS[c], lw=0,
               label="%s (%d pairs)" % (c, len(vals)), zorder=3)
        bottom += counts
    ax.text(0.03, 0.96, "%d of %d pairs (%.0f%%)\nfall inside the reported\n%.1f-%.1f%% gap"
            % (n_gap, n_pairs, 100.0 * n_gap / n_pairs, GAP_LO, GAP_HI),
            transform=ax.transAxes, ha="left", va="top", fontsize=6)
    ax.set_xlim(99.0, 100.0)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_xlabel("Within-clade pairwise ANI (%)")
    ax.set_ylabel("Number of genome pairs")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=9)
    ax.legend(frameon=False, loc="upper right", handlelength=1.3,
              borderpad=0.2, labelspacing=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def panel_b(ax, thresholds, ngv, clades):
    ax.axvspan(GAP_LO, GAP_HI, color=GAP_COLOR, alpha=0.45, lw=0, zorder=0)
    for c in FOCAL_CLADES:
        if c not in ngv:
            continue
        ax.step(thresholds, ngv[c], where="post", color=CLADE_COLORS[c],
                lw=1.1, label="%s (n=%d)" % (c, len(clades[c])), zorder=3)
    ax.set_xlim(99.0, 100.0)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_xlabel("ANI threshold for genomovar delineation (%)")
    ax.set_ylabel("Number of genomovars")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=9)
    ax.legend(frameon=False, loc="upper left", handlelength=1.3,
              borderpad=0.2, labelspacing=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def panel_c(ax, part):
    clades = [c for c in FOCAL_CLADES if c in part]
    thrs = (99.2, 99.5, 99.8, 99.9)
    width = 0.2
    xs = np.arange(len(clades))
    for k, thr in enumerate(thrs):
        vals, labels = [], []
        for c in clades:
            btw, wth, ngv, max_sz = part[c][thr]
            tot = btw + wth
            vals.append(100.0 * btw / tot if tot else 0.0)
            labels.append(str(ngv))
        off = (k - (len(thrs) - 1) / 2) * width
        # 99.9% lies above the proposed 99.2-99.8% genomovar interval; give it a
        # separate colour so it is not read as part of the criterion being tested.
        inside = thr <= GAP_HI
        base = "#3B6FB6" if inside else "#7B5AA6"
        shade = (0.3 + 0.7 * k / (len(thrs) - 2)) if inside else 0.85
        bars = ax.bar(xs + off, vals, width * 0.92, color=base,
                      alpha=shade, lw=0.3, edgecolor="#20334D",
                      label="%.1f%%" % thr if inside
                            else "%.1f%% (outside)" % thr, zorder=3)
        for b, ngv in zip(bars, labels):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5,
                    str(ngv), ha="center", va="bottom", fontsize=5.2,
                    color="#333333")
    ax.set_xticks(xs)
    ax.set_xticklabels(clades)
    # Headroom so the legend sits in a clear band above the tallest bar (100%).
    ax.set_ylim(0, 152)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("Within-clade polymorphic sites\ndifferentially fixed between genomovars (%)")
    ax.set_xlabel("Clade (numerals above bars: genomovars recovered)")
    ax.set_title("c", loc="left", fontweight="bold", fontsize=9)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, frameon=False, title="ANI threshold",
              loc="upper center", ncol=3, handlelength=1.0, borderpad=0.2,
              columnspacing=1.0, labelspacing=0.3, handletextpad=0.5,
              title_fontsize=6.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    clades, multi, thresholds, ngv, part, within = collect()
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6), layout="constrained")
    panel_a(axes[0], within, clades)
    panel_b(axes[1], thresholds, ngv, clades)
    panel_c(axes[2], part)
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.04, wspace=0.04, hspace=0.04)
    for ext, kw in (("pdf", {}), ("png", dict(dpi=600)), ("tif", dict(dpi=600))):
        p = os.path.join(OUT, "Fig_genomovar.%s" % ext)
        fig.savefig(p, bbox_inches="tight", **kw)
        print("wrote %s" % p)
    plt.close(fig)


if __name__ == "__main__":
    main()
