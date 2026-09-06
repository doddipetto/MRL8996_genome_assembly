#!/usr/bin/env python3
"""
10_effector_genomic_map.py
==========================
Rebuilds the linear effector map (Supplementary Image 2 in v1.1) as a
main-text figure, at publication resolution and with the information the
reviewers asked for added:

    panel a  chromosome bars shaded by compartment (core Chr01-11 vs
             accessory Chr12-16), one tick per effector coloured by AMAPEC
             prediction, MRL8996-unique effectors marked above the bar
    panel b  effector density per Mb, core vs accessory, with per-chromosome
             points overlaid

Reimplemented in Python rather than extending plot_map_effector_position.R:
the original script reads effector_map_data.tsv (483 effectors with a
structural model) and has no notion of the unique/shared split, whereas this
version uses the full final effectorome with the unique-effector annotation
produced by script 03. Original outputs are untouched.

Inputs:
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv   (script 03)
    MRL8996_HiC_data/08_circos_features/chrom.sizes

Outputs:
    Revision_2026/figures/Figure_effector_genomic_distribution.pdf
    Revision_2026/figures/Figure_effector_genomic_distribution.png  (600 dpi)
    Revision_2026/tables/effector_density_by_compartment.tsv

Usage: python3 10_effector_genomic_map.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
SIZES = HIC / "08_circos_features" / "chrom.sizes"
ANN = REV / "tables" / "Table_MRL8996_unique_effectors.tsv"

ACCESSORY = {"Chr12", "Chr13", "Chr14", "Chr15", "Chr16"}
AM_C = "#b2182b"        # antimicrobial
NONAM_C = "#2166ac"     # non-antimicrobial
UNIQ_C = "#e08214"      # MRL8996-unique marker
CORE_BAR = "#e6e6e6"
ACC_BAR = "#cfd8dc"


def main():
    sizes = pd.read_csv(SIZES, sep="\t", header=None, names=["Chromosome", "Length"])
    sizes = sizes[sizes["Chromosome"].str.match(r"Chr\d+")].sort_values("Chromosome")
    ann = pd.read_csv(ANN, sep="\t")
    ann = ann[ann["Chromosome"].isin(set(sizes["Chromosome"]))].copy()

    order = list(sizes["Chromosome"])
    ypos = {c: i for i, c in enumerate(order)}

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(9.0, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.0, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])

    for _, r in sizes.iterrows():
        y = ypos[r["Chromosome"]]
        col = ACC_BAR if r["Chromosome"] in ACCESSORY else CORE_BAR
        ax.add_patch(plt.Rectangle((0, y - 0.28), r["Length"] / 1e6, 0.56,
                                   facecolor=col, edgecolor="none", zorder=1))

    for _, e in ann.iterrows():
        y = ypos[e["Chromosome"]]
        x = e["Start"] / 1e6
        c = AM_C if str(e["AMAPEC"]).lower().startswith("anti") else NONAM_C
        ax.plot([x, x], [y - 0.26, y + 0.26], color=c, lw=0.7, zorder=2)
        if bool(e["Unique_to_MRL8996"]):
            ax.plot([x], [y - 0.42], marker="v", markersize=2.6, color=UNIQ_C,
                    zorder=3, clip_on=False)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(-0.05, sizes["Length"].max() / 1e6 * 1.02)
    ax.set_ylim(len(order) - 0.4, -0.9)
    ax.set_xlabel("Genomic position (Mb)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("a  Effector positions across the 16 chromosomes",
                 fontsize=9.5, loc="left")

    handles = [plt.Line2D([], [], color=AM_C, lw=1.4, label="antimicrobial (AMAPEC)"),
               plt.Line2D([], [], color=NONAM_C, lw=1.4, label="non-antimicrobial"),
               plt.Line2D([], [], color=UNIQ_C, marker="v", lw=0, markersize=4,
                          label="MRL8996-unique"),
               plt.Rectangle((0, 0), 1, 1, facecolor=ACC_BAR,
                             label="accessory chromosome")]
    ax.legend(handles=handles, frameon=False, fontsize=7, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.13),
              handletextpad=0.5, columnspacing=1.2)

    # ---- panel b: density ----
    ax2 = fig.add_subplot(gs[0, 1])
    sizes["compartment"] = np.where(sizes["Chromosome"].isin(ACCESSORY),
                                    "accessory", "core")
    n = ann.groupby("Chromosome").size()
    sizes["effectors"] = sizes["Chromosome"].map(n).fillna(0)
    sizes["per_Mb"] = sizes["effectors"] / (sizes["Length"] / 1e6)

    dens = []
    for comp in ("core", "accessory"):
        sub = sizes[sizes["compartment"] == comp]
        dens.append({"compartment": comp,
                     "chromosomes": len(sub),
                     "length_Mb": round(sub["Length"].sum() / 1e6, 2),
                     "effectors": int(sub["effectors"].sum()),
                     "effectors_per_Mb": round(sub["effectors"].sum()
                                               / (sub["Length"].sum() / 1e6), 2)})
    dens = pd.DataFrame(dens)
    dens.to_csv(REV / "tables" / "effector_density_by_compartment.tsv",
                sep="\t", index=False)

    xs = [0, 1]
    ax2.bar(xs, dens["effectors_per_Mb"], color=["#8c9196", AM_C], width=0.6,
            zorder=1)
    rng = np.random.default_rng(0)
    for i, comp in enumerate(("core", "accessory")):
        v = sizes.loc[sizes["compartment"] == comp, "per_Mb"].values
        ax2.scatter(np.full(len(v), i) + rng.uniform(-0.16, 0.16, len(v)), v,
                    s=12, facecolor="white", edgecolor="0.25", linewidth=0.7,
                    zorder=3)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(["core\n(Chr01-11)", "accessory\n(Chr12-16)"], fontsize=8)
    ax2.set_ylabel("Effectors per Mb")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_title("b  Compartment density", fontsize=9.5, loc="left")
    top = max(sizes["per_Mb"].max(), dens["effectors_per_Mb"].max()) * 1.32
    ax2.set_ylim(0, top)
    for i, r in dens.iterrows():
        ax2.text(i, top * 0.99, "%d effectors\n%.1f/Mb"
                 % (r["effectors"], r["effectors_per_Mb"]),
                 ha="center", va="top", fontsize=7)

    fig.savefig(REV / "figures" / "Figure_effector_genomic_distribution.pdf",
                bbox_inches="tight")
    fig.savefig(REV / "figures" / "Figure_effector_genomic_distribution.png",
                dpi=600, bbox_inches="tight")

    print("effectors plotted: %d (unique %d)"
          % (len(ann), int(ann["Unique_to_MRL8996"].sum())))
    print(dens.to_string(index=False))


if __name__ == "__main__":
    main()
