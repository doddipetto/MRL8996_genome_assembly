#!/usr/bin/env python3
"""
02_venn_effectorome.py
======================
Three-set Venn diagram of the MRL8996 / Fol4287 / Fo47 pan-effectorome.

Input : Revision_2026/tables/venn_counts.tsv   (written by 01_pan_effectorome_cdhit.py)
Output: Revision_2026/figures/Venn_Effectorome_MRL8996.pdf
        Revision_2026/figures/Venn_Effectorome_MRL8996.png   (600 dpi)

Counts are CD-HIT clusters (homology groups), not raw proteins: a cluster is
scored as present in a strain if at least one of that strain's effectors falls
into it. The protein counts per strain are given in the caption/table.

No third-party Venn library is required - the three circles are drawn directly
with matplotlib patches so the script runs anywhere matplotlib is available.

Usage: python3 02_venn_effectorome.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

REV = Path(__file__).resolve().parent.parent
COUNTS = REV / "tables" / "venn_counts.tsv"
FIGDIR = REV / "figures"

# MRL8996 is the focal strain (saturated blue); comparators are desaturated.
COLORS = {
    "MRL8996": "#1f5c99",
    "Fol4287": "#c98b3a",
    "Fo47": "#5f9e6e",
}
LABELS = {
    "MRL8996": "MRL8996\n(clinical, this study)",
    "Fol4287": "Fol4287\n(tomato pathogen)",
    "Fo47": "Fo47\n(endophyte)",
}

# Circle geometry: classic three-set arrangement.
R = 1.0
CENTRES = {
    "MRL8996": (0.0, 0.62),
    "Fol4287": (0.60, -0.40),
    "Fo47": (-0.60, -0.40),
}
# Text anchor for each of the seven regions.
REGION_XY = {
    "MRL8996_unique": (0.0, 1.18),
    "Fol4287_unique": (1.12, -0.80),
    "Fo47_unique": (-1.12, -0.80),
    "MRL8996_Fol4287": (0.60, 0.28),
    "MRL8996_Fo47": (-0.60, 0.28),
    "Fol4287_Fo47": (0.0, -0.78),
    "core_all_three": (0.0, -0.06),
}


def read_counts(path):
    clusters, proteins = {}, {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        strains = [h[2:-9] for h in header[2:]]  # n_<strain>_proteins
        for line in fh:
            f = line.rstrip("\n").split("\t")
            clusters[f[0]] = int(f[1])
            proteins[f[0]] = dict(zip(strains, (int(x) for x in f[2:])))
    return clusters, proteins


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    clusters, proteins = read_counts(COUNTS)

    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig, ax = plt.subplots(figsize=(5.2, 5.0))

    for strain, (cx, cy) in CENTRES.items():
        ax.add_patch(Circle((cx, cy), R, facecolor=COLORS[strain], alpha=0.28,
                            edgecolor=COLORS[strain],
                            linewidth=2.2 if strain == "MRL8996" else 1.2, zorder=1))

    # Region counts. The MRL8996-unique region carries the headline number.
    for region, (x, y) in REGION_XY.items():
        n = clusters[region]
        focal = region == "MRL8996_unique"
        ax.text(x, y, str(n), ha="center", va="center", zorder=3,
                fontsize=12 if focal else 10,
                fontweight="bold" if focal else "normal",
                color="#0d3b66" if focal else "0.15")

    # Set labels, placed outside the circles.
    ax.text(0.0, 1.92, LABELS["MRL8996"], ha="center", va="bottom",
            color=COLORS["MRL8996"], fontweight="bold", fontsize=9)
    ax.text(1.52, -1.34, LABELS["Fol4287"], ha="center", va="top",
            color=COLORS["Fol4287"], fontsize=9)
    ax.text(-1.52, -1.34, LABELS["Fo47"], ha="center", va="top",
            color=COLORS["Fo47"], fontsize=9)

    tot = clusters["TOTAL"]
    ax.set_title("135 of 513 MRL8996 effectors form homology groups\n"
                 "absent from both reference strains",
                 fontsize=10, loc="center", pad=14)
    ax.text(0.0, -2.42,
            "CD-HIT clusters (\u2265 50%% identity); %d clusters from "
            "%d + %d + %d effectors" % (tot, proteins["TOTAL"]["MRL8996"],
                                        proteins["TOTAL"]["Fol4287"],
                                        proteins["TOTAL"]["Fo47"]),
            ha="center", va="center", fontsize=7.5, color="0.35")

    ax.set_xlim(-2.35, 2.35)
    ax.set_ylim(-2.6, 2.5)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()

    fig.savefig(FIGDIR / "Venn_Effectorome_MRL8996.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Venn_Effectorome_MRL8996.png", dpi=600, bbox_inches="tight")
    print("wrote %s(.pdf/.png)" % (FIGDIR / "Venn_Effectorome_MRL8996"))
    for region in ("core_all_three", "MRL8996_unique"):
        print("  %-16s %d clusters" % (region, clusters[region]))


if __name__ == "__main__":
    main()
