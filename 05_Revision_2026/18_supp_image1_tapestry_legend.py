#!/usr/bin/env python3
"""
Supplementary Image 1: add the missing read-coverage key to the tapestry
contig plot (reviewer comment R2.21).

The tapestry plot itself is kept exactly as exported (the screenshot in
02_telomeres_centromeres/tapestry_report).  What is added around it:

  * the read-depth key.  Tapestry does not colour the bars by raw depth: it
    assigns each 200-kb window a copy-number level from the Nanopore read
    depth and colours the level on a six-step Greens scale
    (ploidy_colour(): d3.interpolateGreens(0.01 + level/5 * 0.99), levels
    0..5 = none, 0.5x, 1x, 1.5x, 2x, >=2.5x of the genome-median depth).
    Tapestry's own labels are haploid/diploid/triploid/tetraploid/repeat;
    for a haploid fungus they are given here as multiples of the median.
    The swatch colours are the ones tapestry drew (sampled from the plot);
    the two levels absent from the plot (0 and >=2.5x) use the same d3
    scale.  Median depth and window size are read from the report HTML.
  * the telomere key (red bar at a contig end = TTAGGG repeats found;
    opacity proportional to the number of repeats, tapestry convention).
  * chromosome names next to the scaffold names (the assembly names its
    16 pseudomolecules by size rank; PGA_scaffold_3 and _4 swap).
  * the x-axis title, which the screenshot cut off.

Output: Revision_2026/figures/Supplementary_Image_1_tapestry.{png,pdf}
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

HIC = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_HiC_data")
REP = HIC / "02_telomeres_centromeres/tapestry_report"
SHOT = REP / "Screenshot from 2026-06-19 14-41-51.png"
HTML = REP / "tapestry_report.tapestry_report.html"
REV = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026")
OUT = REV / "figures" / "Supplementary_Image_1_tapestry"

# tapestry's Greens at levels 0..5; 1..4 sampled from the exported plot,
# 0 and 5 from d3.interpolateGreens(0.01) and (1.0)
LEVELS = [("0",        (247, 252, 245)),
          ("0.5x",     (209, 237, 202)),
          ("1x",       (149, 211, 146)),
          ("1.5x",     (76, 175, 97)),
          ("2x",       (21, 126, 59)),
          ("\u22652.5x", (0, 68, 27))]
TELO_RED = (220, 20, 60)


def main():
    html = HTML.read_text()
    median = float(re.search(r"median_depth\s*=\s*([\d.]+)", html).group(1))
    win = int(re.search(r'"Window size", "value": (\d+)', html).group(1))

    im = np.array(Image.open(SHOT).convert("RGB"))
    H, W = im.shape[:2]

    # row centres = the red telomere bars; bar x-extent = the green run
    r, g, b = [im[..., i].astype(int) for i in range(3)]
    red = (r > 200) & (g < 80) & (b < 80)
    rows = np.nonzero(red.sum(1) > 2)[0]
    groups = np.split(rows, np.nonzero(np.diff(rows) > 3)[0] + 1)
    ys = [(x[0] + x[-1]) / 2 for x in groups]
    assert len(ys) == 16, len(ys)
    green = (g > r + 20) & (g > b + 20)
    x_right = np.nonzero(green.sum(0) > 0)[0].max()

    # scaffold -> chromosome: pseudomolecules are named by size rank
    names = re.findall(r"PGA_scaffold_(\d+)__\d+_contigs__length_(\d+)", html)
    lengths = {int(s): int(l) for s, l in set(names)}
    ranked = sorted(lengths, key=lambda s: -lengths[s])
    chrom_of = {s: "Chr%02d" % (i + 1) for i, s in enumerate(ranked)}
    # the plot lists scaffolds in the report order (by length, as drawn)
    row_scafs = sorted(lengths, key=lambda s: -lengths[s])

    # canvas: plot at native pixels, extra room right (Chr names) + below (key)
    PAD_R, PAD_B = 125, 165
    dpi = 100
    fig = plt.figure(figsize=((W + PAD_R) / dpi, (H + PAD_B) / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W + PAD_R)
    ax.set_ylim(H + PAD_B, 0)
    ax.axis("off")
    ax.imshow(im, extent=(0, W, H, 0), interpolation="none")

    for y, s in zip(ys, row_scafs):
        ax.text(W + 12, y, chrom_of[s], fontsize=10, va="center", ha="left",
                color="black", family="DejaVu Sans")
    ax.text(W + 12, ys[0] - 22, "chromosome", fontsize=9, ha="left", va="bottom",
            color="#444444", style="italic")
    # the screenshot's bottom edge carries the clipped top of tapestry's own
    # axis title; mask it before writing the full title below
    ax.add_patch(Rectangle((0, H - 14), W, 14, facecolor="white", edgecolor="none", zorder=3))

    # x-axis title (cut off in the screenshot)
    ax.text(380 + (x_right - 380) / 2, H + 14, "Position along scaffold (Mb)",
            fontsize=11, ha="center", va="top")

    # ---- key -------------------------------------------------------------
    ky = H + 62
    ax.text(20, ky - 6, "Read depth (Nanopore reads mapped to the assembly), "
            "%d-kb windows, as multiple of the genome-median depth (%.1f\u00d7):"
            % (win // 1000, median), fontsize=10, ha="left", va="bottom")
    x = 20
    for lab, col in LEVELS:
        ax.add_patch(Rectangle((x, ky), 34, 16, facecolor=np.array(col) / 255,
                               edgecolor="#777777", linewidth=0.6))
        ax.text(x + 40, ky + 8, lab, fontsize=10, va="center", ha="left")
        x += 40 + 12 * len(lab) + 20
    # telomere key
    ty = ky + 40
    ax.add_patch(Rectangle((20, ty - 2), 4, 20, facecolor=np.array(TELO_RED) / 255,
                           edgecolor="none"))
    ax.text(34, ty + 8, "telomere: TTAGGG repeats at the contig end "
            "(opacity proportional to the number of repeats found)",
            fontsize=10, va="center", ha="left")
    ax.text(20, ty + 40, "1\u00d7 = single-copy sequence. Windows above 1\u00d7 carry more reads than one genomic "
            "copy provides: the rDNA array (Chr05 end),\na mitochondrial insertion (Chr05, 4.55 Mb) "
            "and a ~0.5-Mb segment at the start of Chr13 (see Methods).",
            fontsize=8.5, va="center", ha="left", color="#333333", linespacing=1.5)

    fig.savefig(str(OUT) + ".png", dpi=200)
    fig.savefig(str(OUT) + ".pdf")
    print("median depth %.2f, window %d, rows %d -> %s.{png,pdf}" % (median, win, len(ys), OUT))


if __name__ == "__main__":
    main()
