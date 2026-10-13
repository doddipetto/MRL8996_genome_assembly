#!/usr/bin/env python3
"""
Add chromosome labels to Figure 1 panel A (the cropped Juicebox contact map).

The Juicebox export carries only an assembly-coordinate scale in Mb; the 16
pseudomolecules are marked by the blue boxes drawn along the diagonal but are
not named. This script names them.

The label positions are not assumed - they are measured from the figure itself:

  1. the page is rasterised at 1 px = 1 pt and the blue box strokes are
     segmented by hue (B > 120, B - R > 60, B - G > 60);
  2. the row-wise blue profile is grouped into runs, giving 17 horizontal box
     boundaries, i.e. the 16 pseudomolecule limits in page coordinates;
  3. those spacings are checked against the cumulative chromosome lengths in
     chrom.sizes. They agree to < 0.3 % of the map width, which confirms both
     that the map is drawn in the released Chr01..Chr16 order (descending
     length) and that the axis is linear in bp;
  4. the same is done column-wise, matching candidate edges to the positions
     predicted from cumulative length (the column profile also picks up faint
     grey gridlines, so it is matched rather than used raw), and a linear
     bp -> pt fit is taken from the matched edges.

Labels are then drawn in a band added above the map, rotated 90 degrees so that
even the 0.87-Mb Chr16 slot can hold a label, with a tick at every boundary.
The font size is chosen for the panel's FINAL printed width (~85 mm in the
two-column Figure 1), not for the 2434-pt export, hence the large nominal value.

Inputs (read-only):
    Revision_2026/figures/Figure_1A_HiC_original_cropped.pdf
    MRL8996_HiC_data/11_Circos_Plot/chrom.sizes

Output:
    Revision_2026/figures/Figure_1A_HiC_labelled.pdf   (vector; this is the
        panel consumed by 11b_reexport_main_figures.py)

Run:  python3 scripts/16_label_hic_panelA.py     (env: python, needs pymupdf)
"""

import numpy as np
import pymupdf
import pandas as pd

# --- paths ---
REV = "/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026/"
SRC = REV + "figures/Figure_1A_HiC_original_cropped.pdf"
OUT = REV + "figures/Figure_1A_HiC_labelled.pdf"
SIZES = ("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_HiC_data/"
         "11_Circos_Plot/chrom.sizes")

# --- label geometry, in points of the exported panel ---
BAND = 230.0      # height of the added label band
FONTSIZE = 58     # ~7 pt once the panel is scaled to 85 mm
TICK = 12.0
SKIP_TOP_PX = 60  # ignore the blue Mb tick labels above the map


def boundaries(page, cum):
    """17 box boundaries in page pt, along each axis, from the blue strokes."""
    a = np.asarray(page.get_pixmap(dpi=72).pil_image().convert("RGB")).astype(int)
    blue = ((a[:, :, 2] > 120) & (a[:, :, 2] - a[:, :, 0] > 60)
            & (a[:, :, 2] - a[:, :, 1] > 60))[SKIP_TOP_PX:, :]

    def runs(profile, thresh=15, gap=14):
        idx = np.nonzero(profile >= thresh)[0]
        out = []
        for i in idx:
            if out and i - out[-1][-1] <= gap:
                out[-1].append(i)
            else:
                out.append([i])
        return [float(np.average(r, weights=profile[r])) for r in out]

    rows = np.array(runs(blue.sum(1))) + SKIP_TOP_PX
    if len(rows) != 17:
        raise SystemExit("expected 17 row boundaries, found %d" % len(rows))

    # linearity / order check against cumulative length
    obs = np.diff(rows)
    exp = np.diff(cum) / cum[-1] * (rows[-1] - rows[0])
    err = np.abs(obs - exp).max() / (rows[-1] - rows[0])
    if err > 0.005:
        raise SystemExit("box spacing does not match chrom.sizes (%.3f%%)" % (100 * err))

    # columns: match candidate edges to the length-predicted positions
    cand = np.array(runs(blue.sum(0)))
    pred = cand[0] + cum / cum[-1] * (cand[-1] - cand[0])
    matched = np.array([cand[np.argmin(abs(cand - p))] for p in pred])
    cols = np.polyval(np.polyfit(cum, matched, 1), cum)
    return rows, cols, 100 * err


def main():
    cs = pd.read_csv(SIZES, sep="\t", header=None, names=["chrom", "length"])
    cs = cs.sort_values("chrom").reset_index(drop=True)
    cum = np.r_[0, np.cumsum(cs["length"].to_numpy())]

    src = pymupdf.open(SRC)
    rows, cols, err = boundaries(src[0], cum)

    w, h = src[0].rect.width, src[0].rect.height
    out = pymupdf.open()
    page = out.new_page(width=w, height=h + BAND)
    page.show_pdf_page(pymupdf.Rect(0, BAND, w, h + BAND), src, 0)

    grey = (0.35, 0.35, 0.35)
    page.draw_line(pymupdf.Point(cols[0], BAND - TICK - 2),
                   pymupdf.Point(cols[0], BAND - 2), color=grey, width=1.6)
    for i, name in enumerate(cs["chrom"]):
        mid = (cols[i] + cols[i + 1]) / 2
        page.draw_line(pymupdf.Point(cols[i + 1], BAND - TICK - 2),
                       pymupdf.Point(cols[i + 1], BAND - 2), color=grey, width=1.6)
        page.insert_text(pymupdf.Point(mid + FONTSIZE * 0.36, BAND - 26), str(name),
                         fontsize=FONTSIZE, fontname="helv", rotate=90, color=(0, 0, 0))

    out.save(OUT)
    out.close()
    print("box spacing agrees with chrom.sizes to %.3f%% of map width" % err)
    print("labelled %d chromosomes -> %s" % (len(cs), OUT))


if __name__ == "__main__":
    main()
