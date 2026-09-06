#!/usr/bin/env python3
"""
22_figure4_recompose_network.py
===============================
Rebuilds the effectorome figure (Figure 4) with the revised structural
network in panel B.

Panels A and C are taken from the published master figure and are not
redrawn: the master page is placed twice, clipped to the panel A box and to
the panel C strip, so the TM-score heatmap and the ChimeraX renders keep
their original vector/raster content and their panel letters.  Panel B is
replaced by the network from script 07 run with --panel (no title; the
caption carries it), which retains the 131 structural singletons and the 44
sub-threshold communities as outer rings and rings the MRL8996-unique
effectors in red.  The "B" letter is the master's own
glyph, placed from the master page like the panels, so all three letters are
the same face and weight.

Inputs:
    MRL8996_genome_paper/Figure_1_Master_Effectorome.pdf          (panels A, C)
    Revision_2026/figures/Fig_structural_network_singletons_panel.pdf
        (script 07 --panel)

Outputs:
    Revision_2026/figures/Figure_4_Effectorome_recomposed.pdf
    Revision_2026/figures/Figure_4_Effectorome_recomposed.png   (600 dpi)

Usage: python3 22_figure4_recompose_network.py [--dpi 600]
"""

import argparse
from pathlib import Path

import pymupdf as fitz

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
PAPER = ADODDI / "MRL8996_genome_paper"
REV = PAPER / "Revision_2026"
MASTER = PAPER / "Figure_1_Master_Effectorome.pdf"
NETWORK = REV / "figures" / "Fig_structural_network_singletons_panel.pdf"
OUT = REV / "figures" / "Figure_4_Effectorome_recomposed"

# panel boxes on the master page (850.4 x 566.9 pt), measured from its own
# content: panel A drawings span x 38-366 / y 14-357 with the TM-score and
# AMAPEC keys out to x ~430; the box stops at 438 so the master's own "B"
# (x 449-469) is not clipped into panel A. The panel C strip starts at the
# "C" letter.
BOX_A = fitz.Rect(0, 0, 438, 362)
BOX_C = fitz.Rect(0, 362, 850.393677, 566.929138)

PAGE_W = 850.393677       # keep the master's width
TOP_H = 360.0             # height of the panel A / panel B row
GAP = 10.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    src = fitz.open(MASTER)
    net = fitz.open(NETWORK)
    nw, nh = net[0].rect.width, net[0].rect.height

    a_w = TOP_H * BOX_A.width / BOX_A.height
    b_w = TOP_H * nw / nh
    c_w = min(PAGE_W, BOX_C.width)
    c_h = c_w * BOX_C.height / BOX_C.width
    page_h = TOP_H + GAP + c_h

    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=page_h)

    # panel A: master page clipped to the heatmap box
    page.show_pdf_page(fitz.Rect(0, 0, a_w, TOP_H), src, 0, clip=BOX_A)

    # panel B: revised structural network, right-aligned in the remaining space
    b_x = a_w + GAP + max(0.0, (PAGE_W - a_w - GAP - b_w) / 2)
    page.show_pdf_page(fitz.Rect(b_x, 0, b_x + b_w, TOP_H), net, 0)
    page.show_pdf_page(fitz.Rect(b_x + 2, 6, b_x + 2 + 19.6, 43.5),
                       src, 0, clip=fitz.Rect(449.0, 9.4, 468.6, 46.9))

    # panel C: master page clipped to the structure strip, letter included
    page.show_pdf_page(fitz.Rect(0, TOP_H + GAP, c_w, TOP_H + GAP + c_h),
                       src, 0, clip=BOX_C)

    doc.save(str(OUT) + ".pdf", deflate=True, garbage=3)
    page.get_pixmap(dpi=args.dpi).save(str(OUT) + ".png")
    print("page %.1f x %.1f pt ; A %.1f wide, B %.1f wide at x %.1f, C %.1f x %.1f"
          % (PAGE_W, page_h, a_w, b_w, b_x, c_w, c_h))


if __name__ == "__main__":
    main()
