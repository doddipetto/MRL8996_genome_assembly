#!/usr/bin/env python3
"""
20_effector_map_plus_venn.py
============================
Two-panel main figure: the author's own linear effector map (panel a) with
the CD-HIT pan-effectorome Venn diagram added as a panel on the right
(panel b).

Both panels are placed as vector graphics; neither is redrawn, so panel a
is the ggplot output of script 21 (the author's plot_map_effector_position.R
with the unique-effector layer added) and
panel b is the output of 02_venn_effectorome.py.  Only the panel letters
are drawn here.

Inputs:
    Revision_2026/figures/MRL8996_Linear_Effector_Map_unique.pdf  (script 21:
        the author's map with the MRL8996-unique effectors marked)
    Revision_2026/figures/Venn_Effectorome_MRL8996.pdf       (script 02)

Outputs:
    Revision_2026/figures/Figure_effector_map_and_venn.pdf
    Revision_2026/figures/Figure_effector_map_and_venn.png   (600 dpi)

Usage: python3 20_effector_map_plus_venn.py
"""
import io
from pathlib import Path

import pypdfium2 as pdfium
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas as rl_canvas

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
MAP = REV / "figures" / "MRL8996_Linear_Effector_Map_unique.pdf"
VENN = REV / "figures" / "Venn_Effectorome_MRL8996.pdf"
OUT = REV / "figures" / "Figure_effector_map_and_venn"

PAGE_W, PAGE_H = 11.0 * inch, 5.6 * inch
MAP_W = 6.95 * inch          # panel a target width; height follows the aspect
VENN_H = 4.35 * inch         # panel b target height
GAP = 0.10 * inch
PAD_L, PAD_T = 0.20 * inch, 0.26 * inch


def page_box(path):
    p = PdfReader(str(path)).pages[0]
    b = p.mediabox
    return float(b.width), float(b.height)


def letters():
    """A transparent overlay carrying just the two panel letters."""
    buf = io.BytesIO()
    c = rl_canvas.Canvas(buf, pagesize=(PAGE_W, PAGE_H))
    c.setFont("Helvetica-Bold", 13)
    c.drawString(PAD_L, PAGE_H - PAD_T + 4, "a")
    c.drawString(PAD_L + MAP_W + GAP, PAGE_H - PAD_T + 4, "b")
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def main():
    mw, mh = page_box(MAP)
    vw, vh = page_box(VENN)
    s_map = MAP_W / mw
    s_venn = VENN_H / vh

    writer = PdfWriter()
    page = writer.add_blank_page(width=PAGE_W, height=PAGE_H)

    map_h = mh * s_map
    map_y = (PAGE_H - map_h) / 2 - 0.05 * inch
    page.merge_transformed_page(
        PdfReader(str(MAP)).pages[0],
        Transformation().scale(s_map).translate(PAD_L, map_y))

    venn_w = vw * s_venn
    venn_x = PAD_L + MAP_W + GAP + max(0.0, (PAGE_W - PAD_L - MAP_W - GAP
                                             - 0.2 * inch - venn_w) / 2)
    venn_y = (PAGE_H - VENN_H) / 2
    page.merge_transformed_page(
        PdfReader(str(VENN)).pages[0],
        Transformation().scale(s_venn).translate(venn_x, venn_y))

    page.merge_page(letters())
    page.mediabox = RectangleObject((0, 0, PAGE_W, PAGE_H))
    with open(str(OUT) + ".pdf", "wb") as fh:
        writer.write(fh)

    pdf = pdfium.PdfDocument(str(OUT) + ".pdf")
    pdf[0].render(scale=600 / 72).to_pil().save(str(OUT) + ".png")
    print("panel a %.2f x %.2f in (scale %.3f), panel b %.2f x %.2f in (scale %.3f)"
          % (MAP_W / inch, map_h / inch, s_map, venn_w / inch, VENN_H / inch, s_venn))
    print("page %.2f x %.2f in -> %s.pdf/.png" % (PAGE_W / inch, PAGE_H / inch, OUT.name))


if __name__ == "__main__":
    main()
