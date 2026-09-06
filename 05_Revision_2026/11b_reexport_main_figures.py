#!/usr/bin/env python3
"""
11b_reexport_main_figures.py
============================
Re-exports the four main figures at publication resolution for the revision,
answering the reviewers' request for higher-quality figures.

Every panel is taken from the ORIGINAL vector source (PDF) wherever one exists,
so nothing is upsampled from a screenshot: panels are placed with
show_pdf_page(), which keeps them vectorial and zoom-proof. Only the ChimeraX
structure renders inside the effectorome figure are intrinsically raster.

Figure 1 is re-assembled because panel A is new: the Hi-C contact map rendered
by 11a with the candidate centromeric interaction blocks annotated. Figures 2 and 3
already exist as vector PDFs and are re-exported unchanged; Figure 4 is the
recomposition from script 22 (vector PDF copy +
600-dpi PNG for the submission system), so no analysis is repeated.

Inputs (read-only originals):
    Revision_2026/figures/Figure_1A_HiC_labelled.pdf             (script 16:
        the cropped original below, with Chr01-Chr16 labelled)
    Revision_2026/figures/Figure_1A_HiC_original_cropped.pdf     (the original
        Juicebox export 2026.06.17.15.32.23.HiCImage.pdf, cropped to the map
        area only; vector drawing operators preserved, nothing rasterised)
    Revision_2026/figures/Figure_1B_circos_centromeres.pdf       (script 15)
    MRL8996_HiC_data/11_Circos_Plot/MRL8996_Circular_Plot_Final_Gap.pdf
    MRL8996_HiC_data/12_DGenies_Synteny/MRL8996_Synteny_Dotplot_Python.pdf
    MRL8996_HiC_data/12_DGenies_Synteny/MRL8996_Nx_curve.pdf
    MRL8996_HiC_data/06_circos_synteny/Synteny_3Way_MRL_Hub.pdf
    MRL8996_HiC_data/06_circos_synteny/Synteny_MRL8996_Self.pdf
    Revision_2026/figures/Figure_4_Effectorome_recomposed.pdf     (script 22:
        the published master figure with the revised structural network,
        singletons retained, in panel B; panels A and C unchanged)

Outputs (Revision_2026/figures/final/):
    Figure_1_Assembly_Architecture.pdf / .png
    Figure_2_Assembly_Comparison.pdf / .png
    Figure_3_Synteny.pdf / .png
    Figure_4_Effectorome.pdf / .png
    figure_panel_provenance.tsv     source file + type behind every panel

Usage: python3 11b_reexport_main_figures.py [--dpi 600]
"""

import argparse
from pathlib import Path

import pymupdf as fitz

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
PAPER = ADODDI / "MRL8996_genome_paper"
REV = PAPER / "Revision_2026"
FIG = REV / "figures"
OUT = FIG / "final"

MM = 2.834645669          # points per mm
PAGE_W_MM = 180           # two-column width; vector sharpness is resolution-free
MARGIN_MM = 4
GAP_MM = 3

# figure name -> list of rows; each row is a list of (label, source, weight)
FIGURES = {
    "Figure_1_Assembly_Architecture": [
        [("A", FIG / "Figure_1A_HiC_labelled.pdf", 1.0),
         ("B", FIG / "Figure_1B_circos_centromeres.pdf", 1.0)],
    ],
    "Figure_2_Assembly_Comparison": [
        [("A", HIC / "12_DGenies_Synteny" / "MRL8996_Synteny_Dotplot_Python.pdf", 1.0),
         ("B", HIC / "12_DGenies_Synteny" / "MRL8996_Nx_curve.pdf", 1.0)],
    ],
    "Figure_3_Synteny": [
        [("A", HIC / "06_circos_synteny" / "Synteny_3Way_MRL_Hub.pdf", 1.0),
         ("B", HIC / "06_circos_synteny" / "Synteny_MRL8996_Self.pdf", 1.0)],
    ],
    "Figure_4_Effectorome": [
        [("", FIG / "Figure_4_Effectorome_recomposed.pdf", 1.0)],
    ],
}


def fit(outer, sw, sh):
    """Centred sub-rectangle of `outer` preserving the sw:sh aspect ratio."""
    if sw <= 0 or sh <= 0:
        return outer
    ow, oh = outer.width, outer.height
    scale = min(ow / sw, oh / sh)
    w, h = sw * scale, sh * scale
    x = outer.x0 + (ow - w) / 2
    y = outer.y0 + (oh - h) / 2
    return fitz.Rect(x, y, x + w, y + h)


def place(page, rect, src_path):
    src = fitz.open(str(src_path))
    sp = src[0].rect
    page.show_pdf_page(fit(rect, sp.width, sp.height), src, 0)
    src.close()
    return sp.width / sp.height


def build(name, rows, dpi, prov):
    for row in rows:
        for _, p, _ in row:
            if not Path(p).exists():
                raise SystemExit("missing panel source: %s" % p)

    # page height from the panel aspect ratios, so nothing is distorted
    page_w = PAGE_W_MM * MM
    margin, gap = MARGIN_MM * MM, GAP_MM * MM
    row_h = []
    for row in rows:
        wsum = sum(w for _, _, w in row)
        avail = page_w - 2 * margin - gap * (len(row) - 1)
        h = 0.0
        for _, p, w in row:
            src = fitz.open(str(p))
            sp = src[0].rect
            src.close()
            pw = avail * w / wsum
            h = max(h, pw * sp.height / sp.width)
        row_h.append(h)
    page_h = 2 * margin + sum(row_h) + gap * (len(rows) - 1)

    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)
    y = margin
    for row, h in zip(rows, row_h):
        wsum = sum(w for _, _, w in row)
        avail = page_w - 2 * margin - gap * (len(row) - 1)
        x = margin
        for label, p, w in row:
            pw = avail * w / wsum
            place(page, fitz.Rect(x, y, x + pw, y + h), p)
            if label:
                page.insert_text((x, y + 8), label, fontsize=11,
                                 fontname="hebo", color=(0, 0, 0))
            prov.append((name, label or "-", str(p), Path(p).suffix.lstrip("."),
                         "vector"))
            x += pw + gap
        y += h + gap

    OUT.mkdir(parents=True, exist_ok=True)
    pdf_out = OUT / (name + ".pdf")
    doc.save(str(pdf_out), deflate=True, garbage=4)
    pix = doc[0].get_pixmap(dpi=dpi)
    pix.save(str(OUT / (name + ".png")))
    doc.close()
    print("%-34s %5.1f x %5.1f mm  ->  %s (%.1f MB) + %d x %d px PNG"
          % (name, page_w / MM, page_h / MM, pdf_out.name,
             pdf_out.stat().st_size / 1e6, pix.width, pix.height))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=600)
    a = ap.parse_args()
    prov = []
    for name, rows in FIGURES.items():
        build(name, rows, a.dpi, prov)
    with open(OUT / "figure_panel_provenance.tsv", "w") as fh:
        fh.write("figure\tpanel\tsource_file\tsource_format\tembedding\n")
        for r in prov:
            fh.write("\t".join(r) + "\n")
    print("panels re-exported: %d" % len(prov))


if __name__ == "__main__":
    main()
