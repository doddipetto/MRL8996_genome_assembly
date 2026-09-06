#!/usr/bin/env python3
"""Rebuild the supplementary-table workbook after the deposition (editor E7).

The twelve tables of the submitted workbook (`media-1 (3).xlsx`) are now Data
File 1-12 of the repository deposition and are no longer distributed with the
manuscript, so the workbook is rebuilt from scratch and holds only the three
small tables added in this revision, which the reader consults directly.

The submitted workbook is read, never written.

out: Revision_2026/tables/Supplementary_Tables_MRL8996_revised.xlsx
"""
import csv
from pathlib import Path

import openpyxl
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

REV = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper"
           "/Revision_2026")
T = REV / "tables"
OUT = T / "Supplementary_Tables_MRL8996_revised.xlsx"

# (sheet name, title, [(sub-heading or None, tsv)])
SHEETS = [
    ("Supplementary Table 1",
     "Supplementary Table 1. Candidate centromeric regions of the 16 "
     "pseudomolecules of Fusarium oxysporum MRL8996 and the signals supporting "
     "each call (GC minimum, cis-interaction focus in the Hi-C map, "
     "StainedGlass self-identity band).",
     [(None, "Table_centromeric_interaction_blocks.tsv")]),
    ("Supplementary Table 2",
     "Supplementary Table 2. Content of each pseudomolecule of Fusarium "
     "oxysporum MRL8996: length, gene number and density, coding fraction, "
     "repeat-masked fraction, effector candidates, biosynthetic gene clusters "
     "and the assignment to the core or the accessory compartment.",
     [(None, "Table_per_chromosome_content.tsv")]),
    ("Supplementary Table 3",
     "Supplementary Table 3. Three-way comparison of the effector catalogues "
     "of MRL8996, Fol4287 and Fo47 after clustering at 70% identity, and the "
     "candidates found only in MRL8996.",
     [("Clusters and proteins in each region of the comparison",
       "venn_counts.tsv"),
      ("Effector candidates of MRL8996 with cluster, compartment, structural "
       "family and network assignment; the column Unique_to_MRL8996 marks the 135 "
       "candidates found only in this isolate",
       "Table_MRL8996_unique_effectors.tsv")]),
]

BOLD = Font(bold=True)


def main():
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, title, blocks in SHEETS:
        ws = wb.create_sheet(name)
        ws.cell(row=1, column=1, value=title).font = BOLD
        r = 3
        widths = {}
        for sub, fn in blocks:
            path = T / fn
            if not path.exists():
                raise SystemExit("missing %s" % path)
            if sub:
                ws.cell(row=r, column=1, value=sub).font = BOLD
                r += 1
            with open(path, newline="") as fh:
                for i, row in enumerate(csv.reader(fh, delimiter="\t")):
                    for c, val in enumerate(row, start=1):
                        try:
                            val2 = float(val) if val not in ("", "NA") else val
                        except ValueError:
                            val2 = val
                        cell = ws.cell(row=r, column=c, value=val2)
                        if i == 0:
                            cell.font = BOLD
                        widths[c] = min(46, max(widths.get(c, 10),
                                                len(str(val)) + 2))
                    r += 1
            r += 2
        for c, w in widths.items():
            ws.column_dimensions[get_column_letter(c)].width = w
        ws.freeze_panes = "A4"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(OUT))
    print("written", OUT)
    chk = openpyxl.load_workbook(str(OUT), read_only=True)
    for w2 in chk.worksheets:
        print("  %-24s %4d x %2d" % (w2.title, w2.max_row, w2.max_column))


if __name__ == "__main__":
    main()
