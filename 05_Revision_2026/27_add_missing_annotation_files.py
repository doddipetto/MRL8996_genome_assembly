#!/usr/bin/env python3
"""
27_add_missing_annotation_files.py
==================================
Completes the repeat/annotation part of the deposition package.

Two files the package claims or implies but does not deliver:

  1. 03_repeats/MRL8996_simple_repeats.bed is 0 bytes. The simple/low-complexity
     repeat positions were never written out, although the underlying Tandem
     Repeat Finder run exists and its output is what fed the simple-repeat
     circos track. The BED is regenerated here from the TRF .dat with the same
     parsing used by 09_TE_annotation/generate_final_circos_tracks.py (1-based
     inclusive TRF coordinates -> 0-based half-open BED), intervals merged per
     chromosome so the file is non-redundant.

  2. The 20-kb density tracks for transposable elements and for simple repeats
     are absent from 01_assembly_and_annotation/tracks/, which carries only GC,
     gene and exon density. Both are needed to reproduce Figure 1B and the
     within-chromosome TE percentile of the centromere calls (Table 3). They are
     copied from the circos track directory under names matching the tracks
     already deposited.

The MANIFEST rows for the three files are rewritten with the new sizes and md5s.

Inputs (read-only):
    MRL8996_HiC_data/09_TE_annotation/MRL8996_clean_for_EDTA.fasta.2.7.7.80.10.50.500.dat
    MRL8996_HiC_data/11_Circos_Plot/circos_track_TE_20kb.txt
    MRL8996_HiC_data/11_Circos_Plot/circos_track_SR_20kb.txt
    Revision_2026/deposition/01_assembly_and_annotation/tracks/chrom.sizes

Outputs:
    deposition/03_repeats/MRL8996_simple_repeats.bed                 (rewritten)
    deposition/01_assembly_and_annotation/tracks/track_TE_density_20kb.txt
    deposition/01_assembly_and_annotation/tracks/track_simple_repeat_density_20kb.txt
    deposition/MANIFEST.tsv                                          (rows updated)

Run: python3 scripts/27_add_missing_annotation_files.py      (env: python)
"""

import hashlib
import shutil
from pathlib import Path

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
DEP = ADODDI / "MRL8996_genome_paper" / "Revision_2026" / "deposition"

TRF = HIC / "09_TE_annotation" / "MRL8996_clean_for_EDTA.fasta.2.7.7.80.10.50.500.dat"
TRACKS_SRC = HIC / "11_Circos_Plot"
TRACKS_DST = DEP / "01_assembly_and_annotation" / "tracks"

DESCRIPTIONS = {
    "03_repeats/MRL8996_simple_repeats.bed":
        "Positions of the simple and low-complexity (tandem) repeats, from "
        "Tandem Repeat Finder, merged per chromosome.",
    "01_assembly_and_annotation/tracks/track_TE_density_20kb.txt":
        "Transposable-element density in 20-kb windows, as the masked fraction "
        "of each window; the track plotted in Figure 1B and used for the "
        "within-chromosome TE percentile of the centromere calls.",
    "01_assembly_and_annotation/tracks/track_simple_repeat_density_20kb.txt":
        "Simple-repeat density in 20-kb windows, as the masked fraction of each "
        "window; the track plotted in Figure 1B.",
}


def chrom_sizes():
    sizes = {}
    for line in (TRACKS_DST / "chrom.sizes").read_text().splitlines():
        if line.strip():
            c, n = line.split()[:2]
            sizes[c] = int(n)
    return sizes


def trf_intervals(path):
    """TRF .dat -> {chrom: [(start0, end_excl), ...]} in BED coordinates."""
    per_chrom = {}
    chrom = None
    for line in path.read_text().splitlines():
        if line.startswith("Sequence:"):
            chrom = line.split()[1]
        elif chrom and line[:1].isdigit():
            parts = line.split()
            start, end = int(parts[0]) - 1, int(parts[1])
            per_chrom.setdefault(chrom, []).append((start, end))
    return per_chrom


def merge(intervals, limit):
    out = []
    for start, end in sorted(intervals):
        start, end = max(start, 0), min(end, limit)
        if end <= start:
            continue
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_simple_repeats(sizes):
    raw = trf_intervals(TRF)
    dst = DEP / "03_repeats" / "MRL8996_simple_repeats.bed"
    n_int = n_bp = 0
    with open(dst, "w") as out:
        for chrom in sorted(sizes):
            for start, end in merge(raw.get(chrom, []), sizes[chrom]):
                out.write(f"{chrom}\t{start}\t{end}\tsimple_repeat\n")
                n_int += 1
                n_bp += end - start
    skipped = sorted(set(raw) - set(sizes))
    return dst, n_int, n_bp, skipped


def copy_track(src_name, dst_name, sizes):
    src, dst = TRACKS_SRC / src_name, TRACKS_DST / dst_name
    shutil.copy2(src, dst)
    bad = {line.split("\t")[0] for line in dst.read_text().splitlines()
           if line.strip()} - set(sizes)
    return dst, bad


def update_manifest(entries):
    path = DEP / "MANIFEST.tsv"
    lines = path.read_text().splitlines()
    header, rows = lines[0], [ln.split("\t") for ln in lines[1:] if ln.strip()]
    index = {r[0]: r for r in rows}
    for rel, abspath in entries.items():
        row = [rel, str(abspath.stat().st_size), md5(abspath), DESCRIPTIONS[rel]]
        if rel in index:
            index[rel][:] = row
        else:
            rows.append(row)
    rows.sort(key=lambda r: r[0])
    path.write_text(header + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    return len(rows)


def main():
    sizes = chrom_sizes()
    bed, n_int, n_bp, skipped = write_simple_repeats(sizes)
    te, te_bad = copy_track("circos_track_TE_20kb.txt", "track_TE_density_20kb.txt", sizes)
    sr, sr_bad = copy_track("circos_track_SR_20kb.txt",
                            "track_simple_repeat_density_20kb.txt", sizes)
    genome = sum(sizes.values())
    n_rows = update_manifest({
        "03_repeats/MRL8996_simple_repeats.bed": bed,
        "01_assembly_and_annotation/tracks/track_TE_density_20kb.txt": te,
        "01_assembly_and_annotation/tracks/track_simple_repeat_density_20kb.txt": sr,
    })
    print(f"simple_repeats.bed\t{n_int} intervals\t{n_bp} bp"
          f"\t{100 * n_bp / genome:.3f}% of {genome} bp")
    print(f"skipped_contigs\t{skipped if skipped else 'none'}")
    print(f"track_chrom_names_not_in_chrom_sizes\tTE={te_bad or 'none'}\tSR={sr_bad or 'none'}")
    print(f"manifest_rows\t{n_rows}")


if __name__ == "__main__":
    main()
