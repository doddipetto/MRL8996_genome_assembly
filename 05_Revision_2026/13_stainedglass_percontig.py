#!/usr/bin/env python3
"""
13_stainedglass_percontig.py
============================
Re-runs the StainedGlass self-identity analysis one chromosome at a time, at
2-kb windows, so that the centromeric high-identity arrays are resolved.

The genome-wide run in
    MRL8996_HiC_data/02_telomeres_centromeres/StainedGlass/results/small.5000.10000.*
used 5-kb windows over all 16 pseudomolecules at once and retained only 21,818
window pairs; its intra-chromosomal off-diagonal signal is dominated by
subtelomeric arrays and gives no usable centromere position on the core
chromosomes. Running each chromosome separately removes the cross-chromosome
minimizer competition and quadruples the window resolution.

The pipeline reproduces the StainedGlass rules exactly (workflow/Snakefile:
make_windows -> window_fa -> aln -> sort_aln -> identity -> pair_end_bed) with
the tool's own default parameters for W = 2000:

    bedtools makewindows -w 2000
    bedtools getfasta
    minimap2 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx
    samtools sort
    workflow/scripts/samIdentity.py --matches 400 --header
    workflow/scripts/refmt.py --window 2000

so the perID_by_events values are computed by StainedGlass's own code, not a
re-implementation. Nothing is written inside the user's StainedGlass directory.

Outputs (Revision_2026/results/stainedglass/):
    ChrNN.2000.10000.bed.gz        pair-end bed, one row per window pair
    ChrNN.2000.10000.full.tbl.gz   full identity table
    ChrNN.2000.10000.sorted.bam    self-alignment

Usage: python3 13_stainedglass_percontig.py [--threads 12] [--window 2000]
"""

import argparse
import subprocess as sp
from pathlib import Path

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
GENOME = ADODDI / "MRL8996_HiC_data" / "mrl8996-hic.chromosomes.fasta"
SG = ADODDI / "MRL8996_HiC_data" / "02_telomeres_centromeres" / "StainedGlass"
SCRIPTS = SG / "workflow" / "scripts"
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
OUT = REV / "results" / "stainedglass"
F = 10000


def sh(cmd, **kw):
    sp.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def scaffold_to_chr(fai):
    """chrom.sizes in 11_Circos_Plot names chromosomes by descending length,
    which is NOT the scaffold numbering (scaffold_3 and scaffold_4 are swapped).
    Map by length rank so the names match every other track in the paper."""
    rows = [l.split("\t")[:2] for l in open(fai)]
    rows = [(n, int(L)) for n, L in rows]
    order = sorted(rows, key=lambda t: -t[1])
    return {n: "Chr%02d" % (i + 1) for i, (n, L) in enumerate(order)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--window", type=int, default=2000)
    a = ap.parse_args()
    W, S = a.window, a.window // 5

    OUT.mkdir(parents=True, exist_ok=True)
    tmp = OUT / "temp"
    tmp.mkdir(exist_ok=True)

    name2chr = scaffold_to_chr(str(GENOME) + ".fai")
    for scaf, chrom in sorted(name2chr.items(), key=lambda t: t[1]):
        bed = OUT / f"{chrom}.{W}.{F}.bed.gz"
        if bed.exists():
            print(chrom, "done, skipping", flush=True)
            continue
        fa = tmp / f"{chrom}.fasta"
        win = tmp / f"{chrom}.{W}.fasta"
        bam = OUT / f"{chrom}.{W}.{F}.sorted.bam"
        print("==", chrom, scaf, flush=True)
        sh(f"samtools faidx {GENOME} '{scaf}' > {fa} && samtools faidx {fa}")
        sh(f"bedtools makewindows -g {fa}.fai -w {W} > {tmp}/{chrom}.{W}.bed")
        sh(f"bedtools getfasta -fi {fa} -bed {tmp}/{chrom}.{W}.bed > {win}")
        sh(f"minimap2 -t {a.threads} -f {F} -s {S} -ax ava-ont --dual=yes --eqx "
           f"{win} {win} 2> {tmp}/{chrom}.minimap2.log "
           f"| samtools sort -m 2G -@ 4 --write-index -o {bam}")
        sh(f"python {SCRIPTS}/samIdentity.py --threads {a.threads} "
           f"--matches {S} --header {bam} | gzip > {tmp}/{chrom}.tbl.gz")
        sh(f"python {SCRIPTS}/refmt.py --window {W} --fai {fa}.fai "
           f"--full {OUT}/{chrom}.{W}.{F}.full.tbl.gz "
           f"{tmp}/{chrom}.tbl.gz {bed}")
        for f in (fa, win, Path(str(fa) + ".fai"), tmp / f"{chrom}.{W}.bed",
                  tmp / f"{chrom}.tbl.gz"):
            f.unlink(missing_ok=True)
    print("all chromosomes complete")


if __name__ == "__main__":
    main()
