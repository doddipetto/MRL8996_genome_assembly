#!/usr/bin/env python3
"""
14_call_centromeres_stainedglass.py
===================================
Calls one candidate centromeric region per chromosome from the per-chromosome
StainedGlass self-identity data produced by 13_stainedglass_percontig.py, and
scores each call against two independent tracks that were computed for the
circos figure.

Rationale
---------
Fusarium centromeres are 100-200 kb, gene-poor, LTR-retrotransposon-rich blocks
whose constituent elements are mutually similar but appreciably diverged. In a
StainedGlass identity heatmap they therefore appear as a compact off-diagonal
cluster at ~85-93 % identity, distinct both from the near-identical (>95 %)
recent TE arrays that are scattered along the accessory chromosomes and from the
subtelomeric arrays that sit in the terminal few hundred kb.

The caller reproduces that reading quantitatively. For each chromosome it takes
window pairs that are off-diagonal (offset > 20 kb), keeps those with
perID_by_events in [84, 93], counts both members of every pair into 20-kb bins,
smooths over five bins (100 kb) and takes the highest interior bin, excluding the
terminal 200 kb so that subtelomeric arrays cannot win. The block is then grown
outwards while the smoothed density stays above a quarter of the peak.

Support columns:
    gene_pctile  percentile of the block's mean gene density among all 20-kb
                 bins of that chromosome (low is centromere-like)
    TE_pctile    same for TE density (high is centromere-like)
    hic_call_bp  the Hi-C candidate block from
                 Table_centromeric_interaction_blocks.tsv, for comparison

Chromosome naming: 11_Circos_Plot/chrom.sizes orders chromosomes by descending
length, which is not the PGA scaffold numbering (scaffold_3 and scaffold_4 are
swapped). The mapping used here is by length rank and reproduces chrom.sizes and
every 20-kb track exactly (Chr12 differs by 596 bp, a gap-closing difference).

Inputs (read-only):
    Revision_2026/results/stainedglass/ChrNN.2000.10000.bed.gz   (script 13)
    MRL8996_HiC_data/11_Circos_Plot/chrom.sizes
    MRL8996_HiC_data/11_Circos_Plot/track_gene_density_20kb.txt
    MRL8996_HiC_data/11_Circos_Plot/circos_track_TE_20kb.txt
    Revision_2026/tables/Table_centromeric_interaction_blocks.tsv

Outputs:
    Revision_2026/tables/Table_centromeres_stainedglass.tsv
    Revision_2026/results/centromeres_stainedglass.bed   (circos track)
    Revision_2026/figures/Fig_stainedglass_identity_16chr.pdf / .png (600 dpi)

Usage: python3 14_call_centromeres_stainedglass.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                     "pdf.fonttype": 42, "axes.linewidth": 0.6})

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
SGD = REV / "results" / "stainedglass"
CIRC = ADODDI / "MRL8996_HiC_data" / "11_Circos_Plot"

W = 2000
BIN = 20000
ID_LO, ID_HI = 84.0, 93.0
MIN_OFFSET = 20000
TERMINAL = 200000
SMOOTH = 5
EDGE_FRAC = 0.25
CHRS = ["Chr%02d" % i for i in range(1, 17)]


def load_lengths():
    cs = pd.read_csv(CIRC / "chrom.sizes", sep="\t", header=None, names=["chr", "len"])
    return dict(zip(cs["chr"], cs["len"].astype(int)))


def load_bed(chrom):
    b = pd.read_csv(SGD / f"{chrom}.{W}.10000.bed.gz", sep="\t")
    b.columns = [c.lstrip("#") for c in b.columns]
    return b


def call_one(chrom, bed, length):
    n = int(np.ceil(length / BIN)) + 1
    g = bed[((bed["query_start"] - bed["reference_start"]).abs() > MIN_OFFSET)
            & bed["perID_by_events"].between(ID_LO, ID_HI)]
    v = np.zeros(n)
    for s in list(g["query_start"]) + list(g["reference_start"]):
        v[min(int(s // BIN), n - 1)] += 1
    sm = np.convolve(v, np.ones(SMOOTH) / SMOOTH, mode="same")
    keep = np.ones(n, bool)
    k = TERMINAL // BIN
    keep[:k] = False
    keep[-(k + 1):] = False
    if sm[keep].max() == 0:
        return None
    i = np.arange(n)[keep][np.argmax(sm[keep])]
    thr = EDGE_FRAC * sm[i]
    a, z = i, i
    while a > 0 and sm[a - 1] >= thr:
        a -= 1
    while z < n - 1 and sm[z + 1] >= thr:
        z += 1
    start, end = a * BIN, int(min((z + 1) * BIN, length))
    sel = g[(g["query_start"] >= start) & (g["query_start"] < end)]
    return dict(chrom=chrom, start=start, end=end,
                size_kb=(end - start) / 1000,
                peak_bp=int(i * BIN),
                relative_position=round((start + end) / 2 / length, 3),
                n_window_pairs=int(v[a:z + 1].sum()),
                median_identity=round(float(sel["perID_by_events"].median()), 1) if len(sel) else np.nan)


def track_percentile(track, chrom, start, end):
    x = track[track["chr"] == chrom]
    ins = x[(x["s"] >= start) & (x["s"] < end)]["v"]
    if not len(ins) or not len(x):
        return np.nan
    return round(100.0 * float((x["v"] < ins.mean()).mean()), 0)


def figure(beds, lengths, calls, out):
    fig, axes = plt.subplots(4, 4, figsize=(7.09, 7.3))
    for ax, chrom in zip(axes.ravel(), CHRS):
        b = beds[chrom]
        Lm = lengths[chrom] / 1e6
        x = np.r_[b["query_start"], b["reference_start"]] / 1e6
        y = np.r_[b["reference_start"], b["query_start"]] / 1e6
        pid = np.r_[b["perID_by_events"], b["perID_by_events"]]
        sc = ax.scatter(x, y, c=pid, s=1.1, cmap="Spectral_r", vmin=80, vmax=100, lw=0,
                        rasterized=True)
        r = calls.loc[calls["chrom"] == chrom]
        if len(r):
            s, e = float(r["start"].iloc[0]) / 1e6, float(r["end"].iloc[0]) / 1e6
            ax.add_patch(plt.Rectangle((s, s), e - s, e - s, fill=False,
                                       ec="#1a1a1a", lw=0.7, zorder=5))
            ax.axvspan(s, e, color="#1f5c99", alpha=0.10, zorder=0)
            ax.axhspan(s, e, color="#1f5c99", alpha=0.10, zorder=0)
        ax.set_xlim(0, Lm)
        ax.set_ylim(0, Lm)
        ax.set_aspect(1)
        ax.set_title("%s  %.2f Mb" % (chrom, Lm), fontsize=6.5, pad=2)
        ax.tick_params(labelsize=5, length=2, pad=1)
    fig.subplots_adjust(left=0.05, right=0.88, top=0.96, bottom=0.05, wspace=0.35, hspace=0.4)
    cax = fig.add_axes([0.905, 0.35, 0.016, 0.30])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("pairwise identity (%)", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=600)
    plt.close(fig)


def main():
    lengths = load_lengths()
    beds = {c: load_bed(c) for c in CHRS}
    calls = pd.DataFrame([r for r in (call_one(c, beds[c], lengths[c]) for c in CHRS) if r])

    gd = pd.read_csv(CIRC / "track_gene_density_20kb.txt", sep="\t", header=None,
                     names=["chr", "s", "e", "v"])
    td = pd.read_csv(CIRC / "circos_track_TE_20kb.txt", sep="\t", header=None,
                     names=["chr", "s", "e", "v"])
    calls["gene_density_pctile"] = [track_percentile(gd, r.chrom, r.start, r.end)
                                    for r in calls.itertuples()]
    calls["TE_density_pctile"] = [track_percentile(td, r.chrom, r.start, r.end)
                                  for r in calls.itertuples()]

    hic = pd.read_csv(REV / "tables" / "Table_centromeric_interaction_blocks.tsv", sep="\t")
    hmap = dict(zip(hic["chrom"], hic["centromere_call_bp"]))
    calls["hic_candidate_bp"] = calls["chrom"].map(hmap)
    mid = (calls["start"] + calls["end"]) / 2
    calls["hic_within_200kb"] = np.where(
        calls["hic_candidate_bp"].notna(),
        (calls["hic_candidate_bp"] - mid).abs() <= 200000, np.nan)
    calls["support"] = ((calls["gene_density_pctile"] <= 25).astype(int)
                        + (calls["TE_density_pctile"] >= 75).astype(int))

    (REV / "tables").mkdir(exist_ok=True)
    (REV / "figures").mkdir(exist_ok=True)
    calls.to_csv(REV / "tables" / "Table_centromeres_stainedglass.tsv", sep="\t", index=False)
    calls[["chrom", "start", "end"]].to_csv(
        REV / "results" / "centromeres_stainedglass.bed", sep="\t", index=False, header=False)
    figure(beds, lengths, calls, REV / "figures" / "Fig_stainedglass_identity_16chr")

    print(calls[["chrom", "start", "end", "size_kb", "median_identity",
                 "gene_density_pctile", "TE_density_pctile", "support"]].to_string(index=False))
    print("\ncalls with both supports: %d / %d" % (int((calls["support"] == 2).sum()), len(calls)))


if __name__ == "__main__":
    main()
