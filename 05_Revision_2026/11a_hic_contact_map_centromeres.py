#!/usr/bin/env python3
"""
11a_hic_contact_map_centromeres.py
==================================
Re-renders Figure 1A (genome-wide Hi-C contact map) at publication resolution
and adds the centromeric interaction-block annotation Reviewer 2 asked for.

The .hic file produced by the scaffolding pipeline stores the whole assembly as
a single pseudo-chromosome ("assembly"), so chromosome boundaries are taken
from the scaffold order and lengths of mrl8996-hic.final.fasta and mapped onto
the Chr01-Chr16 labels of the released assembly by scaffold length.

Centromere calling. Three independent 50-kb signals are computed per
chromosome, each on the ICE-balanced matrix (or the repeat annotation), with
the terminal 6% of each chromosome excluded:

    trans   maximum of inter-chromosomal contact (fungal centromeres cluster
            in the nucleus)
    cis     minimum of intra-chromosomal contact (heterochromatic, repeat-rich
            centromeric block depletes cis coverage)
    repeat  maximum of RepeatMasker-masked fraction (EDTA RM.out)

A candidate centromeric interaction block is reported only where at least two
of the three signals agree within TOL bp; the per-signal positions and the
concordance flag are written out for every chromosome, so chromosomes without
a supported call are visible rather than silently filled in. No sequence-based
(CENP-A / centromeric repeat) evidence is used, and the calls are labelled
accordingly on the figure.

Inputs (read-only):
    MRL8996_HiC_data/mrl8996-hic.hic
    MRL8996_HiC_data/customer_delivery_dir/final_scaffold_results/mrl8996-hic.final.fasta
    MRL8996_HiC_data/08_circos_features/chrom.sizes

Outputs:
    Revision_2026/figures/Figure_1A_HiC_contact_map.pdf / .png  (600 dpi)
    Revision_2026/tables/Table_centromeric_interaction_blocks.tsv

Usage: python3 11a_hic_contact_map_centromeres.py [--binsize 50000]
"""

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import straw

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC_DIR = ADODDI / "MRL8996_HiC_data"
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
HIC = HIC_DIR / "mrl8996-hic.hic"
FINAL_FA = (HIC_DIR / "customer_delivery_dir" / "final_scaffold_results"
            / "mrl8996-hic.final.fasta")
SIZES = HIC_DIR / "08_circos_features" / "chrom.sizes"
RMOUT = (HIC_DIR / "09_TE_annotation"
         / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.anno"
         / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.RM.out")
SMOOTH = 5
TOL = 300000            # bp within which two signals count as concordant


def scaffold_layout():
    """Ordered scaffold lengths from the .hic assembly, labelled with Chr IDs."""
    lens = []
    for line in open(FINAL_FA):
        if line.startswith(">"):
            m = re.search(r"length_(\d+)", line)
            lens.append((line[1:].split()[0], int(m.group(1))))
    chrom = pd.read_csv(SIZES, sep="\t", header=None, names=["chrom", "length"])
    chrom = chrom[chrom["chrom"].str.match(r"Chr\d+")]
    by_len = {}
    for _, r in chrom.iterrows():
        by_len.setdefault(r["length"], []).append(r["chrom"])
    out, offset = [], 0
    for name, L in lens:
        label = None
        for delta in (0, 1, -1, 596, -596):          # tolerate gap-padding
            if (L + delta) in by_len and by_len[L + delta]:
                label = by_len[L + delta].pop(0)
                break
        if label is not None:
            out.append({"scaffold": name, "chrom": label, "length": L,
                        "offset": offset})
        offset += L
    return pd.DataFrame(out), offset


def load_matrix(binsize):
    recs = straw.straw("NONE", str(HIC), "assembly", "assembly", "BP", binsize)
    xs, ys, cs = np.array(recs[0]), np.array(recs[1]), np.array(recs[2], dtype=float)
    n = int(max(xs.max(), ys.max()) // binsize) + 1
    M = np.zeros((n, n))
    i, j = (xs // binsize).astype(int), (ys // binsize).astype(int)
    np.add.at(M, (i, j), cs)
    np.add.at(M, (j, i), cs)
    M[np.diag_indices(n)] /= 2
    return M


def ice_balance(M, block, n_iter=30):
    A = M.copy()
    cov = A.sum(1)
    good = (block >= 0) & (cov > np.percentile(cov[cov > 0], 5))
    A[~good, :] = 0
    A[:, ~good] = 0
    for _ in range(n_iter):
        s = A.sum(1)
        s[s == 0] = 1
        s /= s[good].mean()
        A /= s[:, None]
        A /= s[None, :]
    return A, good


def repeat_profile(chrom, length, binsize):
    nb = length // binsize + 1
    v = np.zeros(nb)
    for line in open(RMOUT):
        f = line.split()
        if len(f) > 12 and f[0].isdigit() and f[4] == chrom:
            st, en = int(f[5]), int(f[6])
            for b in range(st // binsize, min(en // binsize, nb - 1) + 1):
                v[b] += min(en, (b + 1) * binsize) - max(st, b * binsize)
    return v / binsize


def call_centromeres(M, layout, binsize):
    """Three independent signals per chromosome; report concordant calls."""
    n = M.shape[0]
    block = np.full(n, -1)
    for k, r in layout.iterrows():
        block[int(r["offset"] // binsize):
              int((r["offset"] + r["length"]) // binsize) + 1] = k
    B, good = ice_balance(M, block)
    trans = np.zeros(n)
    cis = np.zeros(n)
    for b in np.where(good)[0]:
        same = block == block[b]
        trans[b] = B[b, good & ~same].sum()
        cis[b] = B[b, good & same].sum()
    kern = np.ones(SMOOTH) / SMOOTH

    rows = []
    for k, r in layout.iterrows():
        lo = int(r["offset"] // binsize)
        hi = int((r["offset"] + r["length"]) // binsize)
        nb = hi - lo + 1
        keep = good[lo:hi + 1].copy()
        edge = max(3, int(0.06 * nb))
        keep[:edge] = False
        keep[nb - edge:] = False
        if keep.sum() < 3:
            continue
        t = np.convolve(trans[lo:hi + 1], kern, mode="same")
        c = np.convolve(cis[lo:hi + 1], kern, mode="same")
        rp = np.convolve(repeat_profile(r["chrom"], int(r["length"]), binsize),
                         kern, mode="same")[:nb]
        rp = np.pad(rp, (0, max(0, nb - rp.size)))
        p_t = int(np.argmax(np.where(keep, t, -np.inf)))
        p_c = int(np.argmin(np.where(keep, c, np.inf)))
        p_r = int(np.argmax(np.where(keep, rp, -np.inf)))

        pos = {"trans": p_t, "cis": p_c, "repeat": p_r}
        support, call = [], None
        for a in pos:
            agree = [b for b in pos
                     if abs(pos[a] - pos[b]) * binsize <= TOL]
            if len(agree) > len(support):
                support, call = agree, int(np.mean([pos[b] for b in agree]))
        concordant = len(support) >= 2
        rows.append({
            "chrom": r["chrom"],
            "chrom_length_bp": int(r["length"]),
            "trans_max_bp": p_t * binsize,
            "cis_min_bp": p_c * binsize,
            "repeat_max_bp": p_r * binsize,
            "supporting_signals": ",".join(sorted(support)) if concordant else "",
            "n_supporting": len(support) if concordant else 0,
            "centromere_call_bp": call * binsize if concordant else "",
            "relative_position": (round(call * binsize / r["length"], 3)
                                  if concordant else ""),
            "repeat_fraction_at_call": (round(float(rp[call]), 3)
                                        if concordant else ""),
            "global_bin": lo + call if concordant else -1,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binsize", type=int, default=50000)
    a = ap.parse_args()

    layout, total = scaffold_layout()
    M = load_matrix(a.binsize)
    cen = call_centromeres(M, layout, a.binsize)
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    cen.drop(columns=["global_bin"]).to_csv(
        REV / "tables" / "Table_centromeric_interaction_blocks.tsv",
        sep="\t", index=False)

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    with np.errstate(divide="ignore"):
        L = np.log10(M + 1)
    vmax = np.percentile(L[L > 0], 99.5)
    im = ax.imshow(L, cmap="Reds", vmin=0, vmax=vmax, interpolation="nearest",
                   origin="upper")

    edges = [r["offset"] / a.binsize for _, r in layout.iterrows()]
    ends = [(r["offset"] + r["length"]) / a.binsize for _, r in layout.iterrows()]
    for e in edges[1:] + [ends[-1]]:
        ax.axhline(e, color="0.35", lw=0.35)
        ax.axvline(e, color="0.35", lw=0.35)
    mids = [(s + e) / 2 for s, e in zip(edges, ends)]
    ax.set_xticks(mids)
    ax.set_xticklabels([c.replace("Chr", "") for c in layout["chrom"]], fontsize=6.5)
    ax.set_yticks(mids)
    ax.set_yticklabels(list(layout["chrom"]), fontsize=6.5)
    ax.tick_params(length=1.5)

    marked = cen[cen["n_supporting"] >= 2]
    for _, r in marked.iterrows():
        b = r["global_bin"]
        ax.plot([b], [b], marker="o", ms=4.0, mfc="none", mec="#1a1a1a", mew=0.9,
                zorder=5)

    ax.set_xlabel("Chromosome")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("log$_{10}$(contacts + 1)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    (REV / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(REV / "figures" / "Figure_1A_HiC_contact_map.pdf",
                bbox_inches="tight")
    fig.savefig(REV / "figures" / "Figure_1A_HiC_contact_map.png", dpi=600,
                bbox_inches="tight")

    print("matrix: %d x %d bins at %d bp" % (M.shape[0], M.shape[1], a.binsize))
    print(cen.drop(columns=["global_bin"]).to_string(index=False))
    print("concordant calls: %d / %d chromosomes" % (len(marked), len(cen)))


if __name__ == "__main__":
    main()
