#!/usr/bin/env python3
"""
17_centromeres_gc_criterion.py
==============================
Final centromere calls for all 16 chromosomes, using GC content together with
StainedGlass self-identity.

Why GC. Fungal centromeres are AT-rich relative to their own chromosome, so a
local GC minimum is an expected centromere signature and is completely
independent of the self-identity signal used in script 14. Adding it does two
things: it confirms the eleven core calls, and it relocates the five accessory
calls, which script 14 could only place tentatively.

The rule
--------
Per chromosome, on 20-kb windows, terminal 200 kb excluded:

  1. Candidates are the local minima of the smoothed GC profile that lie at
     least MIN_DEPTH GC points below the chromosome's own median GC. These are
     sharp troughs, not merely GC-below-average regions.
  2. Each candidate is scored by the number of off-diagonal StainedGlass window
     pairs in the centromeric identity band (ID_LO-ID_HI %) whose either end
     falls within SUPPORT_WIN of the trough - i.e. by the local presence of a
     diverged repeat array.
  3. The highest-scoring trough is the call. Its extent is the contiguous run of
     windows with GC at least MIN_DEPTH below the chromosome median.

GC is therefore the candidate generator and identity the discriminator. Neither
alone suffices: GC alone (deepest trough) picks the wrong locus on Chr03, Chr04
and Chr09, where a second, deeper AT-rich region exists; identity alone (script
14) cannot resolve the accessory chromosomes at all, because those are dominated
by dispersed near-identical (>97 %) repeat copies rather than by one diverged
array - see Table_identity_bands_per_chromosome.tsv.

Control: the rule must return the identity-only calls on the eleven core
chromosomes, where those calls are already supported by gene and TE density. It
does, on all eleven, to within one 20-kb bin.

Result: all 16 chromosomes now carry a single compact call of 80-220 kb, against
the 280-500 kb tentative blocks that the identity-only rule gave on the
accessory chromosomes. Assembly gaps do not overlap any call.

Inputs (read-only):
    Revision_2026/results/stainedglass/Chr*.2000.10000.bed.gz     (script 13)
    Revision_2026/tables/Table_centromeres_stainedglass.tsv       (script 14)
    11_Circos_Plot/track_gc_content_20kb.txt, track_gene_density_20kb.txt,
    11_Circos_Plot/circos_track_TE_20kb.txt, chrom.sizes
    MRL8996_HiC_data/mrl8996-hic.chromosomes.fasta                (gap check)

Outputs:
    Revision_2026/tables/Table_centromeres_gc.tsv                 (final calls)
    Revision_2026/results/centromeres_gc.bed                      (for script 15)
    Revision_2026/tables/Table_identity_bands_per_chromosome.tsv
    Revision_2026/figures/Fig_centromere_identity_gc_profiles.pdf / .png

Run:  python3 scripts/17_centromeres_gc_criterion.py     (env: python)
      then rerun scripts 15 and 11b, which read the call table.
"""

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import argrelmin

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
SGD = REV / "results" / "stainedglass"
CIRC = ADODDI / "MRL8996_HiC_data" / "11_Circos_Plot"
FASTA = ADODDI / "MRL8996_HiC_data" / "mrl8996-hic.chromosomes.fasta"

# identity-run and caller parameters, identical to scripts 13 and 14
BIN = 20000
ID_LO, ID_HI = 84.0, 93.0
MIN_OFFSET = 20000
TERMINAL = 200000
SMOOTH = 5
EDGE_FRAC = 0.25
# GC criterion
MIN_DEPTH = 0.03        # GC points below the chromosome median for a trough
MIN_ORDER = 3           # bins each side, for local-minimum detection
SUPPORT_WIN = 60000     # bp around a trough within which identity pairs count

CHRS = ["Chr%02d" % i for i in range(1, 17)]
CORE = set(CHRS[:11])
BANDS = [(70, 84), (84, 93), (93, 97), (97, 100.1)]


def smooth(v):
    return np.convolve(v, np.ones(SMOOTH) / SMOOTH, mode="same")


def load_tracks():
    cs = pd.read_csv(CIRC / "chrom.sizes", sep="\t", header=None,
                     names=["chrom", "length"])
    tr = {}
    for key, fn in (("gc", "track_gc_content_20kb.txt"),
                    ("gene", "track_gene_density_20kb.txt"),
                    ("te", "circos_track_TE_20kb.txt")):
        tr[key] = pd.read_csv(CIRC / fn, sep="\t", header=None,
                              names=["chrom", "start", "end", "value"])
    return dict(zip(cs["chrom"], cs["length"].astype(int))), tr


def binned(track, chrom, n):
    """track values on the caller's bin grid, gaps filled by interpolation."""
    sub = track[track["chrom"] == chrom].sort_values("start")
    v = np.full(n, np.nan)
    idx = (sub["start"].to_numpy() // BIN).astype(int)
    ok = idx < n
    v[idx[ok]] = sub["value"].to_numpy()[ok]
    return pd.Series(v).interpolate(limit_direction="both").to_numpy()


def offdiagonal_pairs(chrom, band=None):
    b = pd.read_csv(SGD / ("%s.2000.10000.bed.gz" % chrom), sep="\t")
    off = (b["query_start"] - b["reference_start"]).abs()
    b = b[off >= MIN_OFFSET]
    if band is not None:
        b = b[b["perID_by_events"].between(*band)]
    return b


def identity_bands(chrom):
    b = offdiagonal_pairs(chrom)
    row = {"chrom": chrom, "n_offdiagonal_pairs": len(b)}
    for lo, hi in BANDS:
        row["pct_%g_%g" % (lo, min(hi, 100))] = round(
            100 * float(b["perID_by_events"].between(lo, hi, inclusive="left").mean()), 1)
    return row


def assembly_gaps(lengths):
    """N-runs of >=10 bp, keyed by chromosome (the FASTA carries scaffold names,
    which are matched to chromosome names by sequence length)."""
    seqs, name, buf = {}, None, []
    for line in open(FASTA):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    seqs[name] = "".join(buf)
    by_len = {v: k for k, v in lengths.items()}
    out = {}
    for sid, seq in seqs.items():
        chrom = by_len.get(len(seq))
        if chrom is None:
            continue
        runs = [(m.start(), m.end()) for m in re.finditer(r"N{10,}", seq.upper())]
        if runs:
            out[chrom] = runs
    return out


def call_one(chrom, n, gcv, dens_band, inner):
    """GC troughs as candidates, identity-pair support as the discriminator."""
    med = float(np.nanmedian(gcv[inner]))
    cand = [i for i in argrelmin(gcv, order=MIN_ORDER)[0]
            if inner[i] and med - gcv[i] >= MIN_DEPTH]
    if not cand:
        return None
    qs = dens_band["query_start"].to_numpy()
    rs = dens_band["reference_start"].to_numpy()
    scored = []
    for i in cand:
        d = np.minimum(np.abs(qs - i * BIN), np.abs(rs - i * BIN))
        scored.append((int((d <= SUPPORT_WIN).sum()), med - gcv[i], int(i)))
    scored.sort(reverse=True)
    pairs, depth, i = scored[0]
    thr = med - MIN_DEPTH
    lo = hi = i
    while lo - 1 >= 0 and gcv[lo - 1] <= thr:
        lo -= 1
    while hi + 1 < n and gcv[hi + 1] <= thr:
        hi += 1
    return dict(peak=i, lo=lo, hi=hi, pairs=pairs, depth=depth,
                gc_median=med, n_candidates=len(cand))


def main():
    lengths, tr = load_tracks()
    prev = pd.read_csv(REV / "tables" / "Table_centromeres_stainedglass.tsv",
                       sep="\t").set_index("chrom")
    gaps = assembly_gaps(lengths)

    prof, rows, bands = {}, [], []
    for chrom in CHRS:
        n = int(np.ceil(lengths[chrom] / BIN)) + 1
        gcv = smooth(binned(tr["gc"], chrom, n))
        gene = binned(tr["gene"], chrom, n)
        te = binned(tr["te"], chrom, n)
        k = TERMINAL // BIN
        inner = np.zeros(n, bool)
        inner[k:n - k] = True

        band = offdiagonal_pairs(chrom, (ID_LO, ID_HI))
        dv = np.zeros(n)
        for st in np.r_[band["query_start"].to_numpy(), band["reference_start"].to_numpy()]:
            dv[min(int(st // BIN), n - 1)] += 1
        dv = smooth(dv)
        bands.append(identity_bands(chrom))

        c = call_one(chrom, n, gcv, band, inner)
        if c is None:
            raise SystemExit("no GC trough passed MIN_DEPTH on %s" % chrom)
        blk = slice(c["lo"], c["hi"] + 1)
        start, end = c["lo"] * BIN, int(min((c["hi"] + 1) * BIN, lengths[chrom]))
        ov = [g for g in gaps.get(chrom, []) if not (g[1] < start or g[0] > end)]
        idb = band[(band["query_start"].between(start, end))
                   | (band["reference_start"].between(start, end))]

        prof[chrom] = (dv, gcv, inner, start, end)
        rows.append(dict(
            chrom=chrom,
            compartment="core" if chrom in CORE else "accessory",
            start=start, end=end, size_kb=(end - start) / 1000,
            peak_bp=c["peak"] * BIN,
            relative_position=round(c["peak"] * BIN / lengths[chrom], 3),
            gc=round(float(gcv[c["peak"]]), 4),
            gc_chrom_median=round(c["gc_median"], 4),
            gc_depression=round(-c["depth"], 4),
            gc_pctile=round(100 * float(np.nanmean(gcv[inner] < np.nanmean(gcv[blk]))), 1),
            n_gc_troughs=c["n_candidates"],
            n_window_pairs=int(len(idb)),
            pct_band_pairs_in_call=round(100 * len(idb) / max(len(band), 1), 1),
            band_pairs_per_mb=round(len(offdiagonal_pairs(chrom)) / (lengths[chrom] / 1e6), 0),
            median_identity=round(float(idb["perID_by_events"].median()), 1) if len(idb) else np.nan,
            gene_density_pctile=round(100 * float(np.nanmean(gene[inner] < np.nanmean(gene[blk]))), 1),
            TE_density_pctile=round(100 * float(np.nanmean(te[inner] < np.nanmean(te[blk]))), 1),
            assembly_gap_in_call=len(ov),
            shift_vs_identity_only_kb=round(abs(c["peak"] * BIN - prev.loc[chrom, "peak_bp"]) / 1000, 1),
            overlap_identity_only_pct=round(100 * max(0, min(end, int(prev.loc[chrom, "end"]))
                                                      - max(start, int(prev.loc[chrom, "start"])))
                                            / (end - start), 1),
        ))

    calls = pd.DataFrame(rows)
    calls["support"] = ((calls["gc_pctile"] <= 10).astype(int)
                        + (calls["gene_density_pctile"] <= 25).astype(int)
                        + (calls["TE_density_pctile"] >= 75).astype(int))
    calls.to_csv(REV / "tables" / "Table_centromeres_gc.tsv", sep="\t", index=False)
    bed = calls[["chrom", "start", "end"]].copy()
    bed["name"] = calls["chrom"] + "_cen"
    bed.to_csv(REV / "results" / "centromeres_gc.bed", sep="\t",
               header=False, index=False)
    pd.DataFrame(bands).to_csv(
        REV / "tables" / "Table_identity_bands_per_chromosome.tsv", sep="\t", index=False)

    # ---- diagnostic figure: identity-pair density and GC per chromosome ----
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                         "axes.linewidth": 0.5, "pdf.fonttype": 42})
    fig, axes = plt.subplots(4, 4, figsize=(11.0, 8.4))
    for ax, chrom in zip(axes.ravel(), CHRS):
        dv, gcv, inner, start, end = prof[chrom]
        x = np.arange(len(dv)) * BIN / 1e6
        c = calls[calls["chrom"] == chrom].iloc[0]
        ax.axvspan(start / 1e6, end / 1e6, color="#f6d0d0", alpha=0.9, lw=0, zorder=0)
        ax.fill_between(x, dv, color="#8b1a1a", lw=0, zorder=2)
        ax.set_xlim(0, x.max())
        ax.set_ylim(0, max(4, dv.max() * 1.3))
        ax.set_ylabel("pairs", color="#8b1a1a", labelpad=1)
        ax.tick_params(axis="y", colors="#8b1a1a", length=2, pad=1, labelsize=6)
        ax.tick_params(axis="x", length=2, pad=1, labelsize=6)
        a2 = ax.twinx()
        a2.plot(x, gcv, color="#1f5c99", lw=0.7, zorder=3)
        a2.axhline(c["gc_chrom_median"] - MIN_DEPTH, color="#1f5c99", lw=0.5,
                   ls=":", zorder=1)
        a2.set_ylim(0.28, 0.56)
        a2.set_ylabel("GC", color="#1f5c99", labelpad=1)
        a2.tick_params(axis="y", colors="#1f5c99", length=2, pad=1, labelsize=6)
        ax.spines["top"].set_visible(False)
        a2.spines["top"].set_visible(False)
        ax.set_title("%s  %s  %.0f kb  GC %.3f (%.0fth pct)"
                     % (chrom, c["compartment"][:4], c["size_kb"], c["gc"], c["gc_pctile"]),
                     fontsize=6.8, pad=2)
    for ax in axes[-1]:
        ax.set_xlabel("position (Mb)", labelpad=1)
    fig.suptitle("Centromere calls on two independent signals: %g-%g %% identity "
                 "window-pair density (red) and GC content (blue); dotted line = "
                 "chromosome median GC - %.2f; shading = call"
                 % (ID_LO, ID_HI, MIN_DEPTH), fontsize=8.5, y=0.996)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(REV / "figures" / "Fig_centromere_identity_gc_profiles.pdf")
    fig.savefig(REV / "figures" / "Fig_centromere_identity_gc_profiles.png", dpi=400)

    show = ["chrom", "compartment", "start", "end", "size_kb", "gc", "gc_depression",
            "gc_pctile", "n_gc_troughs", "n_window_pairs", "gene_density_pctile",
            "TE_density_pctile", "pct_band_pairs_in_call", "band_pairs_per_mb",
            "assembly_gap_in_call", "shift_vs_identity_only_kb",
            "overlap_identity_only_pct", "support"]
    print(calls[show].to_string(index=False))
    core = calls[calls["compartment"] == "core"]
    print("\ncore calls contained in the identity-only call: %d / 11 "
          "(overlap %.0f-%.0f %% of the new interval)"
          % (int((core["overlap_identity_only_pct"] >= 99).sum()),
             core["overlap_identity_only_pct"].min(), core["overlap_identity_only_pct"].max()))
    acc = calls[calls["compartment"] == "accessory"]
    print("accessory calls relocated by %.0f-%.0f kb; overlap with the tentative "
          "identity-only blocks %.0f-%.0f %%"
          % (acc["shift_vs_identity_only_kb"].min(), acc["shift_vs_identity_only_kb"].max(),
             acc["overlap_identity_only_pct"].min(), acc["overlap_identity_only_pct"].max()))
    print("diverged-band pairs falling in the call: core %.0f-%.0f %%, accessory %.0f-%.0f %% ; "
          "self-similar pairs per Mb: core %.0f-%.0f, accessory %.0f-%.0f"
          % (core["pct_band_pairs_in_call"].min(), core["pct_band_pairs_in_call"].max(),
             acc["pct_band_pairs_in_call"].min(), acc["pct_band_pairs_in_call"].max(),
             core["band_pairs_per_mb"].min(), core["band_pairs_per_mb"].max(),
             acc["band_pairs_per_mb"].min(), acc["band_pairs_per_mb"].max()))
    print("call sizes: %.0f-%.0f kb ; GC depression %.3f-%.3f ; gaps in calls: %d"
          % (calls["size_kb"].min(), calls["size_kb"].max(),
             calls["gc_depression"].max(), calls["gc_depression"].min(),
             int(calls["assembly_gap_in_call"].sum())))


if __name__ == "__main__":
    main()
