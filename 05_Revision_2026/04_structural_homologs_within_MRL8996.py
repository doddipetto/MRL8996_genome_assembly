#!/usr/bin/env python3
"""
04_structural_homologs_within_MRL8996.py
========================================
Asks whether the sequence-unique MRL8996 effectors are also structurally novel,
using the existing all-vs-all TM-align comparison of the 482 AlphaFold2 models
of the MRL8996 effectorome.

For every effector, the best structural partner (excluding self) is recorded,
together with whether that partner is itself MRL8996-unique or shared with
Fol4287/Fo47. Two thresholds are reported: TM >= 0.5 (same fold, standard
criterion) and TM >= 0.4 (permissive, likely same fold).

Inputs (read-only):
    .../Proteins/tm_align_results_all_vs_all.tsv          Protein_A, Protein_B, Max_TM_Score
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv   (script 03)

Outputs:
    Revision_2026/tables/unique_vs_MRL_structural_homologs.tsv    per-effector best hits
    Revision_2026/tables/unique_vs_MRL_structural_summary.tsv     counts per threshold
    Revision_2026/figures/Fig_unique_intra_structural.pdf / .png  (600 dpi)

Usage: python3 04_structural_homologs_within_MRL8996.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
TM = (ADODDI / "MRL8996_genome_paper" / "JGI_reference" / "FoxMRL8996" / "Mycocosm"
      / "Annotation" / "Filtered_Models___best__" / "Proteins"
      / "tm_align_results_all_vs_all.tsv")

FOCAL = "#1f5c99"
OTHER = "#b0b7bf"
FOLD_TM = 0.5
PERMISSIVE_TM = 0.4


def pdb_to_jgi(name):
    """jgi_FoxMRL8996_10088_CE10087_16778_unrelaxed_rank_001_... -> 10088"""
    parts = name.split("_")
    return int(parts[2])


def main():
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    (REV / "figures").mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(REV / "tables" / "Table_MRL8996_unique_effectors.tsv", sep="\t")
    is_unique = dict(zip(ann["JGI_ID"], ann["Unique_to_MRL8996"]))

    tm = pd.read_csv(TM, sep="\t")
    tm["A"] = tm["Protein_A"].map(pdb_to_jgi)
    tm["B"] = tm["Protein_B"].map(pdb_to_jgi)
    tm = tm[tm["A"] != tm["B"]]

    # The all-vs-all run covers the whole modelled secretome (1,419 models);
    # restrict both members of every pair to the 513-protein final effectorome
    # (482 of which have an AlphaFold2 model).
    final = set(ann["JGI_ID"])
    tm = tm[tm["A"].isin(final) & tm["B"].isin(final)]
    print("final effectors with a model in the TM-align run: %d of %d"
          % (len(set(tm["A"]) | set(tm["B"])), len(final)))

    # symmetrise: every ordered pair appears once in the TM-align output
    pairs = pd.concat([
        tm[["A", "B", "Max_TM_Score"]].rename(columns={"A": "query", "B": "target"}),
        tm[["B", "A", "Max_TM_Score"]].rename(columns={"B": "query", "A": "target"}),
    ], ignore_index=True)
    pairs["target_unique"] = pairs["target"].map(is_unique)

    best = (pairs.sort_values("Max_TM_Score", ascending=False)
                 .drop_duplicates("query")
                 .rename(columns={"Max_TM_Score": "best_TM",
                                  "target": "best_partner",
                                  "target_unique": "best_partner_is_unique"}))

    def counts(g, thr):
        return pd.Series({
            "n_effectors": len(g),
            "n_with_hit": int((g["best_TM"] >= thr).sum()),
            "pct_with_hit": round(100 * (g["best_TM"] >= thr).mean(), 1),
            "n_hit_is_unique": int(((g["best_TM"] >= thr)
                                    & g["best_partner_is_unique"]).sum()),
            "median_best_TM": round(float(g["best_TM"].median()), 3),
        })

    best["query_unique"] = best["query"].map(is_unique)
    # number of structural partners above each threshold
    for thr, tag in ((FOLD_TM, "n_partners_TM50"), (PERMISSIVE_TM, "n_partners_TM40")):
        n = pairs[pairs["Max_TM_Score"] >= thr].groupby("query").size()
        best[tag] = best["query"].map(n).fillna(0).astype(int)

    out = best[["query", "query_unique", "best_TM", "best_partner",
                "best_partner_is_unique", "n_partners_TM50", "n_partners_TM40"]]
    out = out.rename(columns={"query": "JGI_ID", "query_unique": "Unique_to_MRL8996"})
    out.sort_values(["Unique_to_MRL8996", "best_TM"],
                    ascending=[False, False]).to_csv(
        REV / "tables" / "unique_vs_MRL_structural_homologs.tsv", sep="\t", index=False)

    rows = []
    for thr in (FOLD_TM, PERMISSIVE_TM):
        for label, g in (("MRL8996-unique", best[best["query_unique"]]),
                         ("shared", best[~best["query_unique"]])):
            r = counts(g, thr)
            r["Threshold_TM"] = thr
            r["Set"] = label
            rows.append(r)
    summary = pd.DataFrame(rows)[["Set", "Threshold_TM", "n_effectors", "n_with_hit",
                                  "pct_with_hit", "n_hit_is_unique", "median_best_TM"]]
    summary.to_csv(REV / "tables" / "unique_vs_MRL_structural_summary.tsv",
                   sep="\t", index=False)

    u = best.loc[best["query_unique"], "best_TM"].to_numpy()
    s = best.loc[~best["query_unique"], "best_TM"].to_numpy()
    stat, p = mannwhitneyu(u, s, alternative="two-sided")

    # ---------------- figure ----------------
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    bins = np.arange(0.15, 1.001, 0.025)
    ax.hist(s, bins=bins, color=OTHER, label="shared (n=%d)" % len(s), density=True)
    ax.hist(u, bins=bins, histtype="step", color=FOCAL, linewidth=1.8,
            label="MRL8996-unique (n=%d)" % len(u), density=True)
    ax.axvline(FOLD_TM, color="0.25", linestyle="--", linewidth=1)
    ax.text(FOLD_TM, ax.get_ylim()[1] * 0.06, " TM = 0.5\n same fold",
            fontsize=7, va="bottom", color="0.25")
    ax.set_xlabel("Best intra-effectorome TM-score")
    ax.set_ylabel("Density")
    ax.set_title("Most sequence-unique effectors share a fold\n"
                 "with another MRL8996 effector",
                 fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    thrs = [FOLD_TM, PERMISSIVE_TM]
    width = 0.36
    xs = np.arange(len(thrs))
    for off, (label, g, col) in enumerate((
            ("shared", best[~best["query_unique"]], OTHER),
            ("MRL8996-unique", best[best["query_unique"]], FOCAL))):
        vals = [100 * (g["best_TM"] >= t).mean() for t in thrs]
        ax.bar(xs + (off - 0.5) * width, vals, width, color=col, label=label)
        for x, v in zip(xs + (off - 0.5) * width, vals):
            ax.text(x, v, "%.0f%%" % v, ha="center", va="bottom", fontsize=7,
                    color="0.25")
    ax.set_xticks(xs)
    ax.set_xticklabels(["TM \u2265 0.5", "TM \u2265 0.4"])
    ax.set_ylabel("Effectors with a structural partner (%)")
    ax.set_ylim(0, 108)
    ax.set_title("Fold-level coverage (colours as in a)", fontsize=9, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(REV / "figures" / "Fig_unique_intra_structural.pdf", bbox_inches="tight")
    fig.savefig(REV / "figures" / "Fig_unique_intra_structural.png", dpi=600,
                bbox_inches="tight")

    print(summary.to_string(index=False))
    print("\nMann-Whitney U on best_TM (unique vs shared): U=%.0f  p=%.4g" % (stat, p))
    print("median best TM  unique=%.3f  shared=%.3f"
          % (float(np.median(u)), float(np.median(s))))


if __name__ == "__main__":
    main()
