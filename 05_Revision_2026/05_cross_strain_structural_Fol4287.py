#!/usr/bin/env python3
"""
05_cross_strain_structural_Fol4287.py
=====================================
Structural comparison of the MRL8996-unique effectors against the Fol4287
effectorome: every MRL8996-unique AlphaFold2 model is aligned with TM-align to
every Fol4287 effector model, and the best cross-strain hit per query is kept.

This answers the second half of the reviewer's question: effectors with no
sequence homolog in the reference strains may still adopt a fold that is
present there. Effectors whose best cross-strain TM-score stays below 0.5 are
structurally, not only sequence-wise, particular to MRL8996.

HEAVY STEP - run this yourself; it is CPU-bound, not sandbox-friendly.
    ~122 queries x ~500 targets = ~60,000 TM-align calls.
    With --threads 20 on a 24-core machine expect roughly 10-20 minutes.

Inputs (read-only):
    .../Dali_Effectorome_Run/PDB_clean/*.pdb                MRL8996 models
    Fol4287_genome_paper/colabfold_results/*rank_001*.pdb   Fol4287 models
    F012_HiC_data/final_results/Effectorome_Final_Results/Fol4287_final_effectors.fasta
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv (script 03)

Outputs:
    Revision_2026/results/tmalign_MRL8996unique_vs_Fol4287.tsv   all pairs
    Revision_2026/tables/cross_strain_best_hits.tsv              best hit per query
    Revision_2026/figures/Fig_cross_strain_structural.pdf / .png (600 dpi)

Usage:
    python3 05_cross_strain_structural_Fol4287.py --threads 20
    python3 05_cross_strain_structural_Fol4287.py --threads 20 --plot-only
"""

import argparse
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
MRL_PDB = (ADODDI / "MRL8996_genome_paper" / "JGI_reference" / "FoxMRL8996" / "Mycocosm"
           / "Annotation" / "Filtered_Models___best__" / "Proteins"
           / "Dali_Effectorome_Run" / "PDB_clean")
FOL_PDB = ADODDI / "Fol4287_genome_paper" / "colabfold_results"
FOL_FASTA = (ADODDI / "F012_HiC_data" / "final_results" / "Effectorome_Final_Results"
             / "Fol4287_final_effectors.fasta")

PAIRS_OUT = REV / "results" / "tmalign_MRL8996unique_vs_Fol4287.tsv"
BEST_OUT = REV / "tables" / "cross_strain_best_hits.tsv"

FOLD_TM = 0.5
PERMISSIVE_TM = 0.4
UNIQUE_C = "#1f5c99"
OTHER_C = "#b0b7bf"

TM_RE = re.compile(r"TM-score=\s*([0-9.]+)")


def fol_accession(name):
    """XP_018231541.1_hypothetical_protein_... -> XP_018231541.1"""
    m = re.match(r"([A-Z]P_\d+\.\d+)", name)
    return m.group(1) if m else None


def collect_targets():
    """One rank_001 model per Fol4287 final effector (unrelaxed preferred)."""
    wanted = set()
    with open(FOL_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                wanted.add(line[1:].split()[0])
    best = {}
    for p in FOL_PDB.glob("*rank_001*.pdb"):
        acc = fol_accession(p.name)
        if acc is None or acc not in wanted:
            continue
        # prefer the unrelaxed model for consistency with the MRL8996 set
        if acc not in best or ("unrelaxed" in p.name and "unrelaxed" not in best[acc].name):
            best[acc] = p
    missing = len(wanted) - len(best)
    print("Fol4287 effectors: %d ; with a model: %d ; missing: %d"
          % (len(wanted), len(best), missing), file=sys.stderr)
    return best


def collect_queries():
    ann = pd.read_csv(REV / "tables" / "Table_MRL8996_unique_effectors.tsv", sep="\t")
    uniq_ids = set(ann.loc[ann["Unique_to_MRL8996"], "JGI_ID"].astype(int))
    # PDB_clean holds the DALI-cleaned models named by Dali_ID (P001.pdb),
    # so the JGI protein IDs come from the mapping dictionary.
    mapping = pd.read_csv(MRL_PDB.parent / "dali_mapping_dict.tsv", sep="\t")
    dali2jgi = dict(zip(mapping["Dali_ID"], mapping["JGI_Protein_ID"].astype(int)))
    out = {}
    for p in sorted(MRL_PDB.glob("P*.pdb")):
        jgi = dali2jgi.get(p.stem)
        if jgi is not None and jgi in uniq_ids:
            out[jgi] = p
    print("MRL8996-unique effectors with a model: %d of %d"
          % (len(out), len(uniq_ids)), file=sys.stderr)
    return out


def tm_pair(args):
    q, qp, t, tp = args
    try:
        r = subprocess.run(["TMalign", str(qp), str(tp)], capture_output=True,
                           text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return (q, t, float("nan"))
    scores = [float(x) for x in TM_RE.findall(r.stdout)]
    # TM-align prints TM-scores normalised by each chain; take the larger,
    # matching the Max_TM_Score convention of the intra-MRL8996 table.
    return (q, t, max(scores) if scores else float("nan"))


def run_alignments(threads):
    queries = collect_queries()
    targets = collect_targets()
    jobs = [(q, qp, t, tp) for q, qp in queries.items() for t, tp in targets.items()]
    print("TM-align calls: %d" % len(jobs), file=sys.stderr)
    PAIRS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(PAIRS_OUT, "w") as out, ProcessPoolExecutor(max_workers=threads) as ex:
        out.write("MRL8996_JGI_ID\tFol4287_accession\tTM_score\n")
        for i, (q, t, s) in enumerate(ex.map(tm_pair, jobs, chunksize=64), 1):
            out.write("%s\t%s\t%.5f\n" % (q, t, s))
            if i % 5000 == 0:
                print("  %d / %d" % (i, len(jobs)), file=sys.stderr, flush=True)
    print("wrote %s" % PAIRS_OUT, file=sys.stderr)


def summarise_and_plot():
    pairs = pd.read_csv(PAIRS_OUT, sep="\t")
    best = (pairs.sort_values("TM_score", ascending=False)
                 .drop_duplicates("MRL8996_JGI_ID")
                 .rename(columns={"TM_score": "best_cross_strain_TM",
                                  "Fol4287_accession": "best_Fol4287_hit"}))
    intra = pd.read_csv(REV / "tables" / "unique_vs_MRL_structural_homologs.tsv", sep="\t")
    best = best.merge(intra[["JGI_ID", "best_TM"]].rename(
        columns={"JGI_ID": "MRL8996_JGI_ID", "best_TM": "best_intra_MRL8996_TM"}),
        on="MRL8996_JGI_ID", how="left")
    best["fold_present_in_Fol4287"] = best["best_cross_strain_TM"] >= FOLD_TM
    best.sort_values("best_cross_strain_TM").to_csv(BEST_OUT, sep="\t", index=False)

    n = len(best)
    n_fold = int(best["fold_present_in_Fol4287"].sum())
    n_novel = n - n_fold

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))

    ax = axes[0]
    bins = np.arange(0.15, 1.001, 0.025)
    ax.hist(best["best_cross_strain_TM"].dropna(), bins=bins, color=UNIQUE_C)
    ax.axvline(FOLD_TM, color="0.25", linestyle="--", linewidth=1)
    ax.text(FOLD_TM - 0.01, ax.get_ylim()[1] * 0.98, "TM = 0.5 ", fontsize=7,
            va="top", ha="right", color="0.25")
    ax.set_xlabel("Best TM-score against the Fol4287 effectorome")
    ax.set_ylabel("MRL8996-unique effectors")
    ax.set_title("%d of %d sequence-unique effectors adopt a fold\nalso present "
                 "in Fol4287" % (n_fold, n), fontsize=9, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    ax.scatter(best["best_intra_MRL8996_TM"], best["best_cross_strain_TM"],
               s=14, c=UNIQUE_C, alpha=0.75, edgecolors="none")
    ax.axhline(FOLD_TM, color="0.6", linestyle="--", linewidth=0.8)
    ax.axvline(FOLD_TM, color="0.6", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Best TM-score within MRL8996")
    ax.set_ylabel("Best TM-score vs Fol4287")
    n_both = int(((~best["fold_present_in_Fol4287"])
                  & (best["best_intra_MRL8996_TM"] < FOLD_TM)).sum())
    ax.set_title("%d effectors have no fold match in either\nset (lower-left "
                 "quadrant)" % n_both, fontsize=9, loc="left")
    ax.set_xlim(0.2, 1.02)
    ax.set_ylim(0.2, 1.02)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(REV / "figures" / "Fig_cross_strain_structural.pdf", bbox_inches="tight")
    fig.savefig(REV / "figures" / "Fig_cross_strain_structural.png", dpi=600,
                bbox_inches="tight")

    print("unique effectors compared: %d" % n)
    print("fold present in Fol4287 (TM >= %.1f): %d" % (FOLD_TM, n_fold))
    print("no Fol4287 fold match: %d" % n_novel)
    print("median best cross-strain TM: %.3f"
          % float(best["best_cross_strain_TM"].median()))
    print("also without an intra-MRL8996 partner at TM >= %.1f: %d"
          % (FOLD_TM, int(((~best["fold_present_in_Fol4287"])
                           & (best["best_intra_MRL8996_TM"] < FOLD_TM)).sum())))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true",
                    help="skip the TM-align run and re-summarise existing pairs")
    a = ap.parse_args()
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    (REV / "figures").mkdir(parents=True, exist_ok=True)
    if not a.plot_only:
        run_alignments(a.threads)
    summarise_and_plot()


if __name__ == "__main__":
    main()
