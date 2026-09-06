#!/usr/bin/env python3
"""
09_accessory_chromosome_content.py
==================================
Per-chromosome gene, TE, BGC and effector content of the MRL8996 assembly,
split into the core (Chr01-Chr11) and accessory (Chr12-Chr16) compartments.
Written for Reviewer 2's request for a description of the genes and TEs on the
accessory chromosomes, and for the auN values missing from Table 2.

Metrics per chromosome:
    length, gene count and genes/Mb, mean exons per gene, CDS bp and coding
    fraction, repeat-masked bp and % by TE class (RepeatMasker output of the
    EDTA run, overlapping hits merged per class), count of intact TE elements
    by superfamily (EDTA TEanno GFF3), BGC count (antiSMASH region files),
    effector count and effectors/Mb, MRL8996-unique effector count.

auN (area under the Nx curve) is computed for both assembly versions from the
FASTA files, giving the length-weighted contiguity statistic Reviewer 1 asked
for alongside N50.

Inputs (read-only):
    MRL8996_HiC_data/08_circos_features/MRL8996_chromosome_level.gff3
    MRL8996_HiC_data/08_circos_features/chrom.sizes
    MRL8996_HiC_data/09_TE_annotation/...EDTA.anno/...EDTA.RM.out
    MRL8996_HiC_data/09_TE_annotation/...EDTA.TEanno.gff3
    MRL8996_HiC_data/10_BGC_annotation/MRL8996_BGC/Chr*.region*.gbk
    MRL8996_HiC_data/FoxMRL8996_AssemblyScaffolds.fasta   (JGI scaffolds)
    MRL8996_HiC_data/MRL8996_chromosome_level.fasta       (this study)
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv (script 03)

Outputs:
    Revision_2026/tables/Table_per_chromosome_content.tsv
    Revision_2026/tables/Table_TE_classes_core_vs_accessory.tsv
    Revision_2026/tables/assembly_auN.tsv
    Revision_2026/figures/Fig_core_vs_accessory_content.pdf / .png (600 dpi)

Usage: python3 09_accessory_chromosome_content.py
"""

import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"

GFF = HIC / "08_circos_features" / "MRL8996_chromosome_level.gff3"
SIZES = HIC / "08_circos_features" / "chrom.sizes"
TE_DIR = HIC / "09_TE_annotation"
RM_OUT = (TE_DIR / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.anno"
          / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.RM.out")
TEANNO = TE_DIR / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.TEanno.gff3"
BGC_DIR = HIC / "10_BGC_annotation" / "MRL8996_BGC"
ASSEMBLIES = {"JGI_scaffolds_v1": HIC / "FoxMRL8996_AssemblyScaffolds.fasta",
              "MRL8996_chromosome_level": HIC / "MRL8996_chromosome_level.fasta"}

ACCESSORY = {"Chr12", "Chr13", "Chr14", "Chr15", "Chr16"}
CORE_C = "#8c9196"
ACC_C = "#b2182b"


def merged_bp(intervals):
    """Total bp covered by a list of (start, end) 1-based inclusive intervals."""
    if not intervals:
        return 0
    iv = sorted(intervals)
    total, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    return total + ce - cs + 1


def read_sizes():
    s = pd.read_csv(SIZES, sep="\t", header=None, names=["chrom", "length"])
    return s[s["chrom"].str.match(r"Chr\d+")].set_index("chrom")["length"].to_dict()


def gene_metrics(chroms):
    genes = defaultdict(int)
    exons = defaultdict(int)
    cds = defaultdict(list)
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in chroms:
                continue
            if f[2] == "gene":
                genes[f[0]] += 1
            elif f[2] == "exon":
                exons[f[0]] += 1
            elif f[2] == "CDS":
                cds[f[0]].append((int(f[3]), int(f[4])))
    return genes, exons, {c: merged_bp(v) for c, v in cds.items()}


def te_metrics(chroms):
    """Masked bp per (chrom, TE class) from the RepeatMasker output."""
    by_class = defaultdict(list)
    with open(RM_OUT) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[4].startswith("Chr"):
                continue
            chrom = f[4]
            if chrom not in chroms:
                continue
            cls = f[10].split("/")[0]
            by_class[(chrom, cls)].append((int(f[5]), int(f[6])))
    bp = {k: merged_bp(v) for k, v in by_class.items()}
    total = defaultdict(int)
    for (chrom, _), v in bp.items():
        total[chrom] += v
    # counts of intact elements by superfamily
    intact = defaultdict(int)
    with open(TEANNO) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[0] not in chroms:
                continue
            if f[2] in ("repeat_fragment", "target_site_duplication",
                        "long_terminal_repeat", "repeat_region"):
                continue
            intact[(f[0], f[2])] += 1
    return bp, total, intact


def bgc_counts(chroms):
    n = defaultdict(int)
    for p in BGC_DIR.glob("Chr*.region*.gbk"):
        chrom = p.name.split(".")[0]
        if chrom in chroms:
            n[chrom] += 1
    return n


def effector_metrics(chroms):
    ann = pd.read_csv(REV / "tables" / "Table_MRL8996_unique_effectors.tsv", sep="\t")
    ann = ann[ann["Chromosome"].isin(chroms)]
    total = ann.groupby("Chromosome").size().to_dict()
    uniq = ann[ann["Unique_to_MRL8996"]].groupby("Chromosome").size().to_dict()
    return total, uniq


def auN(fasta):
    lens = []
    n = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    lens = np.array(sorted(lens, reverse=True), dtype=float)
    tot = lens.sum()
    return float((lens ** 2).sum() / tot), float(tot), len(lens)


def main():
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    (REV / "figures").mkdir(parents=True, exist_ok=True)

    sizes = read_sizes()
    chroms = set(sizes)
    genes, exons, cdsbp = gene_metrics(chroms)
    te_bp, te_total, te_intact = te_metrics(chroms)
    bgc = bgc_counts(chroms)
    eff, eff_u = effector_metrics(chroms)

    rows = []
    for c in sorted(chroms):
        L = sizes[c]
        mb = L / 1e6
        rows.append({
            "chrom": c,
            "compartment": "accessory" if c in ACCESSORY else "core",
            "length_bp": L,
            "genes": genes.get(c, 0),
            "genes_per_Mb": round(genes.get(c, 0) / mb, 1),
            "exons_per_gene": round(exons.get(c, 0) / genes[c], 2) if genes.get(c) else np.nan,
            "CDS_bp": cdsbp.get(c, 0),
            "coding_pct": round(100 * cdsbp.get(c, 0) / L, 2),
            "repeat_bp": te_total.get(c, 0),
            "repeat_pct": round(100 * te_total.get(c, 0) / L, 2),
            "intact_TE_elements": sum(v for (cc, _), v in te_intact.items() if cc == c),
            "BGCs": bgc.get(c, 0),
            "effectors": eff.get(c, 0),
            "effectors_per_Mb": round(eff.get(c, 0) / mb, 2),
            "unique_effectors": eff_u.get(c, 0),
        })
    tab = pd.DataFrame(rows)
    tab.to_csv(REV / "tables" / "Table_per_chromosome_content.tsv", sep="\t", index=False)

    # TE classes, core vs accessory
    cls_rows = []
    classes = sorted({k[1] for k in te_bp})
    for cls in classes:
        r = {"TE_class": cls}
        for comp, members in (("core", chroms - ACCESSORY), ("accessory", ACCESSORY)):
            bp = sum(v for (c, cl), v in te_bp.items() if cl == cls and c in members)
            tot = sum(sizes[c] for c in members)
            r["%s_bp" % comp] = bp
            r["%s_pct" % comp] = round(100 * bp / tot, 3)
        cls_rows.append(r)
    cls_tab = pd.DataFrame(cls_rows).sort_values("accessory_pct", ascending=False)
    cls_tab.to_csv(REV / "tables" / "Table_TE_classes_core_vs_accessory.tsv",
                   sep="\t", index=False)

    # auN
    aun_rows = []
    for name, f in ASSEMBLIES.items():
        a, tot, n = auN(f)
        aun_rows.append({"assembly": name, "n_sequences": n,
                         "total_bp": int(tot), "auN_bp": round(a, 1),
                         "auN_Mb": round(a / 1e6, 3)})
    aun = pd.DataFrame(aun_rows)
    aun.to_csv(REV / "tables" / "assembly_auN.tsv", sep="\t", index=False)

    # ---------------- figure ----------------
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3))
    colours = [ACC_C if c in ACCESSORY else CORE_C for c in tab["chrom"]]
    xs = np.arange(len(tab))
    xlab = [c.replace("Chr", "") for c in tab["chrom"]]

    for ax, col, ylab, title in (
            (axes[0], "genes_per_Mb", "Genes per Mb", "Gene density"),
            (axes[1], "repeat_pct", "Repeat-masked (% of length)", "Repeat content"),
            (axes[2], "effectors_per_Mb", "Effectors per Mb", "Effector density")):
        ax.bar(xs, tab[col], color=colours, width=0.78)
        ax.set_xticks(xs)
        ax.set_xticklabels(xlab, fontsize=7)
        ax.set_xlabel("Chromosome")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9, loc="left")
        ax.spines[["top", "right"]].set_visible(False)

    core = tab[tab["compartment"] == "core"]
    acc = tab[tab["compartment"] == "accessory"]
    axes[2].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=CORE_C,
                                          label="core (Chr01-11)"),
                            plt.Rectangle((0, 0), 1, 1, color=ACC_C,
                                          label="accessory (Chr12-16)")],
                   frameon=False, fontsize=7, loc="upper right")
    fig.suptitle("Accessory chromosomes are repeat-rich and gene-poor: "
                 "%.0f vs %.0f genes/Mb and %.1f%% vs %.1f%% repeats"
                 % (core["genes_per_Mb"].mean(), acc["genes_per_Mb"].mean(),
                    core["repeat_pct"].mean(), acc["repeat_pct"].mean()),
                 fontsize=9.5, y=1.02)
    fig.tight_layout()
    fig.savefig(REV / "figures" / "Fig_core_vs_accessory_content.pdf",
                bbox_inches="tight")
    fig.savefig(REV / "figures" / "Fig_core_vs_accessory_content.png", dpi=600,
                bbox_inches="tight")

    print(tab.to_string(index=False))
    print()
    print("core vs accessory means")
    print(tab.groupby("compartment")[["genes_per_Mb", "coding_pct", "repeat_pct",
                                      "effectors_per_Mb", "exons_per_gene"]]
          .mean().round(2).to_string())
    print()
    print("accessory totals: %d genes, %d effectors (%d MRL8996-unique), %d BGCs, "
          "%.2f Mb"
          % (acc["genes"].sum(), acc["effectors"].sum(),
             acc["unique_effectors"].sum(), acc["BGCs"].sum(),
             acc["length_bp"].sum() / 1e6))
    print()
    print(aun.to_string(index=False))
    print()
    print(cls_tab.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
