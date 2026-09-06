#!/usr/bin/env python3
"""
03_characterise_unique_effectors.py
===================================
Annotates the MRL8996-unique effectors (CD-HIT clusters with no Fol4287 or Fo47
member) with genomic position, compartment (core vs accessory chromosome),
AMAPEC antimicrobial prediction, model confidence (pLDDT) and pre-existing
TM-score structural family, and tests whether the unique set is enriched on
accessory chromosomes or among antimicrobial-predicted effectors.

Inputs (read-only):
    Revision_2026/results/lists/4_MRL8996_unique.txt      (script 01)
    MRL8996_HiC_data/effector_map_data.tsv                position + AMAPEC + pLDDT
    MRL8996_genome_paper/Heatmap_effectorome/Final_Effectorome_Families_Complex.tsv
    MRL8996_genome_paper/structural_network/OCE_struct_network_nodes.tsv
    MRL8996_HiC_data/mrl8996-hic.final.bed                chromosome lengths

Outputs:
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv
    Revision_2026/tables/unique_effector_enrichment.tsv
    Revision_2026/figures/Fig_unique_effector_distribution.pdf / .png (600 dpi)

Usage: python3 03_characterise_unique_effectors.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
HIC = ADODDI / "MRL8996_HiC_data"
PAPER = ADODDI / "MRL8996_genome_paper"

# Chromosomes 12-16 are the accessory (lineage-specific) compartment of MRL8996.
CORE_CHR = {"Chr%02d" % i for i in range(1, 12)}

FOCAL = "#1f5c99"       # MRL8996-unique
OTHER = "#b0b7bf"       # shared effectors


def effector_coordinates(ids):
    """Coordinates of the final effectorome on the Hi-C assembly.

    Taken from the miniprot mapping of the JGI Filtered Models onto the
    chromosome-scale assembly (MRL8996_HiC_proteins.gff), whose Target
    attribute carries the same jgi|FoxMRL8996|<proteinId>|<model> headers as the
    effectorome FASTA. Only the top-ranked (Rank=1) mRNA mapping is kept.
    """
    want = set(ids)
    rows = {}
    with open(HIC / "MRL8996_HiC_proteins.gff") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            attrs = f[8]
            i = attrs.find("Target=jgi|FoxMRL8996|")
            if i < 0:
                continue
            target = attrs[i + len("Target=jgi|FoxMRL8996|"):].split(";")[0]
            pid = target.split("|")[0]
            if not pid.isdigit() or int(pid) not in want:
                continue
            pid = int(pid)
            score = float(f[5]) if f[5] not in (".", "") else 0.0
            if pid not in rows or score > rows[pid][3]:
                rows[pid] = (f[0], int(f[3]), int(f[4]), score, f[6])
    return pd.DataFrame(
        [{"JGI_ID": k, "Chromosome": v[0], "Start": v[1], "End": v[2],
          "Strand": v[4]} for k, v in rows.items()])


def load():
    uniq = pd.read_csv(REV / "results" / "lists" / "4_MRL8996_unique.txt", sep="\t")
    uniq["JGI_ID"] = uniq["MRL8996_effector"].str.split("_").str[0].astype(int)

    # Full final effectorome (513 IDs), not only the 482 with AF2 models.
    ids_path = (ADODDI / "MRL8996_genome_paper" / "JGI_reference" / "FoxMRL8996"
                / "Mycocosm" / "Annotation" / "Filtered_Models___best__" / "Proteins"
                / "FoxMRL8996_Final_Effectorome_IDs.txt")
    all_ids = [int(l.split("|")[2]) for l in open(ids_path) if l.strip()]

    coords = effector_coordinates(all_ids)

    amapec = pd.read_csv(HIC / "effector_map_data.tsv", sep="\t")
    amapec = amapec.rename(columns={
        "Protein ID": "JGI_ID",
        "Probability of antimicrobial activity": "AM_probability",
        "Prediction": "AMAPEC",
    })[["JGI_ID", "pLDDT", "AM_probability", "AMAPEC"]]
    amapec["JGI_ID"] = amapec["JGI_ID"].astype(int)

    emap = pd.DataFrame({"JGI_ID": all_ids}).merge(
        coords, on="JGI_ID", how="left").merge(amapec, on="JGI_ID", how="left")

    fam = pd.read_csv(PAPER / "Heatmap_effectorome" /
                      "Final_Effectorome_Families_Complex.tsv", sep="\t")
    fam = fam[["JGI_ID", "Family_ID", "Family_Rank"]].copy()
    fam["JGI_ID"] = fam["JGI_ID"].astype(int)

    net = pd.read_csv(PAPER / "structural_network" /
                      "OCE_struct_network_nodes.tsv", sep="\t")
    net = net[["JGI_Protein_ID", "family", "degree"]].rename(
        columns={"JGI_Protein_ID": "JGI_ID", "family": "Network_community",
                 "degree": "Network_degree"})
    net = net[pd.to_numeric(net["JGI_ID"], errors="coerce").notna()]
    net["JGI_ID"] = net["JGI_ID"].astype(int)

    return uniq, emap, fam, net


def build_table(uniq, emap, fam, net):
    df = emap.merge(fam, on="JGI_ID", how="left").merge(net, on="JGI_ID", how="left")
    df["Unique_to_MRL8996"] = df["JGI_ID"].isin(set(uniq["JGI_ID"]))
    df["Compartment"] = df["Chromosome"].apply(
        lambda c: "core" if c in CORE_CHR else "accessory")
    df["Structurally_modelled"] = df["Family_ID"].notna()
    df["In_previous_network"] = df["Network_community"].notna()
    df = df.merge(uniq[["JGI_ID", "Cluster_ID", "Cluster_size"]].drop_duplicates("JGI_ID"),
                  on="JGI_ID", how="left")
    return df


def enrichment(df):
    rows = []

    def add(name, table, note):
        odds, p = fisher_exact(table)
        rows.append({"Test": name, "a": table[0][0], "b": table[0][1],
                     "c": table[1][0], "d": table[1][1],
                     "odds_ratio": round(odds, 3), "p_value": p, "Note": note})

    u = df["Unique_to_MRL8996"]
    acc = df["Compartment"] == "accessory"
    add("unique_vs_accessory_chromosome",
        [[int((u & acc).sum()), int((u & ~acc).sum())],
         [int((~u & acc).sum()), int((~u & ~acc).sum())]],
        "rows: unique / shared; cols: accessory / core")

    am = df["AMAPEC"] == "Antimicrobial"
    add("unique_vs_antimicrobial",
        [[int((u & am).sum()), int((u & ~am).sum())],
         [int((~u & am).sum()), int((~u & ~am).sum())]],
        "rows: unique / shared; cols: antimicrobial / non-antimicrobial")

    mod = df["Structurally_modelled"]
    add("unique_vs_structurally_modelled",
        [[int((u & mod).sum()), int((u & ~mod).sum())],
         [int((~u & mod).sum()), int((~u & ~mod).sum())]],
        "rows: unique / shared; cols: has AF2 model in family set / not")

    return pd.DataFrame(rows)


def figure(df, out):
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    order = sorted(df["Chromosome"].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2),
                             gridspec_kw={"width_ratios": [1.55, 1]})

    # --- panel a: effectors per chromosome, unique fraction stacked
    ax = axes[0]
    uniq_n = [int(((df["Chromosome"] == c) & df["Unique_to_MRL8996"]).sum()) for c in order]
    shar_n = [int(((df["Chromosome"] == c) & ~df["Unique_to_MRL8996"]).sum()) for c in order]
    x = range(len(order))
    ax.bar(x, shar_n, color=OTHER, label="shared with Fol4287 and/or Fo47")
    ax.bar(x, uniq_n, bottom=shar_n, color=FOCAL, label="MRL8996-unique")
    ax.set_xticks(list(x))
    ax.set_xticklabels([c.replace("Chr", "") for c in order], fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Predicted effectors")
    ax.margins(y=0.18)
    n_acc_chr = [c for c in order if c not in CORE_CHR]
    if n_acc_chr:
        left = order.index(n_acc_chr[0]) - 0.5
        ax.axvspan(left, len(order) - 0.5, color="0.93", zorder=0)
        ax.text((left + len(order) - 0.5) / 2, ax.get_ylim()[1] * 0.34,
                "accessory\nchromosomes",
                ha="center", va="top", fontsize=7, color="0.35")
    ax.set_title("Unique effectors track total effector load per chromosome",
                 fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper center", ncol=1,
              bbox_to_anchor=(0.42, 1.02))
    ax.spines[["top", "right"]].set_visible(False)

    # --- panel b: effector density per Mb, core vs accessory
    ax = axes[1]
    lens = {}
    with open(HIC / "08_circos_features" / "chrom.sizes") as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and f[0].startswith("Chr"):
                lens[f[0]] = int(f[1])
    dens_all, dens_uni, labels, colors = [], [], [], []
    for comp in ("core", "accessory"):
        chrs = [c for c in order if (c in CORE_CHR) == (comp == "core")]
        mb = sum(lens.get(c, 0) for c in chrs) / 1e6
        sub = df[df["Chromosome"].isin(chrs)]
        dens_all.append(len(sub) / mb if mb else float("nan"))
        dens_uni.append(int(sub["Unique_to_MRL8996"].sum()) / mb if mb else float("nan"))
        labels.append("%s\n(%.1f Mb)" % (comp, mb))
        colors.append(OTHER if comp == "core" else "#8fa8bf")
    xx = [0, 1]
    ax.bar(xx, dens_all, color=OTHER, width=0.6, label="all effectors")
    ax.bar(xx, dens_uni, color=FOCAL, width=0.6, label="MRL8996-unique")
    for i, (a, b) in enumerate(zip(dens_all, dens_uni)):
        ax.text(i, a, " %.1f" % a, ha="center", va="bottom", fontsize=7, color="0.25")
        ax.text(i, b, "%.1f" % b, ha="center", va="bottom", fontsize=7, color="white")
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Effectors per Mb")
    ax.set_title("Effector density by compartment", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=0.16)

    fig.tight_layout()
    fig.savefig(str(out) + ".pdf", bbox_inches="tight")
    fig.savefig(str(out) + ".png", dpi=600, bbox_inches="tight")
    return fig


def main():
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    (REV / "figures").mkdir(parents=True, exist_ok=True)

    uniq, emap, fam, net = load()
    df = build_table(uniq, emap, fam, net)

    cols = ["JGI_ID", "Chromosome", "Start", "Compartment", "Unique_to_MRL8996",
            "Cluster_ID", "Cluster_size", "AMAPEC", "AM_probability", "pLDDT",
            "Family_ID", "Family_Rank", "Network_community", "Network_degree",
            "Structurally_modelled", "In_previous_network"]
    df[cols].sort_values(["Unique_to_MRL8996", "Chromosome", "Start"],
                         ascending=[False, True, True]).to_csv(
        REV / "tables" / "Table_MRL8996_unique_effectors.tsv", sep="\t", index=False)

    enr = enrichment(df)
    enr.to_csv(REV / "tables" / "unique_effector_enrichment.tsv", sep="\t", index=False)

    figure(df, REV / "figures" / "Fig_unique_effector_distribution")

    u = df[df["Unique_to_MRL8996"]]
    print("effectors with position data: %d (unique: %d)" % (len(df), len(u)))
    print("unique on accessory chromosomes: %d / %d"
          % (int((u["Compartment"] == "accessory").sum()), len(u)))
    print("unique predicted antimicrobial: %d / %d"
          % (int((u["AMAPEC"] == "Antimicrobial").sum()), len(u)))
    print("unique with AF2 model / structural family: %d / %d"
          % (int(u["Structurally_modelled"].sum()), len(u)))
    print("unique absent from previous (degree>=3) network: %d"
          % int((u["Structurally_modelled"] & ~u["In_previous_network"]).sum()))
    print(enr[["Test", "odds_ratio", "p_value"]].to_string(index=False))


if __name__ == "__main__":
    main()
