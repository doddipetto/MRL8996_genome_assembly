#!/usr/bin/env python3
"""
07_structural_network_with_singletons.py
========================================
Rebuilds the DALI structural-similarity network of the MRL8996 effectorome
WITHOUT the degree >= 3 filter used in the published figure, so that
low-connectivity nodes and structural singletons are visible, and colours every
node by whether its effector is MRL8996-unique (absent from the Fol4287 and
Fo47 effectoromes) or shared.

Node fill carries the Louvain family (one colour per family) and the node
outline carries provenance (red ring = MRL8996-unique), so that unique
effectors sitting inside a family remain identifiable as members of it.

Method follows structural_network/dali_structure_network.R: the DaliLite
all-vs-all 'ordered' matrix is read, edges are kept at Z >= 5.2, and families
are Louvain communities on the weighted graph. The only deliberate departures
are (i) no degree filter - all 482 modelled effectors are retained - and (ii)
Louvain from networkx rather than igraph, so that the script has no R
dependency. Communities with fewer than MIN_FAMILY members are treated as
unassigned and drawn in the outer ring together with the degree-0 singletons.

Inputs (read-only):
    .../Dali_Effectorome_Run/ordered                 482 x 482 Z-score matrix
    .../Dali_Effectorome_Run/dali_mapping_dict.tsv   Dali_ID -> JGI ID
    Revision_2026/tables/Table_MRL8996_unique_effectors.tsv   (script 03)

Outputs:
    Revision_2026/tables/network_singletons_nodes.tsv
    Revision_2026/tables/network_family_composition.tsv
    Revision_2026/figures/Fig_structural_network_singletons.pdf / .png (600 dpi)

Usage: python3 07_structural_network_with_singletons.py
"""

import argparse
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
DALI = (ADODDI / "MRL8996_genome_paper" / "JGI_reference" / "FoxMRL8996" / "Mycocosm"
        / "Annotation" / "Filtered_Models___best__" / "Proteins" / "Dali_Effectorome_Run")

Z_THRESHOLD = 5.2          # as in the published network
MIN_FAMILY = 3             # communities smaller than this go to the outer ring
FAM_PAD = 0.55             # clearance between family discs, in mean-family-radius units
SEED = 42
PANEL = False              # set by --panel; drops the title for the Figure 4 panel

# --- colour scheme -------------------------------------------------------
# Node FILL encodes structural family: one colour per Louvain family, a single
# grey for everything outside the families (structural singletons and the
# sub-MIN_FAMILY communities). Node OUTLINE encodes provenance: MRL8996-unique
# effectors get a red ring. Provenance is drawn as the outline rather than the
# fill because 67 of the 122 unique effectors sit INSIDE families - recolouring
# them would erase the family membership the figure exists to show.
RING_C = "#d62728"        # MRL8996-unique outline
NOFAM_C = "#c3c8cd"       # singletons and small communities
EDGE_C = "#9aa3ab"


def family_palette(n):
    """n visually distinct fills, red-free so the unique ring stays readable."""
    src = [mpl.colormaps["tab20"], mpl.colormaps["tab20b"],
           mpl.colormaps["tab20c"]]
    out = []
    for cm in src:
        for i in range(cm.N):
            r, g, b = cm(i)[:3]
            if r == max(r, g, b) and r - g > 0.13 and r - b > 0.13:  # drop red hues
                continue
            if min(r, g, b) > 0.80:                      # drop near-whites
                continue
            if max(r, g, b) - min(r, g, b) < 0.16:       # drop greys: reserved
                continue                                 # for the non-family fill
            c = (round(r, 3), round(g, 3), round(b, 3))
            if c not in out:
                out.append(c)
    if len(out) < n:
        raise SystemExit("palette too small: %d < %d" % (len(out), n))
    # spread the picks across the concatenated palette so neighbours differ
    idx = np.linspace(0, len(out) - 1, n).round().astype(int)
    return [out[i] for i in idx]


def read_ordered(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    n = int(lines[0].strip())
    labels, rows = [], []
    for line in lines[1:]:
        f = line.split("\t")
        labels.append(f[0][:-1])           # strip trailing chain letter: P086A -> P086
        rows.append([float(x) for x in f[1:n + 1]])
    M = np.array(rows)
    assert M.shape == (n, n), M.shape
    M = np.maximum(M, M.T)
    np.fill_diagonal(M, 0.0)
    return labels, M


def sunflower_disc(sizes, spread=1.30):
    """Family centres on a phyllotaxis disc, radius ~ sqrt(size)."""
    n = len(sizes)
    radii = np.sqrt(np.asarray(sizes, float) / math.pi)
    rmax = math.sqrt(sum(sizes) / math.pi)
    ga = math.pi * (3 - math.sqrt(5))
    i = np.arange(1, n + 1)
    r = rmax * np.sqrt(i / n) * spread
    return np.column_stack([r * np.cos(i * ga), r * np.sin(i * ga)]), radii


def main():
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    (REV / "figures").mkdir(parents=True, exist_ok=True)

    labels, M = read_ordered(DALI / "ordered")
    mapping = pd.read_csv(DALI / "dali_mapping_dict.tsv", sep="\t")
    dali2jgi = dict(zip(mapping["Dali_ID"], mapping["JGI_Protein_ID"]))

    ann = pd.read_csv(REV / "tables" / "Table_MRL8996_unique_effectors.tsv", sep="\t")
    is_unique = dict(zip(ann["JGI_ID"], ann["Unique_to_MRL8996"]))

    G = nx.Graph()
    G.add_nodes_from(labels)
    iu = np.triu_indices_from(M, k=1)
    keep = M[iu] >= Z_THRESHOLD
    for a, b, w in zip(np.array(iu[0])[keep], np.array(iu[1])[keep], M[iu][keep]):
        G.add_edge(labels[a], labels[b], weight=float(w))
    print("nodes %d  edges (Z>=%.1f) %d" % (G.number_of_nodes(), Z_THRESHOLD,
                                            G.number_of_edges()))

    comms = nx.community.louvain_communities(G, weight="weight", seed=SEED)
    comms = sorted((c for c in comms), key=len, reverse=True)
    fam_of = {}
    fam_id = 0
    for c in comms:
        members = [n for n in c]
        if len(members) >= MIN_FAMILY and any(G.degree(n) > 0 for n in members):
            fam_id += 1
            for n in members:
                fam_of[n] = fam_id
    unassigned = [n for n in G.nodes if n not in fam_of]
    print("families (>=%d members) %d ; outer-ring nodes %d"
          % (MIN_FAMILY, fam_id, len(unassigned)))

    # ---------------- layout ----------------
    rng = np.random.default_rng(SEED)
    pos = {}
    fams = list(range(1, fam_id + 1))
    sizes = [sum(1 for n in fam_of if fam_of[n] == f) for f in fams]
    centres, radii = sunflower_disc(sizes)

    # Each family is laid out on its own first, so its true extent is known
    # before the families are placed relative to one another.
    members_of, local, extent = {}, {}, np.zeros(len(fams))
    for k, (f, rad) in enumerate(zip(fams, radii)):
        members = [n for n in G.nodes if fam_of.get(n) == f]
        sub = G.subgraph(members)
        loc = nx.spring_layout(sub, weight="weight", seed=SEED, k=None)
        arr = np.array([loc[n] for n in members])
        arr -= arr.mean(axis=0)
        norm = np.abs(arr).max() or 1.0
        arr = arr / norm * rad * 0.82
        members_of[f], local[f] = members, arr
        extent[k] = np.hypot(arr[:, 0], arr[:, 1]).max() if len(arr) else rad * 0.1

    # The sunflower disc spaces the centres by sqrt(size) only, which lets a
    # small family land on the rim of a large one (families 2 and 7 in the
    # first version). Relax the centres until no two enclosing discs overlap:
    # deterministic, and it keeps the sunflower as the starting arrangement.
    keep_out = extent + FAM_PAD * extent.mean()
    for _ in range(600):
        moved = 0.0
        for i in range(len(fams)):
            for j in range(i + 1, len(fams)):
                dv = centres[j] - centres[i]
                dist = float(np.hypot(*dv)) or 1e-9
                gap = keep_out[i] + keep_out[j] - dist
                if gap > 0:
                    step = 0.5 * gap * dv / dist
                    centres[i] -= step
                    centres[j] += step
                    moved = max(moved, gap)
        centres -= centres.mean(axis=0)
        if moved < 1e-4:
            break
    print("family-centre relaxation: residual overlap %.4g" % moved)

    for f, c in zip(fams, centres):
        for n, xy in zip(members_of[f], local[f]):
            pos[n] = (c[0] + xy[0], c[1] + xy[1])

    disc_r = max(np.hypot(*np.array(list(pos.values())).T)) if pos else 1.0
    ring_r = disc_r * 1.22
    # outer ring: singletons (degree 0) on the outermost circle, small
    # communities just inside it
    singles = [n for n in unassigned if G.degree(n) == 0]
    smallcomm = [n for n in unassigned if G.degree(n) > 0]
    for group, rr in ((smallcomm, ring_r), (singles, ring_r * 1.12)):
        for i, n in enumerate(sorted(group)):
            th = 2 * math.pi * i / max(len(group), 1) + rng.uniform(-0.004, 0.004)
            pos[n] = (rr * math.cos(th), rr * math.sin(th))

    # ---------------- node table ----------------
    rows = []
    for n in G.nodes:
        jgi = dali2jgi.get(n)
        rows.append({
            "Dali_ID": n,
            "JGI_ID": jgi,
            "Unique_to_MRL8996": is_unique.get(jgi),
            "family": fam_of.get(n, pd.NA),
            "degree": G.degree(n),
            "is_singleton": G.degree(n) == 0,
            "x": pos[n][0], "y": pos[n][1],
        })
    nodes = pd.DataFrame(rows)
    nodes.to_csv(REV / "tables" / "network_singletons_nodes.tsv", sep="\t", index=False)

    comp = (nodes.dropna(subset=["family"])
            .groupby("family")
            .agg(n_members=("Dali_ID", "size"),
                 n_unique=("Unique_to_MRL8996", "sum"),
                 median_degree=("degree", "median"))
            .reset_index())
    comp["pct_unique"] = (100 * comp["n_unique"] / comp["n_members"]).round(1)
    comp.to_csv(REV / "tables" / "network_family_composition.tsv", sep="\t", index=False)

    # are unique effectors over-represented among singletons / unassigned?
    u = nodes["Unique_to_MRL8996"].fillna(False).astype(bool)
    ring = nodes["family"].isna()
    odds, p = fisher_exact([[int((u & ring).sum()), int((u & ~ring).sum())],
                            [int((~u & ring).sum()), int((~u & ~ring).sum())]])
    print("unique in outer ring %d/%d ; shared %d/%d ; OR=%.2f p=%.4g"
          % (int((u & ring).sum()), int(u.sum()),
             int((~u & ring).sum()), int((~u).sum()), odds, p))

    # ---------------- figure ----------------
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(7.2, 7.2))

    for a, b, d in G.edges(data=True):
        same = fam_of.get(a) is not None and fam_of.get(a) == fam_of.get(b)
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
                color=EDGE_C, linewidth=0.28 if same else 0.16,
                alpha=0.55 if same else 0.25, zorder=1)

    fam_cols = dict(zip(fams, family_palette(fam_id)))
    fill = [fam_cols[f] if f in fam_cols else NOFAM_C
            for f in nodes["family"].map(lambda v: int(v) if pd.notna(v) else None)]
    nodes = nodes.assign(_fill=fill)

    # shared effectors first, then the unique ones on top with the red ring
    shared = nodes[~nodes["Unique_to_MRL8996"].fillna(False).astype(bool)]
    uniq = nodes[nodes["Unique_to_MRL8996"].fillna(False).astype(bool)]
    ax.scatter(shared["x"], shared["y"], s=17, c=list(shared["_fill"]),
               edgecolors="0.30", linewidths=0.25, zorder=2)
    ax.scatter(uniq["x"], uniq["y"], s=27, c=list(uniq["_fill"]),
               edgecolors=RING_C, linewidths=0.85, zorder=3)

    for f, (cx, cy) in zip(fams, centres):
        ax.text(cx, cy, str(f), ha="center", va="center", fontsize=6.5,
                fontweight="bold", zorder=4,
                bbox=dict(boxstyle="circle,pad=0.18", facecolor="white",
                          edgecolor="0.4", linewidth=0.4, alpha=0.9))

    ax.text(0, ring_r * 1.12 + disc_r * 0.08,
            "outermost ring: %d structural singletons (no DALI Z \u2265 5.2 partner)\n"
            "inner ring: %d nodes in communities of < %d members    \u00b7    "
            "numbered balls: %d Louvain families, one fill colour each"
            % (len(singles), len(smallcomm), MIN_FAMILY, fam_id),
            ha="center", va="bottom", fontsize=7, color="0.35")

    n_u_ring = int((u & ring).sum())
    if not PANEL:
        ax.set_title("Structural-similarity network of the MRL8996 effectorome, "
                     "singletons retained\n"
                     "%d of %d MRL8996-unique effectors fall outside the "
                     "Louvain families" % (n_u_ring, int(u.sum())),
                     fontsize=10, pad=30)
    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=6.5,
               markerfacecolor=NOFAM_C, markeredgecolor="0.30",
               label="outside the families (singleton or community < %d)" % MIN_FAMILY),
        Line2D([], [], marker="o", linestyle="", markersize=7.5,
               markerfacecolor="white", markeredgecolor=RING_C, markeredgewidth=1.4,
               label="MRL8996-unique (red outline, any fill)"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower center",
              bbox_to_anchor=(0.5, -0.055 if PANEL else -0.085),
              ncol=2 if PANEL else 1, handletextpad=0.6, columnspacing=1.8)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    stem = ("Fig_structural_network_singletons_panel" if PANEL
            else "Fig_structural_network_singletons")
    fig.savefig(REV / "figures" / (stem + ".pdf"), bbox_inches="tight")
    fig.savefig(REV / "figures" / (stem + ".png"), dpi=600, bbox_inches="tight")
    print("families: %d ; singletons: %d ; small communities: %d"
          % (fam_id, len(singles), len(smallcomm)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", action="store_true",
                    help="panel version for Figure 4: no title, two-column "
                         "legend, written as *_panel.pdf/.png")
    PANEL = ap.parse_args().panel
    main()
