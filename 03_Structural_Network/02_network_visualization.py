#!/usr/bin/env python3
"""
Stage 03.2 - Python rendering of the DALI structural-similarity network.
Reads artifacts written by 01_network_construction.R:
    OCE_struct_network.gml         (nodes: 'family'; edges: 'weight','etype')
    OCE_struct_network_nodes.tsv   (Dali_ID, family, degree, x, y, ...)
Writes:
    OCE_struct_network_py.png/.pdf (colored by family; reuses the R layout if present)
    OCE_struct_network.graphml     (Cytoscape-ready)
Same data as the R/ggraph figure - no re-clustering, just an alternative view.
"""
import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

GML, NODES, OUT = "OCE_struct_network.gml", "OCE_struct_network_nodes.tsv", "OCE_struct_network_py"

if not os.path.exists(GML):
    sys.exit(f"{GML} not found - run 01_network_construction.R first.")

G = nx.read_gml(GML, label="id")
fam = {n: int(float(d.get("family", 0))) for n, d in G.nodes(data=True)}

# reuse R layout coordinates for figure consistency, else spring layout
pos = {}
coords = None
if os.path.exists(NODES):
    nt = pd.read_csv(NODES, sep="\t")
    if {"x", "y"}.issubset(nt.columns):
        key = "name" if "name" in nt.columns else nt.columns[0]
        coords = {str(r[key]): (float(r["x"]), float(r["y"])) for _, r in nt.iterrows()}
if coords:
    rng = np.random.default_rng(42)
    for n, d in G.nodes(data=True):
        pos[n] = coords.get(str(d.get("name", n)), tuple(rng.random(2)))
else:
    pos = nx.spring_layout(G, weight="weight", seed=42)

fams = sorted(set(fam.values()))
cmap = plt.get_cmap("tab20", max(len(fams), 1))
cmap_by = {f: cmap(i) for i, f in enumerate(fams)}
node_colors = [cmap_by[fam[n]] for n in G.nodes()]

intra = [(u, v) for u, v, d in G.edges(data=True) if str(d.get("etype", "")) == "intra"]
inter = [(u, v) for u, v, d in G.edges(data=True) if str(d.get("etype", "")) == "inter"]

fig, ax = plt.subplots(figsize=(9, 9))
nx.draw_networkx_edges(G, pos, edgelist=inter, edge_color="0.85", width=0.3, alpha=0.4, ax=ax)
nx.draw_networkx_edges(G, pos, edgelist=intra, edge_color="0.55", width=0.5, alpha=0.6, ax=ax)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=28,
                       edgecolors="0.2", linewidths=0.3, ax=ax)
ax.axis("off"); ax.set_aspect("equal"); plt.tight_layout()
plt.savefig(f"{OUT}.png", dpi=400, bbox_inches="tight")
plt.savefig(f"{OUT}.pdf", bbox_inches="tight")
nx.write_graphml(G, "OCE_struct_network.graphml")
print(f"Wrote {OUT}.png / .pdf and OCE_struct_network.graphml")
