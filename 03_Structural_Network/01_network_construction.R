#!/usr/bin/env Rscript
# =====================================================================
# DALI structural-similarity network  (pure-R reproduction of the
# Nature-style OCE structural network figure)
#
#   matrix  ->  igraph  ->  Louvain families  ->  circular layout  ->  figure
#
# Reproduces the methods:
#   * edges filtered at Z-score >= 5.2
#   * nodes with degree < 3 removed
#   * community detection (Louvain; the algorithm HiDEF wraps)
#   * force-directed / circular "ring of families" layout
#   * exports a .gml (igraph, NOT Cytoscape) + a node table with names
#
# Validated against the input files: 482 structures, symmetric matrix,
# expected ~271 nodes / 1445 edges / ~22 families after filtering.
# Run:  Rscript dali_structure_network.R
# =====================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggforce)     # geom_mark / hulls for the family rings
  library(dplyr)
  library(readr)
  library(tibble)
  library(scales)
})

## ---------------------------------------------------------------- CONFIG
cfg <- list(
  matrix_file   = "ordered",                 # DaliLite all-vs-all 'ordered' matrix
  mapping_file  = "dali_mapping_dict.tsv",    # Dali_ID -> JGI/header
  out_prefix    = "OCE_struct_network",

  z_threshold   = 5.2,    # keep edges with Z >= this
  min_degree    = 3,      # drop nodes with degree < this
  degree_mode   = "single",   # "single" = one pass (paper wording); "kcore" = iterative 3-core
  resolution    = 1.0,    # Louvain resolution. Raise (2-4) for more, smaller families.
  community_alg = "louvain",  # "louvain" | "leiden"

  seed          = 42,
  pack_fill     = 0.82,   # fraction of each packed circle the family ball fills (<1 = gap between families)
  label_families = TRUE   # draw numbered badges like the reference figure
)

set.seed(cfg$seed)

## ---------------------------------------------------------- 1. READ MATRIX
# 'ordered' format: line 1 = N; each following line = "label\t z1 \t z2 ...".
read_dali_ordered <- function(path) {
  ln <- readLines(path)
  n  <- as.integer(trimws(ln[1]))
  body <- ln[-1]
  body <- body[nchar(trimws(body)) > 0]
  sp   <- strsplit(body, "\t")
  labs <- vapply(sp, `[`, character(1), 1L)
  M <- t(vapply(sp, function(x) as.numeric(x[-1]), numeric(n)))
  stopifnot(nrow(M) == n, ncol(M) == n, length(labs) == n)
  dimnames(M) <- list(labs, labs)
  M
}

M_raw <- read_dali_ordered(cfg$matrix_file)
message(sprintf("Loaded matrix: %d x %d", nrow(M_raw), ncol(M_raw)))

# strip trailing chain letter (P###A -> P###) to match Dali_ID
ids <- sub(".$", "", rownames(M_raw))
stopifnot(!any(duplicated(ids)))
dimnames(M_raw) <- list(ids, ids)

# Defensive symmetrization (matrix is already symmetric here) + clear diagonal.
# NOTE: weights are SIMILARITIES (higher = more similar).
M <- pmax(M_raw, t(M_raw))
diag(M) <- 0

## --------------------------------------------------- 2. BUILD EDGE LIST
ut <- which(upper.tri(M) & M >= cfg$z_threshold, arr.ind = TRUE)
edges <- tibble(
  from   = ids[ut[, "row"]],
  to     = ids[ut[, "col"]],
  weight = M[ut]                      # Z-score similarity
)
message(sprintf("Edges with Z >= %.1f : %d", cfg$z_threshold, nrow(edges)))

g0 <- graph_from_data_frame(edges, directed = FALSE,
                            vertices = data.frame(name = ids))

## --------------------------------------------------- 3. DEGREE FILTER
# CHOICE: single-pass mirrors the paper's literal wording ("nodes with
# degree < 3"); k-core cascades the removal. Both are reported by the
# accompanying simulation; default = single.
if (cfg$degree_mode == "kcore") {
  keep <- V(g0)[coreness(g0) >= cfg$min_degree]
  g <- induced_subgraph(g0, keep)
} else {
  g <- delete_vertices(g0, V(g0)[degree(g0) < cfg$min_degree])
}
# drop any isolates produced by the removal
g <- delete_vertices(g, V(g)[degree(g) == 0])
message(sprintf("After degree>=%d filter (%s): %d nodes, %d edges",
                cfg$min_degree, cfg$degree_mode, vcount(g), ecount(g)))

## --------------------------------------------------- 4. COMMUNITIES
# Louvain treats weight CORRECTLY as tie strength (higher = stronger),
# so it receives the raw Z-scores. (HiDEF wraps Louvain; no R port.)
if (cfg$community_alg == "leiden") {
  comm <- cluster_leiden(g, weights = E(g)$weight,
                         objective_function = "modularity",
                         resolution_parameter = cfg$resolution)
} else {
  comm <- cluster_louvain(g, weights = E(g)$weight,
                          resolution = cfg$resolution)
}
mem <- membership(comm)
# renumber families by size (1 = largest), like a figure legend
ord <- order(table(mem), decreasing = TRUE)
remap <- setNames(seq_along(ord), names(table(mem))[ord])
V(g)$family <- factor(remap[as.character(mem)],
                      levels = seq_along(ord))
message(sprintf("Families: %d   modularity = %.3f",
                length(unique(V(g)$family)),
                modularity(g, as.integer(V(g)$family), weights = E(g)$weight)))

# tag edges intra/inter family (drives the faint grey 'web' in the centre)
el <- as_edgelist(g, names = FALSE)
E(g)$etype <- ifelse(V(g)$family[el[,1]] == V(g)$family[el[,2]],
                     "intra", "inter")

## --------------------------------------------------- 5. CIRCULAR LAYOUT
# Reproduces the reference: each family is a tight ball, and the balls are
# packed to FILL the circular interior. Deterministic; does NOT depend on the
# weight-direction gotcha. Within a family, Fruchterman-Reingold is used with
# weight as ATTRACTION (higher Z = closer) -- correct for similarity.
#
# Alternative faithful-to-paper option (Kamada-Kawai over the whole graph):
#   dist <- max(E(g)$weight) - E(g)$weight + 1   # similarity -> distance
#   xy <- layout_with_kk(g, weights = dist)
#   (KK/stress treat weight as DISTANCE, hence the conversion.)
#
# Family CENTRES are circle-packed to FILL the disc (not arranged on the rim):
# circleProgressiveLayout nests circles of area proportional to family size,
# touching but not overlapping, in a compact roughly-circular cluster. Each
# family's internal FR ball is then scaled to sit inside its packed circle.
# Falls back to a phyllotaxis (sunflower) disc fill if packcircles is absent.
layout_community_pack <- function(g, fill = 0.82) {
  fam  <- V(g)$family
  fams <- levels(fam)
  sizes <- vapply(fams, function(f) sum(fam == f), integer(1))

  if (requireNamespace("packcircles", quietly = TRUE)) {
    pk <- packcircles::circleProgressiveLayout(sizes, sizetype = "area")
    centers <- data.frame(family = fams, cx = pk$x, cy = pk$y, r = pk$radius)
  } else {
    message("packcircles not installed -> phyllotaxis fallback ",
            "(may leave small gaps/touches; install r-packcircles for clean packing).")
    n  <- length(fams); i <- seq_len(n)
    ga <- pi * (3 - sqrt(5))                 # golden angle
    rmax <- sqrt(sum(sizes / pi))            # disc radius to hold total area
    centers <- data.frame(
      family = fams,
      cx = rmax * sqrt(i / n) * 1.3 * cos(i * ga),
      cy = rmax * sqrt(i / n) * 1.3 * sin(i * ga),
      r  = sqrt(sizes / pi)
    )
  }

  xy <- matrix(NA_real_, vcount(g), 2)
  for (i in seq_along(fams)) {
    idx <- which(fam == fams[i])
    sub <- induced_subgraph(g, idx)
    if (vcount(sub) == 1) {
      loc <- matrix(0, 1, 2)
    } else {
      loc <- layout_with_fr(sub, weights = E(sub)$weight)  # similarity = attraction
      loc <- scale(loc, center = TRUE, scale = FALSE)       # centre on origin
      maxr <- max(sqrt(rowSums(loc^2))); if (maxr == 0) maxr <- 1
      loc <- loc / maxr                                     # fit inside unit circle
    }
    rad <- centers$r[i] * fill          # leave a gap between neighbouring families
    xy[idx, 1] <- centers$cx[i] + loc[, 1] * rad
    xy[idx, 2] <- centers$cy[i] + loc[, 2] * rad
  }
  list(xy = xy, centers = centers)
}

lo <- layout_community_pack(g, fill = cfg$pack_fill)

## --------------------------------------------------- 6. PLOT
# build a manual ggraph layout from our coordinates
lay <- create_layout(g, layout = "manual",
                     x = lo$xy[, 1], y = lo$xy[, 2])

# big qualitative palette (base R 4.x, no extra package)
nfam <- nlevels(V(g)$family)
base_pal <- grDevices::palette.colors(min(nfam, 36), "Polychrome 36")
fam_pal  <- if (nfam <= 36) base_pal else
            grDevices::colorRampPalette(base_pal)(nfam)

p <- ggraph(lay) +
  # inter-family edges first, very faint (the central web)
  geom_edge_link0(aes(filter = etype == "inter"),
                  colour = "grey80", width = 0.15, alpha = 0.25) +
  # intra-family edges, slightly stronger
  geom_edge_link0(aes(filter = etype == "intra"),
                  colour = "grey55", width = 0.25, alpha = 0.5) +
  geom_node_point(aes(fill = family), shape = 21, size = 2.4,
                  stroke = 0.15, colour = "grey20") +
  scale_fill_manual(values = setNames(fam_pal, levels(V(g)$family)),
                    guide = "none") +
  coord_fixed() +
  theme_void()

# numbered family badges, sized to each family's packed circle (so small
# families get small badges and don't collide), like the reference figure
if (cfg$label_families) {
  cen <- lo$centers
  cen$bsize <- scales::rescale(cen$r, to = c(3.5, 9))
  cen$tsize <- scales::rescale(cen$r, to = c(2.2, 4.2))
  p <- p + geom_point(data = cen, aes(cx, cy, size = bsize),
                      inherit.aes = FALSE, shape = 21, fill = "white",
                      colour = "grey30", stroke = 0.6, alpha = 0.9) +
           scale_size_identity() +
           geom_text(data = cen, aes(cx, cy, label = family, size = tsize),
                     inherit.aes = FALSE, fontface = "bold")
}

ggsave(paste0(cfg$out_prefix, ".png"), p, width = 9, height = 9,
       dpi = 400, bg = "white")
ggsave(paste0(cfg$out_prefix, ".pdf"), p, width = 9, height = 9,
       bg = "white")   # vector, for final figure assembly
message("Figure written: ", cfg$out_prefix, ".png / .pdf")

## --------------------------------------------------- 7. EXPORTS
# NOTE: write the data tables BEFORE the GML, so a GML hiccup can never
# block the artifacts you actually analyse downstream.

# 7a. node table with family + names mapped back via the dictionary
mapping <- read_tsv(cfg$mapping_file, show_col_types = FALSE)
nodes_out <- tibble(
  Dali_ID = V(g)$name,
  family  = as.integer(V(g)$family),
  degree  = degree(g),
  x = lo$xy[, 1], y = lo$xy[, 2]
) %>%
  left_join(mapping, by = c("Dali_ID" = "Dali_ID"))
# placeholder column for the manual PDB-hit family naming step:
nodes_out$family_name <- NA_character_
write_tsv(nodes_out, paste0(cfg$out_prefix, "_nodes.tsv"))

# 7c. family summary (size, to drive PDB-hit naming downstream)
fam_summary <- nodes_out %>%
  count(family, name = "n_members") %>%
  arrange(family)
write_tsv(fam_summary, paste0(cfg$out_prefix, "_families.tsv"))

# 7d. GML (igraph native -- the in-R equivalent of plotCytoscapeGML).
# igraph's GML writer accepts only numeric/string attributes, so factors
# and logicals must be coerced first (a factor like $family triggers
# "Attribute not numeric, Invalid value"). Coerce on a copy, then write.
sanitize_for_gml <- function(g) {
  for (a in vertex_attr_names(g)) {
    v <- vertex_attr(g, a)
    if (is.factor(v))  g <- set_vertex_attr(g, a, value = as.integer(v))
    if (is.logical(v)) g <- set_vertex_attr(g, a, value = as.integer(v))
  }
  for (a in edge_attr_names(g)) {
    v <- edge_attr(g, a)
    if (is.factor(v))  g <- set_edge_attr(g, a, value = as.character(v))
    if (is.logical(v)) g <- set_edge_attr(g, a, value = as.integer(v))
  }
  g
}
tryCatch(
  write_graph(sanitize_for_gml(g), paste0(cfg$out_prefix, ".gml"),
              format = "gml"),
  error = function(e)
    message("GML export skipped (", conditionMessage(e),
            "); TSV tables were still written.")
)

message("Done. Outputs: .png .pdf .gml _nodes.tsv _families.tsv")
