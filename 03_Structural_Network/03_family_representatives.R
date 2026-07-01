#!/usr/bin/env Rscript
# =====================================================================
# family_representatives.R
# For the top-N families of the DALI structural network:
#   * pick a REPRESENTATIVE structure (medoid = highest summed within-
#     family Z-score = most structurally central member)
#   * report compactness QC (tight vs loose family) so you know whether a
#     single consensus annotation is even justified
#   * resolve each member to a structure file on disk
#   * emit a per-family table + an auto-generated PyMOL script (paths and
#     colours already filled in, colours matched to the network figure)
#
# Inputs : OCE_struct_network_nodes.tsv  (Dali_ID, family)   <- from main pipeline
#          ordered                       (DALI all-vs-all matrix)
#          dali_mapping_dict.tsv         (Dali_ID -> JGI/header)
# Run    : Rscript family_representatives.R
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble) })

cfg <- list(
  nodes_file   = "OCE_struct_network_nodes.tsv",
  matrix_file  = "ordered",
  mapping_file = "dali_mapping_dict.tsv",
  top_n        = 5,            # annotate the 5 biggest families
  loose_cutoff = 5.2,         # mean intra-Z below this => flag family as loose

  # ---- structure files: directory holding your cleaned .pdb models ----
  # The scheme (P001 / P001A / JGI id) is auto-detected from the filenames.
  struct_dir   = "PDB_clean",   # relative to where you run this, or use the absolute path
  out_table    = "family_representatives.tsv",
  out_pymol    = "render_representatives.pml"
)

## ---- read DALI 'ordered' matrix ----
read_dali_ordered <- function(path) {
  ln <- readLines(path); n <- as.integer(trimws(ln[1]))
  body <- ln[-1]; body <- body[nchar(trimws(body)) > 0]
  sp <- strsplit(body, "\t"); labs <- vapply(sp, `[`, character(1), 1L)
  M <- t(vapply(sp, function(x) as.numeric(x[-1]), numeric(n)))
  ids <- sub(".$", "", labs); dimnames(M) <- list(ids, ids); M
}
M <- read_dali_ordered(cfg$matrix_file)
M <- pmax(M, t(M)); diag(M) <- 0

nodes   <- read_tsv(cfg$nodes_file,   show_col_types = FALSE)
mapping <- read_tsv(cfg$mapping_file, show_col_types = FALSE)
stopifnot(all(c("Dali_ID","family") %in% names(nodes)))

## ---- structure-file index (auto-detect naming scheme from PDB_clean) ----
build_struct_index <- function(dir, mapping) {
  files <- list.files(dir, pattern = "\\.(pdb|cif|ent|mmcif)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files))
    stop("No .pdb/.cif structures found in '", dir, "'. Check cfg$struct_dir.")
  stem <- tools::file_path_sans_ext(basename(files))
  dali <- mapping$Dali_ID
  jgi  <- as.character(mapping$JGI_Protein_ID)
  # strip a trailing chain letter from stems like 'P001A' -> 'P001'
  stem_strip <- ifelse(grepl("^P[0-9]{3}[A-Za-z]$", stem),
                       substr(stem, 1, nchar(stem) - 1), stem)
  schemes <- list(
    dali        = match(dali, stem),        # files named P001
    `dali+chain`= match(dali, stem_strip),  # files named P001A
    jgi         = match(jgi,  stem)         # files named 10088
  )
  cover <- vapply(schemes, function(m) mean(!is.na(m)), numeric(1))
  best  <- names(which.max(cover))
  idx   <- schemes[[best]]
  paths <- files[idx]; names(paths) <- dali
  message(sprintf("Structure naming scheme: '%s'  (%.1f%% of proteins matched directly)",
                  best, 100 * max(cover)))
  # fuzzy fallback for any still-unresolved (substring contains JGI or Dali id)
  un <- which(is.na(paths))
  for (i in un) {
    cand <- files[grepl(jgi[i], basename(files), fixed = TRUE) |
                  grepl(dali[i], basename(files), fixed = TRUE)]
    if (length(cand)) paths[i] <- cand[1]
  }
  if (anyNA(paths))
    message(sprintf("[!] %d/%d proteins have no structure file (left as NA).",
                    sum(is.na(paths)), length(paths)))
  paths
}
struct_index <- build_struct_index(cfg$struct_dir, mapping)
resolve_struct <- function(dali_id, jgi_id = NULL) unname(struct_index[dali_id])

## ---- per-family representative + QC ----
fams <- sort(unique(nodes$family))
fams <- head(fams, cfg$top_n)
rep_rows <- list(); member_rows <- list()

for (f in fams) {
  members <- nodes$Dali_ID[nodes$family == f]
  members <- members[members %in% rownames(M)]
  S <- M[members, members, drop = FALSE]
  sumZ   <- rowSums(S)
  medoid <- members[which.max(sumZ)]
  off    <- S[upper.tri(S)]
  meanZ_all   <- mean(off)
  meanZ_edges <- if (any(off >= cfg$loose_cutoff)) mean(off[off >= cfg$loose_cutoff]) else NA_real_

  jgi_med <- mapping$JGI_Protein_ID[match(medoid, mapping$Dali_ID)]
  rep_rows[[as.character(f)]] <- tibble(
    family = f, n_members = length(members),
    representative_DaliID = medoid,
    representative_JGI = jgi_med,
    medoid_sumZ = round(max(sumZ), 1),
    mean_intraZ_all = round(meanZ_all, 2),
    mean_intraZ_edges = round(meanZ_edges, 2),
    compactness = ifelse(meanZ_all < cfg$loose_cutoff, "LOOSE (consider split)", "tight"),
    representative_file = resolve_struct(medoid,
                            mapping$JGI_Protein_ID[match(medoid, mapping$Dali_ID)])
  )
  member_rows[[as.character(f)]] <- tibble(
    family = f, Dali_ID = members,
    JGI_Protein_ID = mapping$JGI_Protein_ID[match(members, mapping$Dali_ID)],
    is_representative = members == medoid,
    sumZ_within_family = round(sumZ, 1),
    struct_file = vapply(seq_along(members), function(i)
      resolve_struct(members[i],
        mapping$JGI_Protein_ID[match(members[i], mapping$Dali_ID)]),
      character(1))
  ) %>% arrange(desc(sumZ_within_family))
}

rep_tbl    <- bind_rows(rep_rows)
member_tbl <- bind_rows(member_rows)

write_tsv(rep_tbl,    cfg$out_table)
write_tsv(member_tbl, sub("\\.tsv$", "_members.tsv", cfg$out_table))

cat("\n==== Top", cfg$top_n, "family representatives ====\n")
print(as.data.frame(rep_tbl), row.names = FALSE)
n_missing <- sum(is.na(member_tbl$struct_file))
if (n_missing) cat(sprintf(
  "\n[!] %d/%d member structures NOT found under '%s'. Check cfg$struct_dir.\n",
  n_missing, nrow(member_tbl), cfg$struct_dir))

## ---- auto-generate PyMOL render script (colours matched to network figure) ----
# same palette as the main pipeline so PyMOL families match the network plot
nf_total <- length(unique(nodes$family))
pal <- grDevices::palette.colors(min(nf_total, 36), "Polychrome 36")
if (nf_total > 36) pal <- grDevices::colorRampPalette(pal)(nf_total)
rgb_of <- function(hex) { v <- col2rgb(hex)/255; sprintf("[%.3f,%.3f,%.3f]", v[1], v[2], v[3]) }

pml <- c(
  "# Auto-generated by family_representatives.R",
  "# Renders each family's representative (opaque, spectrum) with members",
  "# superposed (grey, thin) to show fold conservation. cealign = structure-",
  "# based, works at low sequence identity.",
  "bg_color white", "set ray_opaque_background, 0", "set cartoon_transparency, 0",
  "set ray_shadows, 0", "set antialias, 2", ""
)
for (i in seq_len(nrow(rep_tbl))) {
  f   <- rep_tbl$family[i]
  med <- rep_tbl$representative_DaliID[i]
  medf<- rep_tbl$representative_file[i]
  if (is.na(medf)) { pml <- c(pml, sprintf("# family %d: representative file missing, skipped", f)); next }
  mem <- member_tbl %>% filter(family == f, !is.na(struct_file), Dali_ID != med)
  col <- sprintf("fam%d", f)
  pml <- c(pml,
    sprintf("# ---------------- family %d (rep %s) ----------------", f, med),
    "delete all",
    sprintf("set_color %s, %s", col, rgb_of(pal[f])),
    sprintf("load %s, rep%d", medf, f))
  for (j in seq_len(nrow(mem)))
    pml <- c(pml, sprintf("load %s, m_%d_%d", mem$struct_file[j], f, j))
  for (j in seq_len(nrow(mem)))
    pml <- c(pml, sprintf("cealign rep%d, m_%d_%d", f, f, j))
  pml <- c(pml,
    "hide everything",
    sprintf("show cartoon, rep%d", f),
    "show cartoon, m_*",
    sprintf("spectrum count, rainbow, rep%d", f),
    "color grey70, m_*",
    "set cartoon_transparency, 0.6, m_*",
    sprintf("orient rep%d", f),
    "ray 1400, 1400",
    sprintf("png family%d_representative.png, dpi=300", f),
    "")
}
writeLines(pml, cfg$out_pymol)
cat(sprintf("\nPyMOL script written: %s   (run: pymol -cq %s)\n",
            cfg$out_pymol, cfg$out_pymol))
