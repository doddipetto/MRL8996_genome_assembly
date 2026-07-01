# ==============================================================================
# Effectorome Structural Heatmap (ComplexHeatmap)
#  - dendrogrammi (cladogrammi) attivi
#  - TUTTE le famiglie (cluster >= 5) riquadrate e numerate per dimensione
#  - riga AMAPEC in basso: nero = Antimicrobial, bianco = Non-antimicrobial,
#    grigio = nessuna predizione  (+ legenda dedicata)
#  - legenda TM-score: barra continua, gradiente liscio, senza tacche intermedie
# Metrica: TM-score (Max_TM_Score). ID normalizzati con canon() -> 482 proteine.
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(ComplexHeatmap); library(circlize)
})

# ----------------------------- CONFIG ----------------------------------------
ids_file    <- "FoxMRL8996_Final_Effectorome_IDs.txt"
tm_file     <- "tm_align_results_all_vs_all.tsv"
amapec_file <- "AMAPEC_Effectorome/results/prediction.csv"   
BOX_MIN     <- 5    # Dimensione minima per definire una "famiglia" (riquadro)
LABEL_MIN   <- 5    # Dimensione minima per stampare il numero sulla mappa
OUT_PDF     <- "ComplexHeatmap_Final_Effectorome_with_dendro.pdf"
OUT_TABLE   <- "Final_Effectorome_Families_Complex.tsv"
# -----------------------------------------------------------------------------

canon <- function(x) {
  x <- sub("\\.pdb$", "", x); x <- sub("_unrelaxed_rank.*$", "", x)
  x <- gsub("[^A-Za-z0-9.]+", "_", x); x <- gsub("_+", "_", x); gsub("^_|_$", "", x)
}
jgi_num <- function(x) str_match(x, "FoxMRL8996_(\\d+)")[, 2]   # ID numerico JGI

cat("1. IDs + matrice TM (normalizzati)...\n")
valid_ids <- readLines(ids_file) |> str_trim()
valid_ids <- canon(valid_ids[valid_ids != ""]) |> unique()
df_clean <- read_tsv(tm_file, show_col_types = FALSE) |>
  mutate(Protein_A = canon(Protein_A), Protein_B = canon(Protein_B)) |>
  filter(Protein_A %in% valid_ids & Protein_B %in% valid_ids)
active_ids <- sort(unique(c(df_clean$Protein_A, df_clean$Protein_B)))
n <- length(active_ids)
cat(sprintf("   proteine nella heatmap: %d\n", n))

cat("2. Matrice densa simmetrica...\n")
full_mat <- matrix(0, n, n, dimnames = list(active_ids, active_ids))
mA <- match(df_clean$Protein_A, active_ids); mB <- match(df_clean$Protein_B, active_ids)
ok <- !is.na(mA) & !is.na(mB)
full_mat[cbind(mA[ok], mB[ok])] <- df_clean$Max_TM_Score[ok]
full_mat[cbind(mB[ok], mA[ok])] <- df_clean$Max_TM_Score[ok]
diag(full_mat) <- 1.0

cat("3. Clustering gerarchico (average) + taglio a TM > 0.5...\n")
hc <- hclust(as.dist(1 - full_mat), method = "average")
clusters <- cutree(hc, h = 0.5)
fam_size <- sort(table(clusters), decreasing = TRUE)
fam_ids  <- as.integer(names(fam_size))[fam_size >= BOX_MIN]      # Solo famiglie >= 5 membri
fam_rank <- setNames(seq_along(fam_ids), fam_ids)                  # cid -> rango (1=più grande)
cat(sprintf("   famiglie (>= %d membri): %d  |  escluse/singleton (< %d membri): %d\n",
            BOX_MIN, length(fam_ids), BOX_MIN, sum(fam_size < BOX_MIN)))

cat("4. AMAPEC (antimicrobico)...\n")
if (!file.exists(amapec_file)) stop("prediction.csv non trovato: ", amapec_file)
amapec <- read_csv(amapec_file, show_col_types = FALSE)
pred <- amapec$Prediction[match(jgi_num(active_ids), as.character(amapec$`Protein ID`))]
prob <- amapec$`Probability of antimicrobial activity`[
          match(jgi_num(active_ids), as.character(amapec$`Protein ID`))]
plddt <- amapec$pLDDT[match(jgi_num(active_ids), as.character(amapec$`Protein ID`))]
pred[is.na(pred)] <- "Unknown"
cat(sprintf("   Antimicrobial: %d | Non-antimicrobial: %d | senza predizione: %d\n",
            sum(pred == "Antimicrobial"), sum(pred == "Non-antimicrobial"),
            sum(pred == "Unknown")))

# Tabella master esportata per tracciabilità biologica nel paper
tibble(Protein = active_ids, JGI_ID = jgi_num(active_ids),
       Family_ID = clusters, Family_Rank = fam_rank[as.character(clusters)],
       Antimicrobial = pred, AM_probability = prob, pLDDT = plddt) |>
  write_tsv(OUT_TABLE)

cat("5. Disegno...\n")
col_fun <- colorRamp2(seq(0, 1, length.out = 7),
                      c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))

# Riga annotazione AMAPEC in basso
levs <- intersect(c("Antimicrobial","Non-antimicrobial","Unknown"), unique(pred))
am_col <- c("Antimicrobial"="black","Non-antimicrobial"="white","Unknown"="grey85")[levs]
am_lab <- c("Antimicrobial"="Antimicrobial","Non-antimicrobial"="Non-antimicrobial",
            "Unknown"="No prediction")[levs]
bottom_ann <- HeatmapAnnotation(
  "AMAPEC" = pred,
  col = list(AMAPEC = am_col),
  gp = gpar(col = "grey70"),                 # Bordo celle per far risaltare il bianco
  simple_anno_size = unit(4.5, "mm"),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(AMAPEC = list(
    title = "Antimicrobial activity", at = levs, labels = am_lab,
    border = "grey50"))
)

pdf(OUT_PDF, width = 12, height = 11)
ht <- Heatmap(
  full_mat, name = "TM-score", col = col_fun,
  cluster_rows = hc, cluster_columns = hc,
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_dend_width = unit(20, "mm"), column_dend_height = unit(20, "mm"),
  show_row_names = FALSE, show_column_names = FALSE,
  border = TRUE, row_title = NULL, column_title = NULL,
  bottom_annotation = bottom_ann,
  heatmap_legend_param = list(
    title = "TM-score", title_position = "topcenter", direction = "vertical",
    legend_height = unit(4.5, "cm"), grid_width = unit(6, "mm"),
    at = c(0, 1), labels = c("0", "1"),      # Solo estremi -> gradiente continuo senza tacche intermedie
    border = "grey40")
)
ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
           merge_legend = TRUE, 
           padding = unit(c(2, 2, 2, 15), "mm")) # <-- Aumentato l'ultimo valore a 15 mm

# Riquadri + numeri per TUTTE le famiglie filtrate (contigue sul cladogramma)
ro <- row_order(ht); ord <- clusters[ro]
decorate_heatmap_body("TM-score", {
  for (cid in fam_ids) {
    idx <- which(ord == cid); if (!length(idx)) next
    st <- min(idx); en <- max(idx); sz <- en - st + 1
    xc <- (st - 1 + sz/2)/n; yc <- 1 - (st - 1 + sz/2)/n; wh <- sz/n
    grid.rect(x = xc, y = yc, width = wh, height = wh,
              gp = gpar(lwd = 1.2, col = "black", fill = "transparent"))
    if (fam_size[as.character(cid)] >= LABEL_MIN)
      grid.text(fam_rank[as.character(cid)],
                x = min(xc + wh/2 + 0.012, 0.985), y = yc,
                gp = gpar(fontsize = 7, fontface = "bold", col = "black"))
  }
})
invisible(dev.off())
cat(sprintf("\n--- FATTO: %s  (+ %s) ---\n", OUT_PDF, OUT_TABLE))
