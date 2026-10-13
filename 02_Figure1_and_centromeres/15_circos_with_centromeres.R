# 15_circos_with_centromeres.R
# ============================================================================
# Redraws the MRL8996 feature circos (panel B of Figure 1) with the candidate
# centromeric regions from StainedGlass added as a dedicated track.
#
# This is a copy of MRL8996_HiC_data/11_Circos_Plot/MRL_circular_plot.r with
# three changes only:
#   1. all inputs are read from 11_Circos_Plot/ by absolute path, so nothing in
#      the original analysis directory is touched;
#   2. a new track (Track A2) draws the 16 centromere calls from
#      Revision_2026/results/centromeres_gc.bed (script 17: GC trough +
#      StainedGlass identity support; supersedes the script 14 calls);
#   3. the output goes to Revision_2026/figures/.
#
# The karyotype in chrom.sizes is ordered by descending length, which is NOT
# the PGA scaffold numbering; the centromere bed is written in the same
# chromosome namespace as chrom.sizes, so no renaming is applied to it.
#
# Usage: Rscript 15_circos_with_centromeres.R
# ============================================================================

library(circlize)

CIRC <- "/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_HiC_data/11_Circos_Plot"
REV  <- "/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026"
p <- function(...) file.path(CIRC, ...)

rename_scaffold <- function(x) {
  nums <- gsub("scaffold_", "", x)
  ifelse(grepl("^\\d+$", nums), sprintf("Chr%02d", as.numeric(nums)), x)
}

karyotype <- read.table(p("chrom.sizes"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(karyotype) <- c("chr", "end")
karyotype$end <- as.numeric(karyotype$end)
karyotype$start <- 0
karyotype <- karyotype[, c("chr", "start", "end")]
karyotype$chr <- rename_scaffold(karyotype$chr)
karyotype <- karyotype[karyotype$chr %in% sprintf("Chr%02d", 1:16), ]
karyotype <- na.omit(karyotype)

read_density <- function(file) {
  if (!file.exists(file)) return(NULL)
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "value")
  df$start <- as.numeric(df$start)
  df$end   <- as.numeric(df$end)
  df$value <- as.numeric(df$value)
  df$chr   <- rename_scaffold(df$chr)
  df <- df[df$chr %in% karyotype$chr, ]
  na.omit(df)
}

bgc_data <- read.table(p("circos_bgc_highlights.bed"), sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE)[, 1:3]
colnames(bgc_data) <- c("chr", "start", "end")
bgc_data$start <- as.numeric(bgc_data$start)
bgc_data$end   <- as.numeric(bgc_data$end)
bgc_data$chr   <- rename_scaffold(bgc_data$chr)
bgc_data       <- na.omit(bgc_data)

# --- candidate centromeres: GC trough + identity support (script 17) --------
cen <- read.table(file.path(REV, "results", "centromeres_gc.bed"),
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(cen) <- c("chr", "start", "end")
cen$start <- as.numeric(cen$start)
cen$end   <- as.numeric(cen$end)
cen <- cen[cen$chr %in% karyotype$chr, ]

gc_data   <- read_density(p("track_gc_content_20kb.txt"))
gene_data <- read_density(p("track_gene_density_20kb.txt"))
exon_data <- read_density(p("track_exon_density_20kb.txt"))
te_data   <- read_density(p("circos_track_TE_20kb.txt"))
sr_data   <- read_density(p("circos_track_SR_20kb.txt"))

chr_cols   <- rainbow(nrow(karyotype), s = 0.5, v = 0.9)
my_palette <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#f781bf")
CEN_COL    <- "#1a3a6b"  # navy, kept distinct from the black BGC marks

dir.create(file.path(REV, "figures"), showWarnings = FALSE, recursive = TRUE)
pdf(file.path(REV, "figures", "Figure_1B_circos_centromeres.pdf"), width = 12, height = 12)
circos.clear()

gaps <- c(rep(2, nrow(karyotype) - 1), 10)
circos.par(start.degree = 90, gap.degree = gaps,
           track.margin = c(0.002, 0.002), cell.padding = c(0, 0, 0, 0))

par(font = 2, cex = 1.1)
circos.genomicInitialize(karyotype, tickLabelsStartFromZero = TRUE, major.by = 1000000)

# Track A: karyotype + BGC clusters
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  idx <- which(karyotype$chr == chr)
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
              col = chr_cols[idx], border = "black", lwd = 1.5)
  chr_bgc <- bgc_data[bgc_data$chr == chr, ]
  if (nrow(chr_bgc) > 0) {
    for (j in 1:nrow(chr_bgc)) {
      circos.rect(chr_bgc$start[j], 0, chr_bgc$end[j], 1,
                  col = "black", border = "black", lwd = 2)
    }
  }
}, track.height = 0.04, bg.border = NA)

# Track A2 (NEW): candidate centromeric regions
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  r <- cen[cen$chr == chr, ]
  if (nrow(r) > 0) {
    for (j in 1:nrow(r)) {
      circos.rect(r$start[j], 0, r$end[j], 1,
                  col = CEN_COL, border = CEN_COL, lwd = 0.5)
    }
  }
}, track.height = 0.028, bg.border = "grey70", bg.col = "grey96")

# Track B: GC content
if (!is.null(gc_data) && nrow(gc_data) > 0) {
  circos.genomicTrack(gc_data, numeric.column = 4, ylim = c(0, 1), track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", col = my_palette[1], lwd = 0.6)
    }, bg.border = "black")
}

# Track C: gene density
if (!is.null(gene_data) && nrow(gene_data) > 0) {
  circos.genomicTrack(gene_data, numeric.column = 4, ylim = c(0, 1), track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", col = my_palette[2], lwd = 0.6)
    }, bg.border = "black")
}

# Track D: exon density
if (!is.null(exon_data) && nrow(exon_data) > 0) {
  circos.genomicTrack(exon_data, numeric.column = 4, ylim = c(0, 1), track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", col = my_palette[3], lwd = 0.6)
    }, bg.border = "black")
}

# Track E: transposable elements
if (!is.null(te_data) && nrow(te_data) > 0) {
  circos.genomicTrack(te_data, numeric.column = 4, ylim = c(0, 1), track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", col = my_palette[4], lwd = 0.6)
    }, bg.border = "black")
}

# Track F: simple repeats
if (!is.null(sr_data) && nrow(sr_data) > 0) {
  circos.genomicTrack(sr_data, numeric.column = 4, ylim = c(0, 0.2), track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", col = my_palette[5], lwd = 0.6)
    }, bg.border = "black")
}

legend("bottomright", bty = "n", cex = 0.85, pch = 15, pt.cex = 1.4,
       col = c("black", CEN_COL, my_palette[1:5]),
       legend = c("BGC clusters", "candidate centromere", "GC content",
                  "gene density", "exon density", "TE density", "simple repeats"))

dev.off()
circos.clear()
cat("wrote figures/Figure_1B_circos_centromeres.pdf\n")
