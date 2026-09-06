#!/usr/bin/env Rscript
# 21_effector_map_unique.R
# =========================
# The author's linear effector map (MRL8996_HiC_data/plot_map_effector_position.R)
# with the MRL8996-unique effectors marked.  Everything else -- geometry, sizes,
# colours, theme, axis treatment, output dimensions -- is copied verbatim from
# the original script, so panel a keeps its published appearance; the only
# addition is a black triangle above the bar at each unique effector and the
# corresponding legend key.
#
# "Unique" is the CD-HIT definition used throughout the revision: the effector
# falls in a cluster (>= 50% identity) that contains no Fo47 and no Fol4287
# protein (script 01/03).  135 of the 513 effectors are unique; 122 of them are
# among the 482 with a structural model that the map plots.
#
# Inputs : MRL8996_HiC_data/effector_map_data.tsv
#          MRL8996_HiC_data/MRL8996_chromosome_level.fasta.fai
#          Revision_2026/tables/Table_MRL8996_unique_effectors.tsv
# Output : Revision_2026/figures/MRL8996_Linear_Effector_Map_unique.pdf
#
# Usage: Rscript 21_effector_map_unique.R

library(ggplot2)
library(dplyr)
library(stringr)

HIC <- "/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_HiC_data"
REV <- "/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026"

data <- read.delim(file.path(HIC, "effector_map_data.tsv"), check.names = FALSE)
fai <- read.delim(file.path(HIC, "MRL8996_chromosome_level.fasta.fai"), header = FALSE)
colnames(fai) <- c("Chromosome", "Length", "Offset", "LineBases", "LineWidth")

uniq <- read.delim(file.path(REV, "tables", "Table_MRL8996_unique_effectors.tsv"))
uniq_ids <- as.character(uniq$JGI_ID[uniq$Unique_to_MRL8996 == "True"])
data$Unique <- as.character(data[["Protein ID"]]) %in% uniq_ids
message(sprintf("unique effectors on the map: %d of %d plotted (%d unique in total)",
                sum(data$Unique), nrow(data), length(uniq_ids)))

chr_lengths <- fai %>%
  filter(grepl("Chr", Chromosome)) %>%
  select(Chromosome, Length) %>%
  mutate(chr_num = as.numeric(str_extract(Chromosome, "\\d+"))) %>%
  arrange(chr_num)

data <- data %>%
  mutate(chr_num = as.numeric(str_extract(Chromosome, "\\d+"))) %>%
  arrange(chr_num)

chr_lengths$Chromosome <- factor(chr_lengths$Chromosome,
                                 levels = unique(chr_lengths$Chromosome))
data$Chromosome <- factor(data$Chromosome, levels = levels(chr_lengths$Chromosome))
uniq_pts <- data %>% filter(Unique)

p <- ggplot() +
  geom_segment(data = chr_lengths,
               aes(x = 0, xend = Length / 1e6,
                   y = as.numeric(Chromosome),
                   yend = as.numeric(Chromosome)),
               color = "#E0E0E0", size = 6, lineend = "butt") +

  geom_segment(data = data,
               aes(x = Start / 1e6, xend = Start / 1e6,
                   y = as.numeric(Chromosome) - 0.25,
                   yend = as.numeric(Chromosome) + 0.25,
                   color = Prediction),
               size = 0.8) +

  # added layer: MRL8996-unique effectors, marked above the bar
  geom_point(data = uniq_pts,
             aes(x = Start / 1e6, y = as.numeric(Chromosome) - 0.42,
                 shape = "MRL8996-unique"),
             color = "black", fill = "black", size = 1.5, stroke = 0) +

  scale_color_manual(values = c("Antimicrobial" = "red",
                                "Non-antimicrobial" = "blue")) +
  scale_shape_manual(values = c("MRL8996-unique" = 25)) +

  scale_y_reverse(breaks = 1:nrow(chr_lengths),
                  labels = levels(chr_lengths$Chromosome)) +

  theme_classic() +
  labs(x = "Genomic Position (Mb)", y = NULL,
       color = "AMAPEC Activity", shape = NULL) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2,
                              override.aes = list(size = 2.2))) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

out <- file.path(REV, "figures", "MRL8996_Linear_Effector_Map_unique.pdf")
ggsave(out, plot = p, width = 10, height = 7, dpi = 300)
message("wrote ", out)
