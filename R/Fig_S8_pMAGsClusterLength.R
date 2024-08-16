## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-18
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


library(ggstatsplot)
library(tidyverse)
library(RColorBrewer)

source("customFunctions/CompactLetterDiplay_w_pairwiseComp_ggstatsplot.R")
source("customFunctions/plot_functions.R")

cluster <- read_tsv("data/pMAGS_clusters.tsv")
dat_mags	<- read_tsv("data/pMAGs_checkm.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")
magsRem <- read_tsv("data/pMAGsInRocks.txt")

cluster_sel <- cluster %>%
  filter(!(sub_grp_6 %in% magsRem$MAGs)) %>%
  select(sub_grp_6, d) %>%
  group_by(d) %>%
  mutate(nb = n())

colnames(cluster_sel) <- c("MAGs", "cluster", "nb")

dat_plot <- dat_mags %>%
  left_join(tax) %>%
  left_join(cluster_sel) %>%
  mutate(LengthNorm2 = LengthNorm2 / 1e6) %>%
  mutate(
    cluster = case_when(
      cluster == "1" ~ "4 _ Diverse Taxa",
      cluster == "2" ~ "5 _ Myxcoccota &\nBdellovibrionata",
      cluster == "3" ~ "1 _ Alpha and \nGammaproteobacteria",
      cluster == "4" ~ "2 _ Bacteroidia",
      cluster == "5" ~ "6 _ Pactescibacteria",
      cluster == "6" ~ "3 _ Plantomycetes",
      .default = "Archae"
    )
  ) %>%
  filter(!(cluster == "Archae"))

colors <- c(brewer.pal(6, "Set1"))
names(colors) <- c(
  "1 _ Alpha and \nGammaproteobacteria",
  "2 _ Bacteroidia",
  "4 _ Diverse Taxa",
  "5 _ Myxcoccota &\nBdellovibrionata",
  "6 _ Pactescibacteria",
  "3 _ Plantomycetes"
)

p_length <- ggbetweenstats(
  dat = dat_plot,
  x = cluster,
  y = LengthNorm2,
  type = "nonparametric",
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  # palette = "Set1",
  # point.args = list(size = 3),
  centrality.label.args = list(
    size = 4,
    nudge_x = 0.4,
    nudge_y = 0.4
  ),
  # k = 1,
  xlab = "Functional Cluster",
  ylab = "Estimated genome length (million bp)"
) +
  scale_color_manual(values = colors)+
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.position = "none"
  )

p2 <-
  AddLetters(
    dat_plot,
    plot = p_length,
    x = cluster,
    y = LengthNorm2,
    ytext = 17
  )
p2

ggsave_fitmax("Figures/Fig_S7_MAGs_cluster_length.pdf", p2, maxwidth = 10)
