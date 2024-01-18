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
  # mutate(
  #   cluster = case_when(
  #     cluster == "1" ~ "Diverse Taxa",
  #     cluster == "2" ~ "Myxcoccota &\n Bdellovibrionata",
  #     cluster == "3" ~ "Alpha and \n Gammaproteobacteria",
  #     cluster == "4" ~ "Bacteroidia",
  #     cluster == "5" ~ "Pactescibacteria",
  #     cluster == "6" ~ "Plantomycetes",
  #     .default = "Archae"
  #   )
  # ) %>%
  mutate(
    cluster = case_when(
      cluster == "1" ~ "A4",
      cluster == "2" ~ "A5",
      cluster == "3" ~ "A1",
      cluster == "4" ~ "A2",
      cluster == "5" ~ "A6",
      cluster == "6" ~ "A3",
      .default = "Archae"
    )
  ) %>%
  filter(!(cluster == "Archae"))

colors <- c(brewer.pal(6, "Set1"))
names(colors) <- c("A1", "A2", "A4", "A5", "A6", "A3")

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

ggsave_fitmax("Figures/Fig_SX_MAGs_cluster_length.pdf", p2, maxwidth = 10)
