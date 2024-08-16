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

dat_mags	<- read_tsv("data/pMAGs_checkm.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")
magsRem <- read_tsv("data/pMAGsInRocks.txt")
specgen <- read_tsv("data/pMAGs_spec_gen.tsv")


dat_plot <- dat_mags %>%
  filter(!(MAGs %in% magsRem$MAGs )) %>%
  left_join(tax) %>%
  right_join(specgen) %>%
  filter(sign != "NON SIGNIFICANT")%>%
  mutate(LengthNorm2 = LengthNorm2 / 1e6)



p_length <- ggbetweenstats(
  dat = dat_plot,
  x = sign,
  y = LengthNorm2,
  type = "nonparametric",
  pairwise.comparisons = FALSE,
  pairwise.display = "significant",
  # palette = "Set1",
  # point.args = list(size = 3),
  centrality.label.args = list(
    size = 4,
    nudge_x = 0.4,
    nudge_y = 0.4
  ),
  # k = 1,
  xlab = "Niche Breath",
  ylab = "Estimated genome length (million bp)"
) 

p_length

ggsave_fitmax("Figures/Fig_S6_length_MAGs_niche.pdf", p_length, maxwidth = 10)
