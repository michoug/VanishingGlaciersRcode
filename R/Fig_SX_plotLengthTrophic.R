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
library(ggpubr)
library(statsExpressions)
library(tidyverse)
library(RColorBrewer)
library("ggpubr")


source("customFunctions/CompactLetterDiplay_w_pairwiseComp_ggstatsplot.R")
source("customFunctions/plot_functions.R")

dat_troph <- read_tsv("../Prokaryotes/Microtrait/MAGsTrophicState.txt")
dat_mags	<- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGs_checkm_stats.txt")
magsRem <- read_tsv("../Prokaryotes/MAGsInRocks.txt")

dat_plot <- dat_mags %>%
  filter(!(MAGs %in% magsRem$MAGs))%>%
  left_join(dat_troph)%>%
  mutate(LengthNorm2 = LengthNorm2/1e6)%>%
  mutate(TrophicState = case_when(
    TrophicState == "Autotrophy" ~ "Chemolitoautotrophy",
    TrophicState == "Heterotrophy" ~ "Chemolitoheterotrophy\nMix",
    TrophicState == "Heterotrophy_aerobic_respiration" ~ "Chemoorganotrophy\naerobic respiration",
    TrophicState == "Heterotrophy_anaerobic_respiration" ~ "Chemoorganotrophy\nanaerobic respiration",
    TrophicState == "Heterotrophy_fermentation" ~ "Chemoorganotrophy\nfermentation",
    TrophicState == "Mixotrophy_Autotrophy-Phototrophy" ~ "Mixotrophy\nautotroph/phototroph",
    TrophicState == "Mixotrophy_Heterotrophy-Autotrophy" ~ "Mixotrophy\nheterotroph/autotroph",
    TrophicState == "Mixotrophy_Heterotrophy-Phototrophy" ~ "Mixotrophy\nheterotroph/phototroph",
    TrophicState == "Mixotrophy_Other"~ "Mixotrophy\nother",
    .default = TrophicState
  ))%>%
  filter(!(TrophicState == "Phototrophy"))

p_length <- ggbetweenstats(dat =dat_plot,x = TrophicState,
                     y = LengthNorm2, type = "nonparametric",
                     pairwise.comparisons = FALSE, 
                     pairwise.display = "none",
                     palette = "Paired",
                     # point.args = list(size = 3),
                     centrality.label.args = list(size = 4, nudge_x = 0.4, nudge_y = 0.4),
                     # k = 1,
                     xlab = "Functional Cluster",
                     ylab = "Estimated genome length (million bp)"
)+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.position = "none")

p2 <- AddLetters(dat_plot, plot = p_length, x= TrophicState, y=LengthNorm2, ytext = 17)
p2

ggsave_fitmax("Figures/Fig_SX_MAGs_trophic_length.pdf",p2, maxwidth = 10)

