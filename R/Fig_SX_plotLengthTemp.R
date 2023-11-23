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

dat_ogt <- read_tsv("data/pMAGS_OptimalGrowthTemperature.txt")
dat_mags	<- read_tsv("data/pMAGs_checkm.txt")
magsRem <- read_tsv("data/pMAGsInRocks.txt")

dat_ogt$bins <- cut(dat_ogt$Temperature,breaks = c(0, 10, 20, 30, 40, 50, 60),include.lowest = TRUE)

dat_plot <- dat_mags %>%
  filter(!(MAGs %in% magsRem$MAGs)) %>%
  left_join(dat_ogt) %>%
  mutate(LengthNorm2 = LengthNorm2 / 1e6) %>%
  mutate(
    type = case_when(
      # Temperature > 45 ~ "Thermophiles",
      Temperature > 20 ~ "Mesophiles",
      .default = "Psychrophiles"
    )
  )%>%
  mutate(type = factor(type, levels = c("Psychrophiles", "Mesophiles", "Thermophiles")))


 
p_length <- ggbetweenstats(dat =dat_plot,x = type,
                           y = LengthNorm2, type = "nonparametric",
                           pairwise.comparisons = FALSE, 
                           pairwise.display = "none",
                           # palette = "Paired",
                           # point.args = list(size = 3),
                           centrality.label.args = list(size = 4, nudge_x = 0.4, nudge_y = 0.4),
                           # k = 1,
                           xlab = "Optimal Temperature",
                           ylab = "Estimated genome length (million bp)"
)+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.position = "none")

p2 <- AddLetters(dat_plot, plot = p_length, x= type, y=LengthNorm2, ytext = 17)
p2

ggsave_fitmax("Figures/Fig_SX_MAGs_temperature_length.pdf",p2, maxwidth = 10)

