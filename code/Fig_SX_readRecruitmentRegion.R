## ---------------------------
##
## Script name: readRecruitmentRegion.R
##
## Purpose of script: Plot the coverage of the MAGs per region
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-14
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(statsExpressions)

source("customFunctions/plot_functions.R")

dat <- read_tsv("../Prokaryotes/Coverage/perc_total.txt")
meta <- read_tsv("../metadata_NOMIS.txt")
meta$Sample <- gsub("X","", meta$Sample)

# meta$Site_c <- factor(meta$Site_c)

dat_meta <- dat %>%
  left_join(meta)%>%
  select(Site_c, Percentage)%>%
  # filter(Site_c %in% c("Alps", "Caucasus"))%>%
  na.omit()%>%
  arrange(Site_c)%>%
  mutate(Site_c = if_else(Site_c == "New_Zealand", "New Zealand", Site_c))


p1 <- ggbetweenstats(dat =dat_meta,x = Site_c,
               y = Percentage, type = "nonparametric",
               pairwise.comparisons = FALSE, 
               pairwise.display = "none",
               palette = "Paired",
               point.args = list(alpha = 1, size = 3),
               centrality.label.args = list(size = 4, nudge_x = 0.4, nudge_y = 0.2),
               k = 1,
               xlab = element_blank(),
               ylab = "Percentage of Reads (%)"
               )+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.position = "none")

p1

ggsave_fitmax("Figures/Fig_SX_readRecruitment.pdf",p1, maxwidth = 10)
