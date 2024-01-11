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

dat <- read_tsv("data/pMAGs_reads_mapped.txt")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
meta$Sample <- gsub("X", "", meta$Sample)

# meta$Site_c <- factor(meta$Site_c)

dat_meta <- dat %>%
  left_join(meta) %>%
  select(Site_c, Percentage) %>%
  # filter(Site_c %in% c("Alps", "Caucasus"))%>%
  na.omit() %>%
  arrange(Site_c) %>%
  mutate(region = case_when(
    Site_c == "New_Zealand" ~ "Southern Alps,\n New Zealand",
    Site_c == "Nepal" ~ "Himalayas, Nepal",
    Site_c == "Caucasus" ~ "Caucasus Mountains,\n Russia",
    Site_c == "Kyrgyzstan" ~ "Pamir and Tian Shan,\n Kyrgyzstan",
    Site_c == "Alps" ~ "European Alps",
    Site_c == "Norway" ~ "Scandinavian Mountains,\n Norway",
    Site_c == "Greenland" ~ "Southwest Greenland",
    Site_c == "Uganda" ~ "Rwenzori Mountains,\n Uganda",
    Site_c == "Ecuador" ~ "Ecuadorian Andes"
  ))


p1 <- ggbetweenstats(
  dat = dat_meta,
  x = region,
  y = Percentage,
  type = "nonparametric",
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  palette = "Paired",
  point.args = list(alpha = 1, size = 3),
  centrality.label.args = list(
    size = 4,
    nudge_x = 0.4,
    nudge_y = 0.2
  ),
  k = 1,
  xlab = element_blank(),
  ylab = "Percentage of Reads (%)"
) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.position = "none"
  )

p1

ggsave_fitmax("Figures/Fig_S2_readRecruitment.pdf", p1, maxwidth = 18)
