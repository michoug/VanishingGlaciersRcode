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

dat <- read_tsv("../Prokaryotes/Coverage/perc_total.txt")
meta <- read_tsv("../metadata_NOMIS.txt")
meta$Sample <- gsub("X","", meta$Sample)

meta$Site_c <- factor(meta$Site_c)

dat_meta <- dat %>%
  left_join(meta)%>%
  select(Site_c, Percentage)%>%
  filter(Site_c %in% c("Alps", "Norway"))%>%
  na.omit()%>%
  arrange(Site_c)

ggbetweenstats(dat_meta,x = Site_c,
               y = Percentage,
               type = "nonparametric",
               pairwise.comparisons = FALSE, 
               palette = "Paired",
               point.args = list(alpha = 1))

ggsave("readRecruitment.pdf",width = 7, height = 5)
