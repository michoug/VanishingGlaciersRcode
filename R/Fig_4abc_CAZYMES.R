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

library(tidyverse)
library(ggstatsplot)
library(ggpubr)


source("customFunctions/plot_functions.R")
source("customFunctions/CompactLetterDiplay_w_pairwiseComp_ggstatsplot.R")

dat_select <-
  read_tsv("../Prokaryotes/CAZymes_Sulfatase/CAZYME_egg.txt.gz")
dat_cont_genes <- read_tsv("../Prokaryotes/mags_cont_genes.txt.gz")
cluster <- read_tsv("../Prokaryotes/Cluster/agnes_groups_6.tsv")
element <-
  read_tsv("../Prokaryotes/CAZymes_Sulfatase/CAZyme_Element_Baumgen.txt")
magsRem <- read_tsv("../Prokaryotes/MAGsInRocks.txt")

element_sep <- element %>%
  separate_rows(CAZyme)

fucoidan <- c("GH29", "GH95", "GH141", "GH107")

dat_cazy_mags <- dat_select %>%
  left_join(dat_cont_genes, by = c("#query" = "Genes")) %>%
  separate_rows(CAZy) %>%
  right_join(
    element_sep,
    by = c("CAZy" = "CAZyme"),
    multiple = "all",
    relationship = "many-to-many"
  ) %>%
  filter(!(MAGs %in% magsRem$MAGs)) %>%
  group_by(Element, MAGs) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  complete(Element, MAGs, fill = list(n = 0)) %>%
  left_join(cluster, by = c("MAGs" = "sub_grp_6")) %>%
  na.omit() %>%
  group_by(Element) %>%
  mutate(sum = sum(n), mean = mean(n)) %>%
  ungroup() %>%
  filter(sum > 0) %>%
  filter(mean > 0.2) %>%
  filter(!(Element == "Mannan")) %>%
  mutate(
    d = case_when(
      d == "1" ~ "Diverse Taxa",
      d == "2" ~ "Myxcoccota &\nBdellovibrionata",
      d == "3" ~ "Alpha and \nGammaproteobacteria",
      d == "4" ~ "Bacteroidia",
      d == "5" ~ "Patescibacteria",
      d == "6" ~ "Plantomycetes",
      .default = "Archae"
    )
  ) %>%
  filter(!(d == "Archae")) %>%
  filter(n > 0)


p1 <- grouped_ggbetweenstats(
  data = dat_cazy_mags,
  x = d,
  y = n,
  type = "n",
  grouping.var = Element,
  results.subtitle = T,
  ylab = "Number of Genes per MAGs",
  xlab = "Functional Clusters",
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  ggtheme = theme_classic2(),
  palette = "Set1"
)

p1


ggsave_fitmax("Figures/Fig_4abc_CAZymes_cluster.pdf", p1, maxwidth = 15)


dat_mostAbund <- dat_cazy_mags %>%
  group_by(c, Element) %>%
  summarise(mean = mean(n))
