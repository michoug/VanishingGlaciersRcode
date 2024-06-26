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
library(ggbreak)
library(RColorBrewer)
library(ggpubr)

source("customFunctions/plot_functions.R")


dat_troph <- read_tsv("data/pMAGs_trophicState.txt")
cluster <- read_tsv("data/pMAGS_clusters.tsv")
magsRem <- read_tsv("data/pMAGsInRocks.txt")


dat_summary <- dat_troph %>%
  filter(!(MAGs %in% magsRem$MAGs)) %>%
  inner_join(cluster, by = c("MAGs" = "sub_grp_6")) %>%
  group_by(d, TrophicState) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(d) %>%
  mutate(Frequence = round(n / sum(n), 3)) %>%
  ungroup() %>%
  complete(d, TrophicState, fill = list(n = 0, Frequence = 0)) %>%
  mutate(
    d = case_when(
      d == "1" ~ "4 - Diverse Taxa",
      d == "2" ~ "5 - Myxcoccota &\nBdellovibrionata",
      d == "3" ~ "1 - Alpha and \nGammaproteobacteria",
      d == "4" ~ "2 - Bacteroidia",
      d == "5" ~ "6 - Pactescibacteria",
      d == "6" ~ "3 - Plantomycetes",
      .default = "Archae"
    )
  ) %>%
  filter(!(d == "Archae"))

dat_all <- dat_troph %>%
  filter(!(MAGs %in% magsRem$MAGs)) %>%
  inner_join(cluster, by = c("MAGs" = "sub_grp_6")) %>%
  group_by(TrophicState) %>%
  summarise(n = n()) %>%
  mutate(Frequence = round(n / sum(n), 3))

dat_all$d <- "All"
#write_tsv(dat_summary, "MAGsTrophicSpegGenSummary.txt")

dat_final <- rbind(dat_summary, dat_all)

dat_final <- dat_final %>%
  mutate(
    TrophicState = case_when(
      TrophicState == "Autotrophy" ~ "Chemolitoautotrophy",
      TrophicState == "Heterotrophy" ~ "Chemolitoheterotrophy\nMix",
      TrophicState == "Heterotrophy_aerobic_respiration" ~ "Chemoorganotrophy\naerobic respiration",
      TrophicState == "Heterotrophy_anaerobic_respiration" ~ "Chemoorganotrophy\nanaerobic respiration",
      TrophicState == "Heterotrophy_fermentation" ~ "Chemoorganotrophy\nfermentation",
      TrophicState == "Mixotrophy_Autotrophy-Phototrophy" ~ "Mixotrophy\nautotroph/phototroph",
      TrophicState == "Mixotrophy_Heterotrophy-Autotrophy" ~ "Mixotrophy\nheterotroph/autotroph",
      TrophicState == "Mixotrophy_Heterotrophy-Phototrophy" ~ "Mixotrophy\nheterotroph/phototroph",
      TrophicState == "Mixotrophy_Other" ~ "Mixotrophy\nother",
      .default = TrophicState
    )
  ) %>%
  filter(!(TrophicState == "Phototrophy"))

colors <- c(brewer.pal(6, "Set1"), "black")
names(colors) <- c(
  "4 - Diverse Taxa",
  "5 - Myxcoccota &\nBdellovibrionata",
  "1 - Alpha and \nGammaproteobacteria",
  "2 - Bacteroidia",
  "6 - Pactescibacteria",
  "3 - Plantomycetes",
  "All"
)

p1 <-
  ggplot(dat_final,
         aes(x = TrophicState, y = Frequence * 100, fill = as.factor(d))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(fill = "Functional Cluster", x = "Trophic State", y = "Percentage (%)") +
  scale_y_break(c(50, 85)) +
  theme_classic2() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.2),
    panel.grid.major.y = element_line(color = "black", linewidth = 0.1),
    text = element_text(colour = "black", size = 12)
  )

ggsave_fitmax("Figures/Fig_3c_TrophicStateCluster.pdf", p1, maxwidth = 10)
