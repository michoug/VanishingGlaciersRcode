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

source("customFunctions/plot_functions.R")

dat <- read_tsv("../Prokaryotes/LightGenes/lucifer_heatmap_select.tsv")
cluster <- read_tsv("../Prokaryotes/Cluster/agnes_groups_6.tsv")
covmax <- read_tsv("../Prokaryotes/MAGs_cov_sum.txt")
magsRem <- read_tsv("../Prokaryotes/MAGsInRocks.txt")

cluster <- cluster %>%
  group_by(d)%>%
  mutate(n = n())

dat_merge <- dat %>%
  pivot_longer(cols = !c(GeneName, Description))%>%
  filter(!(name %in% magsRem$MAGs))%>%
  left_join(cluster, by = c("name" = "sub_grp_6"))%>%
  left_join(covmax, by = c("name" = "MAGs"))


dat_temp <- dat_merge %>%
  filter(value > 0)%>%
  mutate(new_value = if_else(value > 1, 1, 1))%>%
  filter(GeneName == "Photo_RC")%>%
  filter(c == "p__Proteobacteria")%>%
  select(name)%>%
  distinct()%>%
  left_join(covmax, by = c("name" = "MAGs"))%>%
  summarise(sum(percentage))

dat_filter <- dat_merge%>%
  filter(value > 0)%>%
  mutate(new_value = if_else(value > 1, 1, 1))%>%
  group_by(Description,d)%>%
  reframe(sum = sum(new_value), Proportion = sum(new_value)*100 / n, Coverage = sum(percentage)*100)%>%
  distinct()%>%
  na.omit()%>%
  complete(Description, d, fill = list(sum =0, Proportion = 0, Coverage = 0))%>%
  mutate(Description = case_when(
    Description == "Chlorophyll A-B binding protein" ~ "Chlorophyll A-B\nbinding protein",
    Description == "Heliorhodopsin_Pfam" ~ "Heliorhodopsin",
    Description == "Photosynthetic reaction centre protein" ~ "Photosynthetic\nreaction centre",
    Description == "Photosystem-II-phosphoprotein" ~ "Photosystem II\nphosphoprotein",
    Description == "Photosystem I psaA/psaB protein" ~ "Photosystem I\npsaA/psaBprotein",
    .default = Description
  ))%>%
  filter(!(Description == "Anabaena sensory rhodopsin transducer"))%>%
  mutate(d = case_when(
    d == "1" ~ "Diverse Taxa",
    d == "2" ~ "Myxcoccota &\nBdellovibrionata",
    d == "3" ~ "Alpha and \nGammaproteobacteria",
    d == "4" ~ "Bacteroidia",
    d == "5" ~ "Pactescibacteria",
    d == "6" ~ "Plantomycetes",
    .default = "Archae"
  ))%>%
  filter(!(d == "Archae"))

(p1 <- ggplot(dat_filter, aes(x = Description, y = Proportion, fill = as.factor(d)))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_brewer(palette = "Set1")+
  labs(fill = "Functional Cluster", x = NULL, y = "Percentage (%)")+
  scale_y_break(c(20,35))+
  theme_classic()+
  theme(axis.text.x = element_text(angle =45, hjust = 1),
        # panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.2),
        panel.grid.major.y = element_line(color = "black",linewidth = 0.1),
        text = element_text(colour = "black", size = 12)))




ggsave_fitmax("Figures/Fig_X_lucifer_cluster.pdf", p1, maxwidth = 10)

