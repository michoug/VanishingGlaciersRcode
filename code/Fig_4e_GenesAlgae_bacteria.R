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
library(ggpubr)
library(vegan)

source("customFunctions/plot_functions.R")

dat_select <- read_tsv("../Prokaryotes/EggNog/allMAGs_egg_KO.txt.gz")
cluster <- read_tsv("../Prokaryotes/Cluster/agnes_groups_6.tsv")
magsRem <- read_tsv("../Prokaryotes/MAGsInRocks.txt")
genes_path <- read_tsv("../Prokaryotes/Algae_Bacteria_Genes/listFunctionsAlgaeBacteria_pathway.txt")
genes_mod <- read_tsv("../Prokaryotes/Algae_Bacteria_Genes/listFunctionsAlgaeBacteria_modules.txt")

genes_path$Description <- NULL 
genes_mod$Category <- gsub("(Biotin biosynthesis),.*","\\1", genes_mod$Category, perl = T)
genes_path$Category <- gsub("(Biofilm formation).*","\\1", genes_path$Category, perl = T)

colnames(genes_path) <- colnames(genes_mod)
all_genes <- rbind(genes_path, genes_mod)

all_genes <- all_genes %>%
  select(-Module)%>%
  distinct()%>%
  group_by(Category)%>%
  mutate(number = n())

dat_filter <- dat_select%>%
  right_join(all_genes, by =c("KEGG_ko"="KO"), multiple = "all",relationship = "many-to-many")%>%
  filter(!(MAGs %in% magsRem$MAGs))%>%
  group_by(MAGs, KEGG_ko, Category, number)%>%
  summarise(n = n())%>%
  ungroup()%>%
  mutate(n = if_else(n > 1, 1,n))%>%
  group_by(MAGs, Category, number)%>%
  reframe(n = n()/number)%>%
  ungroup()%>%
  distinct()%>%
  complete(MAGs, Category, fill= list(n =0))%>%
  right_join(cluster, by = c("MAGs" = "sub_grp_6"))%>%
  na.omit()

dat_plot <- dat_filter %>%
  group_by(d, Category)%>%
  summarise(mean = mean(n))%>%
  ungroup()%>%
  complete(d, Category, fill= list(mean =0))%>%
  mutate(d = case_when(
    d == "1" ~ "Diverse Taxa",
    d == "2" ~ "Myxcoccota &\nBdellovibrionata",
    d == "3" ~ "Alpha and \nGammaproteobacteria",
    d == "4" ~ "Bacteroidia",
    d == "5" ~ "Patescibacteria",
    d == "6" ~ "Plantomycetes",
    .default = "Archae"
  ))%>%
  filter(!(d == "Archae"))


dat_plot_t <- dat_plot %>%
  pivot_wider(names_from = Category, values_from = mean)

dat_plot_t <- as.data.frame(dat_plot_t)
rownames(dat_plot_t) <- dat_plot_t$d
dat_plot_t$d <- NULL

data.dist <- vegdist(t(dat_plot_t) , method = "euclidian",na.rm = T)
col.clus <- hclust(data.dist, "aver")
ord <- col.clus$order
a <- colnames(dat_plot_t)[ord]

dat_plot$Category <- factor(dat_plot$Category, levels = a)

p1 <- ggplot(dat_plot, aes(x = as.character(d), y = Category, fill = mean*100))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "darkgreen")+
  labs(x = NULL, y = NULL, fill = "Proprotion (%)")+
  theme_pubr()
p1
ggsave_fitmax("Figures/Fig_4e_Heatmap_cat_algae_bact.pdf", p1, maxwidth = 12)
