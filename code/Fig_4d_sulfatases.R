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
library(statsExpressions)
library(ggstatsplot)
library(ggpubr)

source("customFunctions/plot_functions.R")
source("customFunctions/CompactLetterDiplay_w_pairwiseComp_ggstatsplot.R")

dat <- read_tsv("../Prokaryotes/CAZymes_Sulfatase/NOMIS_MAGs_sulfatase_hmm_sulfAtlas.tsv")
dat_blast <- read_tsv("../Prokaryotes/CAZymes_Sulfatase/NOMIS_MAGs_sulfatase_60_sulfAtlas.tsv")
dat_cont_genes <- read_tsv("../Prokaryotes/mags_cont_genes.txt.gz")
cluster <- read_tsv("../Prokaryotes/Cluster/agnes_groups_6.tsv")
blast <- read_tsv("../Prokaryotes/CAZymes_Sulfatase/NOMIS_MAGs_sulfatase.txt.gz")
descrition <- read_tsv("../Prokaryotes/CAZymes_Sulfatase/sulfatase_description.txt")
magsRem <- read_tsv("../Prokaryotes/MAGsInRocks.txt")

# descrition_fucan <- descrition %>%
#   filter(Activity == "Fucan")

sulfGenes <- dat$protein
sulfGenes <- as.data.frame(sulfGenes)
sulfGenesblast <- dat_blast$protein
sulfGenesblast <- as.data.frame(sulfGenesblast)

colnames(sulfGenesblast) <- colnames(sulfGenes)

sulfGenes_all <- distinct(rbind(sulfGenes, sulfGenesblast))

sulfGenes_all$MAGs <- gsub("(.*)_\\d+_\\d+$", "\\1", sulfGenes_all$sulfGenes, perl = T)

dat_sulf <- sulfGenes_all %>%
  filter(!(MAGs %in% magsRem$MAGs))%>%
  left_join(blast, by = c("sulfGenes" = "qseqid"))%>%
  separate(sseqid, into = c("bad","id","na"), sep = "\\|")%>%
  mutate(cat = gsub("(.*?)_(.*)","\\2", perl = T, id))%>%
  separate_rows(cat, sep = ";")%>%
  select(MAGs, cat)%>%
  left_join(descrition, by = c("cat" = "Subfamily"), multiple = "all",relationship = "many-to-many")%>%
  na.omit()%>%
  filter(!(Activity == "Other"))%>%
  group_by(MAGs,Activity)%>%
  summarise(n = n())%>%
  ungroup()%>%
  complete(MAGs, Activity, fill = list(n = 0))%>%
  left_join(cluster, by = c("MAGs" = "sub_grp_6"))%>%
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

p1 <-  ggbetweenstats(
  data = dat_sulf,
  x = d,
  y = n,
  type = "n",
  pairwise.comparisons = FALSE, 
  pairwise.display = "none",
  palette = "Set1",
  results.subtitle = T)
# p1

p2 <- AddLetters(dat_sulf, plot = p1, x= d, y=n,ytext = 130, type = "n")+
  theme(text = element_text(colour = "black", size = 12))+
  labs(y = "Number of Genes per MAGs", x = NULL, color = "Functional Clusters")+
  theme_classic2()
p2

ggsave_fitmax("Figures/Fig_4d_sulfatases.pdf", p2, maxwidth = 12)
