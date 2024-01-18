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
library(patchwork)

source("customFunctions/plot_functions.R")
source("customFunctions/CompactLetterDiplay_w_pairwiseComp_ggstatsplot.R")

dat_select <-
  read_tsv("data/pMAGs_CAZYME.txt.gz")
dat_cont_genes <- read_tsv("data/pMAGS_contigs_genes.txt.gz")
cluster <- read_tsv("data/pMAGS_clusters.tsv")
element <- read_tsv("data/Other/CAZyme_Element_Baumgen.txt")
magsRem <- read_tsv("data/pMAGsInRocks.txt")

dat <- read_tsv("data/pMAGs_sulfatase_sulfatlas_hmm.tsv")
dat_blast <- read_tsv("data/pMAGs_sulfatase_sulfatlas_blast_60.tsv")
blast <- read_tsv("data/pMAGs_sulfatase_blast.txt.gz")
description <- read_tsv("data/Other/sulfatase_description.txt")

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
  complete(MAGs, Element, fill = list(n = 0)) %>%
  left_join(cluster, by = c("MAGs" = "sub_grp_6")) %>%
  na.omit() %>%
  group_by(Element) %>%
  mutate(sum = sum(n), mean = mean(n)) %>%
  ungroup() %>%
  filter(sum > 0) %>%
  filter(mean > 0.2) %>%
  filter(!(Element %in% c("Mannan", "Xylan"))) %>%
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
  select(-c(sum, mean)) %>%
  filter(n > 0)

sulfGenes <- dat$protein
sulfGenes <- as.data.frame(sulfGenes)
sulfGenesblast <- dat_blast$protein
sulfGenesblast <- as.data.frame(sulfGenesblast)

colnames(sulfGenesblast) <- colnames(sulfGenes)

sulfGenes_all <- distinct(rbind(sulfGenes, sulfGenesblast))

sulfGenes_all$MAGs <-
  gsub("(.*)_\\d+_\\d+$", "\\1", sulfGenes_all$sulfGenes, perl = T)

dat_sulf <- sulfGenes_all %>%
  filter(!(MAGs %in% magsRem$MAGs)) %>%
  left_join(blast, by = c("sulfGenes" = "qseqid")) %>%
  separate(sseqid, into = c("bad", "id", "na"), sep = "\\|") %>%
  mutate(cat = gsub("(.*?)_(.*)", "\\2", perl = T, id)) %>%
  separate_rows(cat, sep = ";") %>%
  select(MAGs, cat) %>%
  left_join(
    description,
    by = c("cat" = "Subfamily"),
    multiple = "all",
    relationship = "many-to-many"
  ) %>%
  na.omit() %>%
  filter(!(Activity == "Other")) %>%
  group_by(MAGs, Activity) %>%
  summarise(numb = n()) %>%
  ungroup() %>%
  complete(MAGs, Activity, fill = list(numb = 0)) %>%
  group_by(MAGs) %>%
  summarise(n = sum(numb)) %>%
  mutate(Element = "Sulfatase") %>%
  left_join(cluster, by = c("MAGs" = "sub_grp_6")) %>%
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

dat_final <- rbind(dat_cazy_mags, dat_sulf)

plot_list <- list()

for (i in unique(dat_final$Element)) {
  print(i)
  
  dat_temp <- dat_final %>%
    filter(Element == !!i)
  
  p1 <- ggbetweenstats(
    data = dat_temp,
    x = d,
    y = n,
    type = "n",
    results.subtitle = T,
    # ylab = "Number of Genes per MAGs", xlab = "Functional Clusters",
    pairwise.comparisons = FALSE,
    pairwise.display = "none",
    ggtheme = theme_classic2(),
    palette = "Set1",
    title = i
  ) +
    theme(text = element_text(colour = "black", size = 12),
          legend.position = "none") +
    labs(y = "Number of Genes per MAGs", x = NULL, color = NULL)
  
  # textVal <- max(dat_temp$n) + 5
  # print(textVal)
  #
  # p2 <- AddLetters(dat_temp, plot = p1, x= d, y=n,ytext = max(dat_temp$n) + 5, type = "n")+
  #   theme_classic2()+
  #   theme(text = element_text(colour = "black", size = 12),
  #         legend.position="none")+
  #   labs(y = "Number of Genes per MAGs", x = NULL, color = NULL)
  
  plot_list[[i]] <- p1
  
}

p <- wrap_plots(plot_list, 2, 2)

ggsave_fitmax("Figures/Fig_4abcd_CAZymes_Sulfatases_cluster.pdf", p,
              maxwidth = 15)
