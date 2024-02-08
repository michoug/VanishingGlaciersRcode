## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-15
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(ape)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(RColorBrewer)
library(tidyverse)

source("customFunctions/plot_functions.R")

ArTree	<- read.tree("data/pMAGs_arch_gtdtk.tree")
dat <- read_tsv("data/pMAGs_tax.tsv")
dataset <- read_tsv("data/pMAGs_Presence_Datasets.txt")
covmax <- read_tsv("data/pMAGs_cov_sum.txt")

archDat <- dat %>%
  filter(d == "d__Archaea") %>%
  mutate(p = gsub("p__", "", p)) %>%
  mutate(Abundance = 0.2)

archDatset <- dataset %>%
  filter(MAGs %in% archDat$MAGs)

archDatset <- as.data.frame(archDatset)
rownames(archDatset) <- archDatset$MAGs
archDatset$MAGs <- NULL

archDatset$Busi <-
  gsub("Busi", "Busi et al., Nat Com, 2022", archDatset$Busi)
archDatset$ENSEMBLE <-
  gsub("ENSEMBLE", "Michoud et al, L&O, 2023", archDatset$ENSEMBLE)
archDatset$Tibet <-
  gsub("Tibet", "Tibetan Glacier Genome and Gene", archDatset$Tibet)
archDatset$Tara <- gsub("Tara", "Tara Oceans", archDatset$Tara)

archcov <- covmax %>%
  filter(MAGs %in% archDat$MAGs) %>%
  select(MAGs, count)

archcov <- as.data.frame(archcov)
rownames(archcov) <- archcov$MAGs
archcov$MAGs <- NULL

archcov$count <- log10(archcov$count)

a	<- split(archDat$MAGs, archDat$p)
tree	<- groupOTU(ArTree, a)

colourCount	<- 31
getPaletteArc =	colorRampPalette(brewer.pal(9, "Set1"))
archColor	<- getPaletteArc(length(unique(archDat$p)) + 1)
archColor[1]	<- "black"


p1_arch	<-
  ggtree(tree,
         layout = 'circular',
         aes(color = group)) +
  #branch.length = "none") + #
  geom_tree() +
  theme_tree() +
  geom_treescale(width	= 0.1) +
  scale_color_manual(values	= archColor,
                     na.value = "transparent",
                     guide = "none") +
  #geom_text2(aes(subset=!isTip,	label=node), hjust=-.3)+
  theme(legend.position = "right")
p1_arch

p2_arch <- p1_arch +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = archDat,
    pwidth	= 0.01,
    geom = geom_bar,
    mapping = aes(y = MAGs,	x = Abundance, fill = p),
    orientation = "y",
    stat = "identity",
  ) +
  scale_fill_manual(values	= archColor[-1]) +
  labs(fill = "Taxa") +
  new_scale_colour() +
  new_scale_fill()

p2_arch

p3_arch <- gheatmap(
  p2_arch,
  archDatset,
  width = 0.2,
  offset = 0.1,
  colnames = FALSE
) +
  scale_fill_discrete(na.translate = F) +
  labs(fill = "Datasets") +
  new_scale_colour() +
  new_scale_fill()
p3_arch

(
  p4_arch <- gheatmap(
    p3_arch,
    archcov,
    offset = 0.9,
    width = 0.05,
    colnames = FALSE
  ) +
    scale_fill_viridis_b() +
    labs(fill = "Normalized log10\nabundance")
)

ggsave_fitmax("Figures/Fig_S2_ArchaealTree.pdf",
              p4_arch,
              maxwidth = 9)
