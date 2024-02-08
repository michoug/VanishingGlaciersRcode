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

BacTree	<- read.tree("data/pMAGs_bact_gtdtk_midroot.tree")
dat <- read_tsv("data/pMAGS_tax.tsv")
dataset <- read_tsv("data/pMAGS_Presence_Datasets.txt")
covmax <- read_tsv("data/pMAGs_cov_sum.txt")

dat$p_c <- if_else(dat$p == "p__Proteobacteria", dat$c, dat$p)
dat$p_c <- gsub(".__", "", dat$p_c, perl = T)
dat$p_c <- gsub("_.$", "", dat$p_c, perl = T)

bacDat <- dat %>%
  filter(d == "d__Bacteria") %>%
  mutate(Abundance	= 0.2)

mostAbundantTax	<- bacDat %>%
  group_by(p_c) %>%
  summarise(total	= n()) %>%
  arrange(desc(total)) %>%
  slice(1:16)

list <- mostAbundantTax$p_c
list <- c(list, "Nitrospirota")


bacDat	<- bacDat %>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))

bacDatset <- dataset %>%
  filter(MAGs %in% bacDat$MAGs)

bacDatset$Busi <-
  gsub("Busi", "Busi et al., Nat Com, 2022", bacDatset$Busi)
bacDatset$ENSEMBLE <-
  gsub("ENSEMBLE", "Michoud et al, L&O, 2023", bacDatset$ENSEMBLE)
bacDatset$Tibet <-
  gsub("Tibet", "Tibetan Glacier Genome and Gene", bacDatset$Tibet)
bacDatset$Tara <- gsub("Tara", "Tara Oceans", bacDatset$Tara)

bacDatset <- as.data.frame(bacDatset)
rownames(bacDatset) <- bacDatset$MAGs
bacDatset$MAGs <- NULL


bactcov <- covmax %>%
  filter(MAGs %in% bacDat$MAGs) %>%
  select(MAGs, count)

bactcov <- as.data.frame(bactcov)
rownames(bactcov) <- bactcov$MAGs
bactcov$MAGs <- NULL

bactcov$count <- log10(bactcov$count)

a	<- split(bacDat$MAGs, bacDat$p_c)

tree	<- groupOTU(BacTree, a)

getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))

bactColor	<- getPaletteBact(length(unique(bacDat$p_c)) + 1)

bactColor[1]	<- "black"


p1	<-
  ggtree(tree,
         layout = 'circular',
         aes(color = group)) + #
  geom_tree() +
  theme_tree() +
  geom_treescale(width	= 0.1) +
  scale_color_manual(values	= bactColor,
                     na.value = "transparent",
                     guide = "none") +
  #geom_text2(aes(subset=!isTip,	label=node), hjust=-.3)+
  theme(legend.position = "right")
# p1

p2 <- p1 +
  new_scale_colour() +
  new_scale_fill() +
  geom_fruit(
    data = bacDat,
    pwidth	= 0.01,
    geom = geom_bar,
    mapping = aes(y = MAGs, fill = p_c, x = 1),
    # orientation="y",
    stat = "identity",
  ) +
  scale_fill_manual(values	= bactColor[-1]) +
  labs(fill = "Taxa") +
  new_scale_colour() +
  new_scale_fill()

# p2

p3 <-
  gheatmap(
    p2,
    bacDatset,
    width = 0.2,
    offset = 0.1,
    # offset=8, width=0.6,
    colnames = FALSE,
    color = NULL
  ) +
  scale_fill_discrete(na.translate = F) +
  labs(fill = "Datasets") +
  new_scale_colour() +
  new_scale_fill()
# p3

p4 <- gheatmap(
  p3,
  bactcov,
  offset = 0.6,
  width = 0.05,
  colnames = FALSE,
  color = NULL
) +
  scale_fill_viridis_c() +
  labs(fill = "Normalized log10\nabundance")

# p4

ggsave_fitmax("Figures/Fig_1c_bacterialTree.pdf",
              p4,
              maxheight =  15)
