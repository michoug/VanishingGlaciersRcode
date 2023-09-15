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

BacTree	<- read.tree("../Prokaryotes/TaxQual/NOMIS_MAGS_gtdbtk_full/classify/gtdbtk.bac120.classify_user_midpoint.tree")
dat <- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv")
dataset <- read_tsv("../Prokaryotes/TaxQual/MAGS_Presence_Datasets.txt")
covmax <- read_tsv("../Prokaryotes/MAGs_cov_sum.txt")

dat$p_c <- if_else(dat$p == "p__Proteobacteria", dat$c, dat$p)
dat$p_c <- gsub(".__","",dat$p_c, perl = T)
dat$p_c <- gsub("_.$","",dat$p_c, perl = T)

bacDat <- dat %>%
  filter(d == "d__Bacteria")

bacDatset <- dataset%>%
  filter(MAGs %in% bacDat$MAGs)#%>%
# left_join(covmax)%>%
# select(MAGs, percentage, Busi:Tara)

mostAbundantTax	<- bacDat%>%
  group_by(p_c)%>%
  summarise(total	= n())%>%
  arrange(desc(total))%>%
  slice(1:16)

list <- mostAbundantTax$p_c
list <- c(list, "Nitrospirota")


bacDat	<- bacDat%>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))

bactcov <- covmax%>%
  filter(MAGs %in% bacDat$MAGs)%>%
  select(MAGs, count)

bacDatset <- as.data.frame(bacDatset)
rownames(bacDatset) <- bacDatset$MAGs
bacDatset$MAGs <- NULL

bactcov <- as.data.frame(bactcov)
rownames(bactcov) <- bactcov$MAGs
bactcov$MAGs <- NULL

bactcov$count <- log10(bactcov$count)

a	<- split(bacDat$MAGs, bacDat$p_c)

bacDat$Abundance	<- 0.2

tree	<- groupOTU(BacTree, a)

getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))

bactColor	<- getPaletteBact(length(unique(bacDat$p_c)) + 1)

bactColor[1]	<- "black"

#bactColor	<- c("red", "blue", bactColor)

p1	<-
  ggtree(tree,
         layout = 'circular',
         aes(color = group),
         #branch.length = "none"
         ) + #
  geom_tree() +
  theme_tree() +
  geom_treescale(width	= 0.1) +
  scale_color_manual(values	= bactColor,
                     na.value = "transparent",
                     guide = "none") +
  #geom_text2(aes(subset=!isTip,	label=node), hjust=-.3)+
  theme(legend.position = "right", legend.title = element_blank())
# p1

p2 <- p1 +
  new_scale_colour()+
  new_scale_fill()+
  geom_fruit(
    data=bacDat,
    pwidth	= 0.01,
    geom=geom_bar,
    mapping=aes(y=MAGs,	x = Abundance, fill = p_c),
    orientation="y",
    stat="identity",
  )+
  scale_fill_manual(values	= bactColor[-1])+
  new_scale_colour()+
  new_scale_fill()

# p2

p3 <- gheatmap(p2, bacDatset,width = 0.2,offset = 0.1,# offset=8, width=0.6, 
               colnames=FALSE) +
  new_scale_colour()+
  new_scale_fill()
# p3

p4 <- gheatmap(p3, bactcov, offset=0.9, width=0.05,
               colnames=FALSE)+
  scale_fill_viridis_b()

# p4

# bac <- print(open_tree(p4, angle = 2))

ggsave_fitmax("Figures/bacterialTree.pdf",
       p4,
       maxwidth = 18)
