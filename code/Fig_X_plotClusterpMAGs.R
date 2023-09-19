## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-19
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(ggnewscale)
library(RColorBrewer)

hc <- readRDS("../Prokaryotes/Cluster/agnes_cluster_pMAGs.rds")
dat <- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv")
specgen <- read_tsv("../Prokaryotes/EcolUtils/spec_gen_new.tsv")
max <- read_tsv("../Prokaryotes/MAGs_max_cov.txt")

specgen_sel <- specgen%>%
  select(MAGs, sign)

specgen_sel <- as.data.frame(specgen_sel)
rownames(specgen_sel) <- specgen_sel$MAGs
specgen_sel$MAGs <- NULL

dat$p_c <- if_else(dat$p == "p__Proteobacteria", dat$c, dat$p)
dat$p_c <- gsub(".__","",dat$p_c, perl = T)
dat$p_c <- gsub("_.$","",dat$p_c, perl = T)

bacDat <- dat %>%
  filter(d == "d__Bacteria")

mostAbundantTax	<- bacDat%>%
  group_by(p_c)%>%
  summarise(total	= n())%>%
  arrange(desc(total))%>%
  slice(1:10)

bacDat	<- bacDat%>%
  mutate(p_c	= if_else(p_c %in% mostAbundantTax$p_c,  p_c, "Other"))%>%
  select(MAGs, p_c)

bacDat <- as.data.frame(bacDat)
rownames(bacDat) <- bacDat$MAGs
bacDat$MAGs <- NULL

getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))
bactColor	<- getPaletteBact(length(unique(bacDat$p_c)))

clus <- cutree(as.hclust(hc), k = 6)
g <- split(names(clus), clus)

p <- ggtree(hc)
clades <- sapply(g, function(n) MRCA(p, n))

names(clades) <- c("Diverse Taxa","Myxcoccota &\n Bdellovibrionata",
    "Alpha and \n Gammaproteobacteria",
    "Bacteroidia",
    "Pactescibacteria",
    "Plantomycetes")

p <- groupClade(p, clades, group_name='subtree') + 
  aes(color=subtree)

colors <- c(brewer.pal(6, "Set1"), "black")

p1 <- p +
  layout_circular() + 
  scale_color_manual(values = colors, breaks = names(clades)) +
  labs(color = "Functional Cluster")+
  theme_dendrogram(plot.margin=margin(6,6,80,6)) +
  theme(legend.position=c(.9, .6))+
  new_scale_colour()+
  new_scale_fill()

p2 <- gheatmap(p1, bacDat,width = 0.05,offset = 0.1,# offset=8, width=0.6, 
               colnames=FALSE, color = NULL)+
  scale_fill_manual(values	= bactColor)+
  scale_x_ggtree()+
  labs(fill = "Taxa")
# p2

ggsave("Figures/Fig_X_bactClustTree.pdf",p2,width = 18,height = 18)
