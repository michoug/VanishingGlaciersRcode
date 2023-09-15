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

library("ggtree")
library("ape")
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)

source("customFunctions/plot_functions.R")

tree	<- read.tree("../Eukaryotes/Phylofisher/Above30perc/tree_fast.tree")
tax_mags <- read_tsv("../Eukaryotes/Phylofisher/Above30perc/metadata_input.tsv") 
tax_ref <- read_tsv("../Eukaryotes/Phylofisher/Above30perc/metadata.tsv")
cov <- read_tsv("../Eukaryotes/Eukaryotes_MAGs_cov_norm.txt")

tax_mags <- tax_mags %>%
  select(`Unique ID`, `Lower Taxonomy`)%>%
  mutate(type = if_else(grepl("^GC",`Unique ID`),"ref","mags"))

taxToPlot <- tax_mags%>%
  filter(type == "mags")%>%
  filter(`Unique ID`%in% tree$tip.label)

tax_ref <- tax_ref%>%
  select(`Unique ID`, `Lower Taxonomy`)

tax_ref$type <- "ref"

tax <- rbind(tax_mags, tax_ref)
colnames(tax) <- c("MAGs","Tax", "type")

tax <- tax %>%
  mutate(taxa_good = if_else(Tax %in% taxToPlot$`Lower Taxonomy`, Tax, "Other"))

meta	<- split(tax$MAGs, tax$taxa_good)

tree_meta	<- groupOTU(tree, meta)

getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))

# color	<- getPaletteBact(length(unique(tax$taxa_good)) + 1)

color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
           "#FF7F00", "#FFFF33", "#A65628", "black", "#999999")

color[1]	<- "black"

p1	<-
  ggtree(tree_meta,
         layout = 'circular',
         aes(color = group),
         #branch.length = "none"
         ) + #
  geom_tree() +
  theme_tree() +
  # geom_tiplab()+
  geom_treescale(width	= 0.1) +
  scale_color_manual(values	= color,
                     na.value = "transparent",
                     guide = "none") +
  theme(legend.position = "right", legend.title = element_blank())+
  new_scale_colour()+
  new_scale_fill()

tax_temp <- tax %>%
  filter(MAGs %in% tree$tip.label)
# 
# tax_comb <- tax_temp %>%
#   select(MAGs,Tax)
# 
# cov$MAGs <- gsub("_","-",cov$MAGs)
# cov_max <- cov %>%
#   # filter(MAGs %in% tree$tip.label)%>%
#   rowwise()%>%
#   summarise(total = sum(c_across(where(is.numeric))))
# 
# 
# cov_max$MAGs <- cov$MAGs
# 
# cov_max_all <- cov_max %>%
#   filter(MAGs %in% tree$tip.label)%>%
#   right_join(tax_comb, multiple = "all")%>%
#   distinct()%>%
#   select(MAGs, total)%>%
#   mutate_all(~replace_na(.,1))
# 
# cov_max_all <- as.data.frame(cov_max_all)
# 
# rownames(cov_max_all) <- cov_max_all$MAGs
# cov_max_all$MAGs <- NULL
# cov_max_all$count <- log10(cov_max_all$total)
# cov_max_all$total <- NULL

#write_tsv(tax_temp, "metadata_input_above30.tsv")

p2 <- p1 +
  new_scale_colour()+
  new_scale_fill()+
  geom_fruit(
    data=tax,
    pwidth	= 0.05,
    geom=geom_bar,
    mapping=aes(y=MAGs, x = 4, fill = taxa_good),
    orientation="y",
    stat="identity"
  )+
  scale_fill_manual(values	= color[-1])

p2

p3 <- p2 %<+% tax_temp + geom_tippoint(aes(color = type))#+geom_tiplab()
p3

# p4 <- gheatmap(p3, cov_max_all, offset=0.9, width=0.05,
#                     colnames=FALSE)#+
#   scale_fill_viridis_b()
# p4

ggsave_fitmax("Figures/EukaryoticTreeAllAbove30_ochrophyte_name.pdf",
       p3,
       maxwidth = 10)
