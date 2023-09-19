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

library("RColorBrewer")
library("ggplot2")
library(tidyverse)

source("customFunctions/plot_functions.R")

datTax	<- read_tsv("../Eukaryotes/Phylofisher/Eukaryotic_MAGs_taxonomy.txt", col_names = T)

datCov	<- read_tsv("../Eukaryotes/Eukaryotes_MAGs_cov_norm.txt")
map	<- read_tsv("../metadata_NOMIS_sed_chem.txt")
busco <- read_tsv("../Eukaryotes/busco_30.txt")

map_sel <- map %>%
  select(Sample, Site_c)

datToPlot	<- datCov %>%
  pivot_longer(cols = !MAGs)%>%
  right_join(map_sel, by = join_by(name == Sample))%>%
  left_join(datTax, by = join_by(MAGs)) %>%
  filter(MAGs %in% busco$MAGs)%>%
  group_by(Lower_Taxonomy, Site_c) %>%
  summarise(sum	= sum(value)) %>%
  ungroup()%>%
  group_by(Site_c)%>%
  mutate(percentage	= sum / sum(sum))%>%
  mutate(Site_c = if_else(Site_c == "New_Zealand", "New Zealand", Site_c))%>%
  mutate(Lower_Taxonomy = gsub("_", " - ", Lower_Taxonomy))

datToPlot$Tax	<- as.factor(datToPlot$Lower_Taxonomy)
datToPlot$Tax	<-
  factor(datToPlot$Tax, levels = rev(levels(datToPlot$Tax)))
datToPlot$Tax	<-
  factor(datToPlot$Tax, levels = rev(levels(datToPlot$Tax)))

p1	<-
  ggplot(datToPlot, aes(x = Site_c, y = percentage*100, fill = Tax)) +
  geom_bar(stat	= "identity") +
  scale_fill_manual(name = "Most Abundant Taxa", values = colorRampPalette(brewer.pal(10, "Paired"))(length(unique(datToPlot$Tax)))) +
  theme_classic()+
  labs(y = "Percentage of Coverage (%)", x = NULL)+
  theme(
    panel.grid.major	= element_blank(),
    panel.grid.minor = element_blank(),
    panel.background	= element_blank(),
    legend.background = element_blank(),
    text	= element_text(colour = "black",size = 12),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    axis.line	= element_line(colour = "black")
  ) +
  theme(strip.background	= element_rect(fill = "transparent")) +
  theme(strip.text	= element_text(colour = 'black', size = 12))

p1
ggsave_fitmax("Figures/Fig_X_barplot_Eukaryotic_MAGS.pdf", p1, maxwidth = 13)


