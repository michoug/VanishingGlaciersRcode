## ---------------------------
##
## Script name: donutPlot.R
##
## Purpose of script: Create donuts plots based on the the number and abundance (coverage) of MAGs in the
## Vanishing glacier dataset
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-05
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
library("ggrepel")

source("customFunctions/plot_functions.R")

dat	<- read_tsv("data/pMAGs_tax.tsv")
covsum <- read_tsv("data/pMAGs_cov_sum.txt")


dat$p_c <- if_else(dat$p == "p__Proteobacteria", dat$c, dat$p)
dat$p_c <- gsub(".__","",dat$p_c, perl = T)
dat$p_c <- gsub("_.$","",dat$p_c, perl = T)

mostAbundantTax	<- dat%>%
  group_by(p_c)%>%
  summarise(total	= n())%>%
  arrange(desc(total))%>%
  slice(1:16)

list <- mostAbundantTax$p_c
list <- c(list, "Nitrospirota")

datToPlot	<- dat%>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))%>%
  group_by(p_c)%>%
  summarise(sum	= n())

datToPlot$p_c <- gsub(".__","",datToPlot$p_c)

# Compute percentages
datToPlot$fraction = datToPlot$sum / sum(datToPlot$sum)

# Compute the cumulative percentages (top of each rectangle)
datToPlot$ymax = cumsum(datToPlot$fraction)

# Compute the bottom of each rectangle
datToPlot$ymin = c(0, head(datToPlot$ymax, n=-1))

datToPlot$labelPosition <- (datToPlot$ymax + datToPlot$ymin) / 2


p1 <- ggplot(datToPlot, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=p_c)) +
  geom_rect() +
  geom_label_repel( x=4, aes(y=labelPosition, label=p_c), size=6) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(10,"Paired"))(length(unique(datToPlot$p_c)))) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")


ggsave_fitmax("Figures/Fig_SX_donutPlotAbundance.pdf", p1, maxwidth = 10)


datToPlotCov	<- dat%>%
  left_join(covsum)%>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))%>%
  group_by(p_c)%>%
  summarise(sum	= sum(percentage))

datToPlotCov$p_c <- gsub(".__","",datToPlotCov$p_c)

# Compute percentages
datToPlotCov$fraction = datToPlotCov$sum / sum(datToPlotCov$sum)

# Compute the cumulative percentages (top of each rectangle)
datToPlotCov$ymax = cumsum(datToPlotCov$fraction)

# Compute the bottom of each rectangle
datToPlotCov$ymin = c(0, head(datToPlotCov$ymax, n=-1))

datToPlotCov$labelPosition <- (datToPlotCov$ymax + datToPlotCov$ymin) / 2


p2 <- ggplot(datToPlotCov, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=p_c)) +
  geom_rect() +
  geom_label_repel( x=4, aes(y=labelPosition, label=p_c), size=6) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(10,"Paired"))(length(unique(datToPlotCov$p_c)))) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p2

ggsave_fitmax("Figures/Fig_SX_donutPlotCoverage.pdf", p2, maxwidth = 10)

