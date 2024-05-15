library(tidyverse)
library("RColorBrewer")
library(ggrepel)
library(patchwork)

source("customFunctions/plot_functions.R")
specgen <- read_tsv("data/pMAGs_spec_gen.tsv")
tax <- read_tsv("data/pMAGs_tax.tsv")
covsum <- read_tsv("data/pMAGs_cov_sum.txt")

tax$p_c <- if_else(tax$p == "p__Proteobacteria", tax$c, tax$p)
tax$p_c <- gsub(".__","",tax$p_c, perl = T)
tax$p_c <- gsub("_.$","",tax$p_c, perl = T)

mostAbundantTax	<- tax%>%
  group_by(p_c)%>%
  summarise(total	= n())%>%
  arrange(desc(total))%>%
  slice(1:16)

list <- mostAbundantTax$p_c
list <- c(list, "Nitrospirota")

datToPlot	<- tax%>%
  right_join(specgen,join_by(MAGs))%>%
  filter(sign != "NON SIGNIFICANT") %>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))%>%
  group_by(p_c, sign)%>%
  summarise(sum	= n())%>%
  group_by(sign)%>%
  mutate(totsum = sum(sum))%>%
  mutate(fraction = sum / totsum)%>%
  mutate(ymax = cumsum(fraction))%>%
  mutate(ymin = c(0, head(ymax, n=-1)))%>%
  mutate(labelPosition = (ymax +ymin) / 2)


p1 <- ggplot(datToPlot, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=p_c)) +
  geom_rect() +
  geom_label_repel( x=4, aes(y=labelPosition, label=p_c), size=4, max.overlaps = 15) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(10,"Paired"))(length(unique(datToPlot$p_c)))) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  facet_wrap(~sign, drop = TRUE)

# p1

datToPlotCov	<- tax%>%
  right_join(specgen,join_by(MAGs))%>%
  left_join(covsum)%>%
  filter(sign != "NON SIGNIFICANT") %>%
  mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))%>%
  group_by(p_c, sign)%>%
  summarise(sum	= sum(percentage))%>%
  group_by(sign)%>%
  mutate(totsum = sum(sum))%>%
  mutate(fraction = sum / totsum)%>%
  mutate(ymax = cumsum(fraction))%>%
  mutate(ymin = c(0, head(ymax, n=-1)))%>%
  mutate(labelPosition = (ymax +ymin) / 2)

p2 <- ggplot(datToPlotCov, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=p_c)) +
  geom_rect() +
  geom_label_repel( x=4, aes(y=labelPosition, label=p_c), size=4, max.overlaps = 15) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(10,"Paired"))(length(unique(datToPlotCov$p_c)))) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  facet_wrap(~sign, drop = TRUE)

p2

pfin <- wrap_plots(p1,p2, ncol = 1) +
  plot_annotation(tag_levels = "a")
pfin

ggsave_fitmax("NewAnalyses/specgen_abundance_coverage.pdf", pfin, maxwidth = 15)
