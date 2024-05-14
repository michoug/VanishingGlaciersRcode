## ---------------------------
##
## Script name: 
##
## Purpose of script: Plot different parameters describing MAGs
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2024-04-24
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(ggplot2)
library(tidyverse)
library(ggpubr)

source("customFunctions/plot_functions.R")

dat_mags	<- read_tsv("data/pMAGs_checkm.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")
dat_gtdb <- read_tsv("../bac120_metadata_r220.tsv.gz")

dat_family_MAGs <- dat_mags %>%
  left_join(tax,join_by(MAGs))%>%
  arrange(desc(LengthNorm2))%>%
  head(100)%>%
  group_by(f)%>%
  mutate(n = n())%>%
  # summarise(mean = mean(LengthNorm2), sd = sd(LengthNorm2))%>%
  filter(f != "f__")%>%
  arrange(desc(n))%>%
  filter(n > 6)%>%
  select(f, LengthNorm2)%>%
  mutate(type = "MAGs")

dat_class_MAGs <- dat_mags %>%
  left_join(tax,join_by(MAGs))%>%
  arrange(desc(LengthNorm2))%>%
  head(100)%>%
  filter(o == "o__")%>%
  select(c, LengthNorm2)%>%
  mutate(type = "MAGs")

gtdb_select <- dat_gtdb %>%
  select(gtdb_taxonomy, checkm2_completeness, checkm2_contamination, genome_size)%>%
  separate(gtdb_taxonomy, sep = ";", into = c("d","p","c","o","f","g"))%>%
  filter(f %in% dat_family_MAGs$f)%>%
  filter(checkm2_completeness != "none")%>%
  mutate(LengthNorm2 = genome_size * (100/as.numeric(checkm2_completeness)))%>%
  select(f, LengthNorm2)%>%
  mutate(type = "GTDB")

dat_plot <- dat_family_MAGs%>%
  bind_rows(gtdb_select)%>%
  mutate(LengthNorm2 = LengthNorm2 / 1e6)%>%
  rename(Type = type)%>%
  mutate(f = gsub("f__", "",f))%>%
  arrange(f)
  # ggplot(aes(x = f, y = LengthNorm2, fill = type))+
  # geom_violin()


(p <- ggboxplot(dat_plot, x = "f", y = "LengthNorm2",
          color = "Type",
          add = "jitter", shape = "Type", ggtheme = theme_pubr())+
  labs(y = "Normalized Genome Length (Mbp)", x = "Family")+
  theme(text = element_text(colour = "black", size = 12)))
  
ggsave_fitmax("Figures/Fig_X_LargestGenomesFamilies.pdf",p,maxwidth = 10)  

