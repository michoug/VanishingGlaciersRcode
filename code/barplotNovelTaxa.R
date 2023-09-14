## ---------------------------
##
## Script name: barplotNovelTaxa.R
##
## Purpose of script: Create a barplot showing for each taxonomical level, the number of known and unknows taxa
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-14
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
library(ggplot2)
library(reshape2)
library(scales)

dat	<- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv", col_names = T)

dat_m <- melt(dat, id.vars = "MAGs")
dat_m$value <- gsub(".__","", dat_m$value)
dat_m$status <- if_else(dat_m$value == "", "UnKnown", "Known")

dat2Plot <- dat_m %>%
  select(-MAGs)%>%
  group_by(variable, status)%>%
  summarise(sum = n(), perc = sum/2868, label = percent(perc, accuracy = 0.1))%>%
  filter(!(variable %in% c("d","p","s")))%>% # Remove taxa levels if there is no unknown
  mutate(label = if_else(status == "UnKnown",paste(sum, " MAGs (", label, ")", sep = ""),""))%>%
  mutate(variable = case_when(
    variable == "c" ~ "Class",
    variable == "o" ~ "Order",
    variable == "g" ~ "Genus",
    variable == "f" ~ "Family"
    # .default = "other".default = "other"
  ))
  
dat2Plot$variable <- factor(dat2Plot$variable, levels = c("Class", "Order", "Family", "Genus"))

ggplot(dat2Plot, aes(x = perc, y = reorder(variable, desc(variable)), fill = status, label = label))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(nudge_x = 0.25)+
  scale_fill_manual(values = c("Known" = "white", "UnKnown" = "black"), labels = c("Known taxa", "Novel taxa"))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  xlim(0,1)+
  labs(fill = NULL, x = "Percentage of MAGs (%)", y = NULL)+
  theme_classic()+
  theme(legend.position="top")

ggsave("Figures/barplotNovelTaxa.pdf")
