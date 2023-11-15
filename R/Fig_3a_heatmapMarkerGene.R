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

library("ggplot2")
library(tidyverse)
library(scales)
library(ggpubr)
library(RColorBrewer)

dat <- read_tsv("../Prokaryotes/EggNog/listMetabolismMAGs.txt.gz")
taxAll <- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv")
marker_list <- read_tsv("../Prokaryotes/EggNog/MarkerGenesReduced.txt")

marker_list$Gene_name <- factor(marker_list$Gene_name, levels = unique(marker_list$Gene_name))


taxAll_t <- taxAll %>%
  group_by(c)%>%
  summarise(n = n())%>%
  filter(n >=10)

taxAll <- merge(taxAll, taxAll_t, by = "c")
taxAll$c <- paste(taxAll$c, " (",taxAll$n,")", sep = "")

taxAll$c <-  gsub("c__", "", taxAll$c)

dat_merge <- dat%>%
  distinct()%>%
  mutate(value = 1)%>%
  select(-c(Gene_name, Gene_description, Metabolism, Pathway))%>%
  complete(MAGs, KO_ID, fill = list(value = 0))%>%
  left_join(marker_list,by = join_by(KO_ID)) %>%
  left_join(taxAll, by = join_by(MAGs))%>%
  select(-n)%>%
  filter(!(KO_ID == "K14469"))

tax <- dat_merge %>%
  select(p, c)%>%
  distinct()

tax$c <- factor(tax$c)

tax$c <- fct_reorder(tax$c, tax$p)


dat_Pivot<- dat_merge%>%
  group_by(MAGs, KO_ID, Gene_name,Gene_description, Pathway, Pathway_small,Metabolism,c)%>%
  summarise(n = sum(value))%>%
  ungroup()%>%
  mutate(
    number = if_else(n > 1, 1, n)
  )%>%
  group_by(Gene_name, Pathway_small, Metabolism, c)%>%
  summarise(mean_tax = mean(number)*100)%>%
  mutate(c = gsub("c__", "", c))


dat_Pivot$c <- factor(dat_Pivot$c, levels = levels(tax$c))

list <- unique(dat_Pivot$Metabolism)
list <- as.data.frame(list)
list$color <- brewer.pal(length(unique(list$list)), "Dark2")
plot_list <- list()

d <- list$list[6]

for (d in list$list) {
  dat_lol <- dat_Pivot%>%
    filter(Metabolism == !!d)
  
  dat_sum <- dat_lol%>%
    group_by(Gene_name)%>%
    summarise(sum = sum(mean_tax))%>%
    filter(sum == 0)
  
  dat_lol <- dat_lol%>%
    filter(!(Gene_name %in% dat_sum$Gene_name))%>%
    na.omit()
  
  list_temp <- list%>%
    filter(list == !!d)
  
  (p1 <- ggplot(dat_lol , aes(Gene_name,  reorder(c, desc(c)))) + #reorder(variable, desc(variable)))
      geom_tile(aes(fill = mean_tax),colour = "black") +
      scale_fill_gradient(low = "white",high = list_temp$color, limits = c(0,100))+
      #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
      theme_minimal()+
      coord_fixed(ratio=1)+
      labs(x = NULL, y = NULL, fill = "Propotion of MAGs (%)")+
      theme(text = element_text(size = 12, colour = "black"),axis.text.x = element_text(angle =45, hjust = 1))
  )
  name = d
  name <- gsub(" ", "_", name)
  name <- gsub("/", "_", name)
  name <- paste("Figures/temp/Fig_3a_",name, ".pdf", sep = "")
  ggsave(name, p1, height = 7)
  plot_list[[d]] = p1
}

# ggarrange(plotlist=plot_list, common.legend = T, align = "h")
