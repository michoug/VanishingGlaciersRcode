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
library(tidyverse)
library(vegan)

source("customFunctions/plot_functions.R")

dat_host	<- read_tsv("../Virus/genomad_virus_host_cov_norm.txt.gz")
dat_euk <- read_tsv("../Virus/Genomad_Virus_Eukaryotes_cov_norm.txt.gz")
dat_tax_virus <- read_tsv("../Virus/nomis_high_complete_kaiju.names.out", col_names = F)
dat_cov_virus <- read_tsv("../Virus/phamb_viruses_depth_norm.txt.gz")
dat_virus_genomad <- read_tsv("../Virus/phamb_viruses_depth_genomad_norm.txt.gz")
map	<- read_tsv("../metadata_NOMIS.txt")
virusToRemove <- read_tsv("../Virus/virusToRemove.txt")

colnames(dat_tax_virus) <- c("class", "seq_name", "taxId", "taxonomy")

dat_virus_sel <- dat_cov_virus%>%
  pivot_longer(cols = !Contig)%>%
  right_join(dat_tax_virus, by = c("Contig" = "seq_name"))%>%
  left_join(dat_virus_genomad, by = c("Contig" = "Contig"))%>%
  mutate(taxonomy_final = if_else(
    class == "U", taxonomy.y, taxonomy.x
  ))%>%
  select(Contig, taxonomy_final,name,value)%>%
  distinct()%>%
  mutate(type = "Virus")

colnames(dat_virus_sel) <- c("seq_name","taxonomy","name","value", "type")


dat_host_red <- dat_host%>%
  select(!c(length,topology,coordinates,n_genes,
            genetic_code, virus_score,fdr,n_hallmarks,
            marker_enrichment,clean_seq_name, MAGs))%>%
  distinct()

dat_host_m <- dat_host_red %>%
  pivot_longer(cols = !c(seq_name, taxonomy))%>%
  mutate(type = "prok")%>%
  filter(!(seq_name %in% virusToRemove$Virus))

dat_euk_red <- dat_euk %>%
  select(!c(length,topology,coordinates,n_genes,
            genetic_code, virus_score,fdr,n_hallmarks,
            marker_enrichment,MAGs))%>%
  distinct()
  
dat_euk_m	<- dat_euk_red %>%
  pivot_longer(cols = !c(seq_name, taxonomy))%>%
  mutate(type = "euk")

datMerge	<- rbind(dat_virus_sel, dat_host_m, dat_euk_m)


dat_stats <- datMerge %>%
  select(-c(name, taxonomy, value))%>%
  distinct()%>%
  group_by(type)%>%
  summarise(n = n())



datMerge_tax <- datMerge %>%
  separate(taxonomy, into = c("Virus","Realm","Kingdom","Phylum","Class","Order","Family", "Genus", "Species","Name"), sep = ";")%>%
  filter(!(Realm == "unclassifiedviruses"))%>%
  filter(!(Kingdom == "unclassifiedRiboviria"))

datMergefinal <- datMerge_tax %>%
  mutate(goodTax = case_when(
    Class == "Caudoviricetes" & is.na(Order) ~ "Caudoviricetes_Unclassified",
    Class == "Caudoviricetes" & Order == "unclassifiedCaudoviricetes" ~ "Caudoviricetes_Unclassified",
    Class == "Caudoviricetes" ~ paste("Caudoviricetes", Order, sep = "_"),
    is.na(Class) ~ "UnclassifiedClass",
    .default = Class
    )
  )


mostAbundantTax	<- datMergefinal %>%
  group_by(goodTax) %>%
  summarise(total	= sum(value)) %>%
  arrange(desc(total)) %>%
  slice(1:9)


datToPlot	<- datMergefinal %>%
  mutate(Tax	= if_else(goodTax %in% mostAbundantTax$goodTax,  goodTax, "Other")) %>%
  group_by(type, Tax) %>%
  summarise(sum	= sum(value)) %>%
  na.omit()%>%
  mutate(percentage	= sum / sum(sum))%>%
  mutate(Tax = gsub("_", " - ", Tax))%>%
  mutate(Tax = if_else(Tax == "UnclassifiedClass", "Unclassified", Tax))%>%
  mutate(type = case_when(
    type == "euk" ~ "Eukaryotes\nAssociated",
    type == "prok" ~ "Prokaryotes\nAssociated",
    type == "Virus" ~ "Virus\nNon Associated"
  ))

datToPlot$Tax	<- as.factor(datToPlot$Tax)
datToPlot$Tax	<-
  factor(datToPlot$Tax, levels = rev(levels(datToPlot$Tax)))
datToPlot$Tax	<- relevel(datToPlot$Tax, "Other")
datToPlot$Tax	<-
  factor(datToPlot$Tax, levels = rev(levels(datToPlot$Tax)))

colorsName <- unique(datToPlot$Tax)
colors <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(datToPlot$Tax)))

p1	<-
  ggplot(datToPlot, aes(x = type, y = 100*percentage, fill = Tax)) +
  geom_bar(stat	= "identity") +
  scale_fill_manual(name = "Most Abundant Taxa", values = colors) +
  theme_classic()+
  theme(
    panel.grid.major	= element_blank(),
    panel.grid.minor = element_blank(),
    panel.background	= element_blank(),
    legend.background = element_blank(),
    text	= element_text(colour = "black", size = 12),
    #axis.line	= element_line(colour = "black")
  ) +
  theme(strip.background	= element_rect(fill = "transparent")) +
  theme(strip.text	= element_text(colour = 'black', size = 12))+
  labs(x = NULL, y = "Percentage of Coverage (%)")
  

p1
ggsave_fitmax("Figures/Fig_X_barplotVirusTax.pdf", p1)


