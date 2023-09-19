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
library(ggsankey)

source("customFunctions/plot_functions.R")

dat_prok <- read_tsv("../Virus/Genomad_Virus_host/NOMIS_genomad_virus_summary.txt")
dat_prok_vtax <- read_tsv("../Virus/Genomad_Virus_host/allhostviruses.names.out", col_names = F)
dat_prok_ptax <- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv")

colnames(dat_prok_vtax) <- c("class", "seq_name", "taxId", "taxonomy")

dat_prok_vtax$taxonomy <- gsub(" ", "", dat_prok_vtax$taxonomy)
dat_prok_vtax$taxonomy <- gsub(";$", "", dat_prok_vtax$taxonomy)

dat_virus_prok <- dat_prok%>%
  left_join(dat_prok_vtax,by = c("seq_name"),relationship = "many-to-many")%>%
  mutate(VTax = if_else(
    class == "U", taxonomy.x, taxonomy.y
  ))%>%
  select(MAGs, seq_name,VTax)%>%
  inner_join(dat_prok_ptax)%>%
  separate(VTax, into = c("V_Virus","V_Realm","V_Kingdom","V_Phylum","V_Class","V_Order","V_Family", "V_Genus", "V_Species","V_Name"), sep = ";")%>%
  mutate(goodVirusTax = case_when(
    V_Class == "Caudoviricetes" & is.na(V_Order) ~ "Caudoviricetes_Unclassified",
    V_Class == "Caudoviricetes" & V_Order == "unclassifiedCaudoviricetes" ~ "Caudoviricetes_Unclassified",
    V_Class == "Caudoviricetes" ~ paste("Caudoviricetes", V_Order, sep = "_"),
    is.na(V_Class) ~ "UnclassifiedClass",
    .default = V_Class
  ))%>%
  mutate(taxa_good = if_else(p == "p__Proteobacteria", c, p))%>%
  select(seq_name, goodVirusTax, MAGs, taxa_good)



mostabundantVir <- dat_virus_prok %>%
  select(-MAGs)%>%
  distinct()%>%
  group_by(goodVirusTax)%>%
  summarise(n = n())%>%
  arrange(desc(n)) %>%
  slice(1:15)

# Find the most abundant prokaryotic class

mostabundantProk <- dat_virus_prok %>%
  select(-MAGs)%>%
  distinct()%>%
  group_by(taxa_good)%>%
  summarise(n = n())%>%
  arrange(desc(n)) %>%
  slice(1:15)

# Filter the data to get the most abundant Viral and Prokaryotic classes and remove duplicate viruses that are hosted by bacteria of
# the same taxonomy

dat_virus_prok_sel <- dat_virus_prok %>%
  filter(goodVirusTax %in% mostabundantVir$goodVirusTax)%>%
  filter(taxa_good %in% mostabundantProk$taxa_good)%>%
  select(-MAGs)%>%
  distinct()

colnames(dat_virus_prok_sel) <- c("seq_name", "Viral Taxonomy", "Prokaryotic Taxonomy")

dat_virus_sank_p <- dat_virus_prok_sel %>%
  make_long(`Viral Taxonomy`, `Prokaryotic Taxonomy`)%>%
  mutate(node  = case_when(
    str_detect(node, "Caudoviricetes_") ~ gsub("_", " ", node),
    str_detect(node, "__") ~ gsub(".__", "", node, perl = T),
    .default = node
  ))%>%
  mutate(next_node = if_else(str_detect(next_node, "__"),gsub(".__", "", next_node, perl = T), next_node ))

p1 <- ggplot(dat_virus_sank_p, aes(x = x, 
                                 next_x = next_x, 
                                 node = node, 
                                 next_node = next_node,
                                 fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_text(size = 3, color = "black") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5, colour = "black")) 

p1
ggsave_fitmax("Figures/Fig_X_Sankey_Genomad_Prok.pdf", p1)



