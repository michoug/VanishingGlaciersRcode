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

dat_euk <- read_tsv("data/vMAGs_genomadEuk.txt")
dat_euk_vtax <-
  read_tsv("data/vMAGs_genomadEuk_tax.out", col_names = F)
dat_euk_ptax <- read_tsv("data/eMAGs_tax.tsv")

colnames(dat_euk_vtax) <-
  c("class", "seq_name", "taxId", "taxonomy")

dat_euk_ptax$MAGs <- gsub("-", "_", dat_euk_ptax$MAGs)

dat_virus_euk <- dat_euk %>%
  left_join(dat_euk_vtax,
            by = c("seq_name"),
            relationship = "many-to-many") %>%
  mutate(VTax = if_else(class == "U", taxonomy.x, taxonomy.y)) %>%
  select(MAGs, seq_name, VTax) %>%
  inner_join(dat_euk_ptax,by = join_by(MAGs)) %>%
  separate(
    VTax,
    into = c(
      "V_Virus",
      "V_Realm",
      "V_Kingdom",
      "V_Phylum",
      "V_Class",
      "V_Order",
      "V_Family",
      "V_Genus",
      "V_Species",
      "V_Name"
    ),
    sep = ";"
  ) %>%
  mutate(
    goodVirusTax = case_when(
      V_Class == "Caudoviricetes" &
        is.na(V_Order) ~ "Caudoviricetes Unclassified",
      V_Class == "Caudoviricetes" &
        V_Order == "unclassifiedCaudoviricetes" ~ "Caudoviricetes Unclassified",
      V_Class == "Caudoviricetes" ~ paste("Caudoviricetes", V_Order, sep = "_"),
      is.na(V_Class) ~ "UnclassifiedClass",
      .default = V_Class
    )
  ) %>%
  select(seq_name, goodVirusTax, MAGs, Tax)

dat_virus_euk$goodVirusTax <-
  gsub(" ", "", dat_virus_euk$goodVirusTax)

colnames(dat_virus_euk) <-
  c("seq_name", "Viral Taxonomy", "MAGs", "Eukaryotic Taxonomy")

dat_virus_sank <- dat_virus_euk %>%
  make_long(`Viral Taxonomy`, `Eukaryotic Taxonomy`) %>%
  mutate(
    node  = case_when(
      str_detect(node, "Caudoviricetes_") ~ gsub("_", " ", node),
      str_detect(node, "CaudoviricetesU") ~ gsub("sU", "s U", node),
      str_detect(node, "-") ~ gsub("-", " ", node, perl = T),
      str_detect(node, "UnclassifiedClass") ~ "Unclassified",
      .default = node
    )
  ) %>%
  mutate(next_node = if_else(
    str_detect(next_node, "-"),
    gsub("-", " ", next_node, perl = T),
    next_node
  ))

p1 <- ggplot(
  dat_virus_sank,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_text(size = 3, color = "black") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5, colour = "black"))

p1
ggsave_fitmax("Figures/Fig_5c_Sankey_Genomad_Euk.pdf", p1)
