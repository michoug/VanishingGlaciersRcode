## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2024-01-16
##
## Copyright (c) GrÃ©goire Michoud, 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(tidyverse)


dat_egg <- read_tsv("../data/pMAGs_egg_KO.txt.gz")
uv_genes <- read_tsv("../../Prokaryotes/UV_resistance/DNA_repair.txt")
uv <- uv_genes %>%
  group_by(Pathway)%>%
  mutate(nb_genes = n())

dat_select <- dat_egg %>%
  right_join(uv, by = join_by(KEGG_ko),relationship = "many-to-many")%>%
  select(-c(`#query`, Contigs))%>%
  group_by(MAGs, Pathway, KEGG_ko) %>%
  distinct()%>%
  mutate(n= n())%>%
  mutate(n = if_else(n > 1, 1, n))%>%
  ungroup()%>%
  group_by(MAGs, Pathway, nb_genes)%>%
  summarise(nb = sum(n))
