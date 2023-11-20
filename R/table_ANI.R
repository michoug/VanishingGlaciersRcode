## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-11-20
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
library(purrr)

dat_Tara <- read_tsv("../Prokaryotes/MAGs_Comparison/Tara_NOMIS.txt", col_names = F)
dat_Busi <- read_tsv("../Prokaryotes/MAGs_Comparison/Busi_NOMIS.txt", col_names = F)
dat_Ensemble <- read_tsv("../Prokaryotes/MAGs_Comparison/NOMIS_ENSEMBLE.txt", col_names = F)
dat_TG2G <- read_tsv("../Prokaryotes/MAGs_Comparison/TG2G_NOMIS.txt", col_names = F)
mapping <- read_tsv("../metadata_NOMIS.txt")
maxcov <- read_tsv("../Prokaryotes/MAGs_max_cov.txt")
magsinRocks <- read_tsv("../Prokaryotes/MAGsInRocks.txt")


filter_table <- function(dat){

  colnames(dat) <- c("NOMIS", "Data", "ANI", "NOMIS_cont", "Data_cont")
  
  dat_mod <- dat %>%
    filter(ANI >= 95)%>%
    mutate(NOMIS = gsub(".*/", "", NOMIS))%>%
    mutate(NOMIS = gsub(".fa.gz", "", NOMIS))%>%
    mutate(NOMIS = gsub(".fa", "", NOMIS))%>%
    separate(Data, into = c("dataset", "dataset_MAGs"), sep = "/")
  
  return(dat_mod)
}


df_temp <- map_df(list(dat_Tara, dat_Busi, dat_Ensemble, dat_TG2G), filter_table)

final_df <- df_temp %>%
  # filter(!(NOMIS %in% magsinRocks$MAGs))%>%
  mutate(dataset = if_else(dataset == "genomes-fasta", "Tara_MAGs", dataset))%>%
  select(NOMIS, dataset)%>%
  distinct()%>%
  left_join(maxcov, join_by(NOMIS == MAGs))%>%
  left_join(mapping, join_by(row_max == Sample))%>%
  select(NOMIS, dataset, Site_c)%>%
  group_by(dataset, Site_c)%>%
  summarise(n = n())%>%
  na.omit()


