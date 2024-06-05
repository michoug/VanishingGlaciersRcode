library(tidyverse)

dat_egg <- read_tsv("data/pMAGs_egg_KO.txt.gz")
carot <- read_tsv("data/Other/carotenoids_genes.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")
cov <- read_tsv("data/pMAGs_cov_sum.txt")

dat_clean <- dat_egg %>%
  right_join(carot, join_by(KEGG_ko))%>%
  select(MAGs, Pathway, KEGG_ko)%>%
  left_join(tax, join_by(MAGs))%>%
  left_join(cov, join_by(MAGs))%>%
  group_by(Pathway, c)%>%
  summarise(sum = sum(percentage))
  
