library(tidyverse)

dat_egg <- read_tsv("data/pMAGs_egg_KO.txt.gz")
carot <- read_tsv("data/Other/carotenoids_genes.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")
cov <- read_tsv("data/pMAGs_cov_sum.txt")

dat_clean <- dat_egg %>%
  right_join(carot, join_by(KEGG_ko))%>%
  select(MAGs, Pathway, KEGG_ko)%>%
  left_join(tax, join_by(MAGs))%>%
  group_by(Pathway, KEGG_ko, c)%>%
  mutate(n = n())%>%
  ungroup()%>%
  left_join(cov, join_by(MAGs))%>%
  group_by(Pathway, c)%>%
  mutate(sum = sum(percentage))%>%
  select(Pathway, c, n, sum)%>%
  distinct()%>%
  mutate(Percentage = sum *100)%>%
  select(-sum)

write_tsv(dat_clean, "Table_SX_carotenoid_genes.txt")  
