library(tidyverse)

tax <- read_tsv("data/pMAGs_tax.tsv")
check <- read_tsv("data/pMAGs_checkm.txt")
genes <- read_tsv("data/pMAGS_contigs_genes.txt.gz")
metab <- read_tsv("data/pMAGs_egg_KO.txt.gz")
coverage <- read_tsv("data/pMAGs_cov_norm.txt.gz")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

arch <- tax %>%
  filter(d == "d__Archaea")%>%
  left_join(check, join_by(MAGs))%>%
  left_join(metab, join_by(MAGs))%>%
  filter(KEGG_ko == "K00399")

arch_cont_genes <- metab %>%
  filter(Contigs %in% arch$Contigs)

arch_cov <- coverage %>%
  filter(MAGs %in% arch_cont_genes$MAGs)%>%
  as.data.frame()%>%
  column_to_rownames(var = "MAGs")

a <- t(arch_cov)
