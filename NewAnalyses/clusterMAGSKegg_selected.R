library(tidyverse)
library(ggtree)
library(cluster)


dat <- read_tsv("data/pMAGs_egg_KO.txt.gz")
tree <- read.tree("data/pMAGs_bact_gtdtk_midroot.tree")
marker_list <-
  read_tsv("data/Other/DefKeggPathway.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")

marker_met <- marker_list %>%
  filter(Path_1 == "09100 Metabolism")
  
dat_mat <- dat%>%
  select(-c(Contigs, `#query`)) %>%
  filter(KEGG_ko %in% marker_met$ko)%>%
  filter(MAGs %in% tree$tip.label)%>%
  add_count(KEGG_ko, MAGs)%>%
  complete(KEGG_ko, MAGs, fill=list(n=0))%>%
  distinct()%>%
  pivot_wider(names_from = KEGG_ko, values_from = n)%>%
  as.data.frame()%>%
  column_to_rownames("MAGs")

hc <- agnes(dat_mat, method = "ward")
saveRDS(hc, "data/pMAGs_cluster_select.rds")

sub_grp_6 <- cutree(as.hclust(hc), k = 6)
sub_grp_6 <- as.data.frame(sub_grp_6)
sub_grp_6 <- merge(sub_grp_6, tax, by.x = 0, by.y = "MAGs")
write.table(sub_grp_6, "data/pMAGS_clusters_select.tsv", sep = "\t", quote = F)


