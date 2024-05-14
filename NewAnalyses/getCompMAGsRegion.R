library(tidyverse)

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
spec <- read_tsv("NewAnalyses/pMAGS_spec_gen.tsv")

df_max <- dat %>%
  pivot_longer(cols = !MAGs)%>%
  right_join(meta, join_by(name == Sample))%>%
  select(MAGs, name, value, Site_c)%>%
  group_by(MAGs, Site_c)%>%
  filter(value > 0)%>%
  summarise(median = median(value))%>%
  slice(which.max(median))%>%
  select(-median)%>%
  rename(MaxSite = Site_c)


dat_merge <- dat %>%
  pivot_longer(cols = !MAGs) %>%
  right_join(meta, join_by(name == Sample))%>%
  left_join(spec, by = "MAGs")%>%
  left_join(df_max, by = "MAGs")%>%
  filter(value > 0)%>%
  select(name, MAGs, sign, MaxSite)

for (i in unique(dat_merge$MaxSite)) {
  dat_sel <- dat_merge%>%
    filter(MaxSite == !!i)%>%
    select(MAGs)%>%
    distinct()
  
  name = paste("NewAnalyses/", i, "_MAGs.tsv", sep = "")
  write_tsv(dat_sel, name)
  
  list <- dat_merge%>%
    filter(sign == "SPECIALIST")%>%
    filter(MaxSite == !!i)%>%
    select(MAGs)%>%
    distinct()
  
  name2 <- paste("NewAnalyses/", i, "_MAGs_spec_select.tsv", sep = "")
  write_tsv(list, name2)
}
