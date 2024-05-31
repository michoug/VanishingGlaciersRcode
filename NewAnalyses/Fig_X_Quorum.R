library(tidyverse)
# library(janitor)


dat <- read_tsv("NewAnalyses/pMAGs_QuorumSensing_sub.txt")
db <- read_tsv("data/Other/QSAP_categories.txt")


dat_clean <- dat %>% 
  filter(!(GeneName == "GeneName")) %>%
  select(GeneName, Subtype) %>%
  na.omit() %>%
  # left_join(db, join_by("HmmscanResults" == "protein")) %>%
  mutate(MAGs = gsub("(.*)_(\\d+)_(\\d+)$", "\\1", GeneName, perl = T ))%>%
  group_by(MAGs, Subtype) %>%
  summarise(n = n())
