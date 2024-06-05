library(tidyverse)
# library(janitor)


dat <- read_tsv("NewAnalyses/pMAGs_QuorumSensing_sub.txt")
db <- read_tsv("data/Other/QSAP_categories.txt")
tax <- read_tsv("data/pMAGs_tax.tsv")

tax_num <- tax %>%
  group_by(c)%>%
  mutate(numb = n())

dat_clean_numb_type_I <- dat %>% 
  filter(!(GeneName == "GeneName")) %>%
  select(GeneName, Subtype) %>%
  na.omit() %>%
  # left_join(db, join_by("HmmscanResults" == "protein")) %>%
  mutate(MAGs = gsub("(.*)_(\\d+)_(\\d+)$", "\\1", GeneName, perl = T ))%>%
  group_by(MAGs, Subtype) %>%
  summarise(n = n())%>%
  filter(n > 0)%>%
  mutate(Subtype = gsub("AI-2", "AI_2", Subtype))%>%
  separate(Subtype, into = c("type_I", "type_II", "type_III"), sep = "-")%>%
  group_by(type_I, MAGs)%>%
  summarise(sum = sum(n))%>%
  ungroup()%>%
  left_join(tax_num, join_by(MAGs))%>%
  group_by(type_I, c)%>%
  reframe(n = n_distinct(MAGs)/numb)%>%
  distinct()

dat_clean_type_II <-  dat %>% 
  filter(!(GeneName == "GeneName")) %>%
  select(GeneName, Subtype) %>%
  na.omit() %>%
  # left_join(db, join_by("HmmscanResults" == "protein")) %>%
  mutate(MAGs = gsub("(.*)_(\\d+)_(\\d+)$", "\\1", GeneName, perl = T ))%>%
  group_by(MAGs, Subtype) %>%
  summarise(n = n())%>%
  filter(n > 0)%>%
  mutate(Subtype = gsub("AI-2", "AI_2", Subtype))%>%
  separate(Subtype, into = c("type_I", "type_II", "type_III"), sep = "-")%>%
  group_by(type_I, type_II, MAGs)%>%
  summarise(sum = sum(n))%>%
  ungroup()%>%
  left_join(tax_num, join_by(MAGs))%>%
  group_by(type_I, type_II)%>%
  reframe(n = n_distinct(MAGs))%>%
  distinct()
  
  
length(unique(dat_clean$type_I))

p1 <- ggplot(dat_clean, aes(type_I, sum))+
  geom_violin()+
  geom_point()+
  facet_wrap(~sign)
p1
