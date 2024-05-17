library(tidyverse)
dat_troph <- read_tsv("data/pMAGs_trophicState.txt")
specgen <- read_tsv("data/pMAGs_spec_gen.tsv")
maxcov <- read_tsv("data/pMAGs_cov_sum.txt")
metadata <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

dat_merge <- merge(dat_troph, specgen, by = "MAGs")

dat_summary <- dat_troph%>%
  right_join(specgen, join_by("MAGs")) %>%
  group_by(sign, TrophicState)%>%
  summarise(n = n())%>%
  ungroup()%>%
  group_by(sign)%>%
  mutate(Frequence = round(n / sum(n), 3))
  

ggplot(dat_summary, aes(x = TrophicState, y = Frequence, fill = sign))+
  geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle =45, hjust = 1))

ggsave("TrophicStateSpecGen.pdf")
