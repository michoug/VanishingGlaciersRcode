library(tidyverse)

source("customFunctions/plot_functions.R")

dat_troph <- read_tsv("data/pMAGs_trophicState.txt")
specgen <- read_tsv("data/pMAGs_spec_gen.tsv")
maxcov <- read_tsv("data/pMAGs_cov_sum.txt")
metadata <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

dat_merge <- merge(dat_troph, specgen, by = "MAGs")

dat_summary <- dat_troph%>%
  right_join(specgen, join_by("MAGs")) %>%
  filter(sign != "NON SIGNIFICANT")%>%
  group_by(sign, TrophicState)%>%
  summarise(n = n())%>%
  ungroup()%>%
  group_by(sign)%>%
  mutate(Frequence = round(n / sum(n), 3))%>%
  ungroup()%>%
  complete(sign, TrophicState, fill = list(n = 0, Frequence = 0))%>%
  mutate(
    TrophicState = case_when(
      TrophicState == "Autotrophy" ~ "Chemolitoautotrophy",
      TrophicState == "Heterotrophy" ~ "Chemolitoheterotrophy\nMix",
      TrophicState == "Heterotrophy_aerobic_respiration" ~ "Chemoorganotrophy\naerobic respiration",
      TrophicState == "Heterotrophy_anaerobic_respiration" ~ "Chemoorganotrophy\nanaerobic respiration",
      TrophicState == "Heterotrophy_fermentation" ~ "Chemoorganotrophy\nfermentation",
      TrophicState == "Mixotrophy_Autotrophy-Phototrophy" ~ "Mixotrophy\nautotroph/phototroph",
      TrophicState == "Mixotrophy_Heterotrophy-Autotrophy" ~ "Mixotrophy\nheterotroph/autotroph",
      TrophicState == "Mixotrophy_Heterotrophy-Phototrophy" ~ "Mixotrophy\nheterotroph/phototroph",
      TrophicState == "Mixotrophy_Other" ~ "Mixotrophy\nother",
      .default = TrophicState
    ))
  
  

p <- ggplot(dat_summary, aes(x = TrophicState, y = Frequence, fill = sign))+
  geom_bar(stat = "identity",position = "dodge")+
  theme_classic()+
  theme(axis.text.x = element_text(angle =45, hjust = 1))+
  labs(fill = "Niche",x = "Trophic State")

ggsave_fitmax("Figures/Fig_S6_trophic_state_niche.pdf", p, maxwidth = 10)
