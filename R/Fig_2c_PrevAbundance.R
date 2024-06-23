library(tidyverse)

mags_tab <- read_tsv("data/pMAGs_cov_norm.txt.gz")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
specgen <- read_tsv("data/pMAGs_spec_gen.tsv")
magsRocks <- read_tsv("data/pMAGsInRocks.txt")
meta <- read_tsv("data/metadata_NOMIS_sed_chem.txt")


mags_tab[mags_tab < 100] = 0

prevAbund <- mags_tab %>%
  pivot_longer(cols = !MAGs)%>%
  filter(name %in% meta$Sample)%>%
  filter(!(MAGs %in% magsRocks$MAGs))%>%
  group_by(MAGs)%>%
  mutate(abundance = sum(value),
         n_occu =n(),
         n_gt0= sum(value > 0),
         prevalence = n_gt0 / n_occu)%>%
  ungroup()%>%
  select(-c(name, value))%>%
  distinct()%>%
  mutate(abundance_perc = abundance / sum(abundance))%>%
  left_join(specgen, join_by(MAGs))

p1 <- ggplot(prevAbund, aes(y = log10(abundance_perc*100), x = prevalence*100, color = sign))+
  geom_point()+
  labs(x = "Prevalence (%)", y = "Abundance (log10 %)", color = "Niche Breadth")+
  # scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3"))+
  scale_color_manual(values = c("#66c2a5","grey","#fc8d62"))+
  theme_classic()+
  theme(
    axis.title = element_text(
      size = 14,
      colour = "black"
    ),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    legend.text = element_text(size = 10, colour = "black"),
    legend.position = 'top',
    legend.box = 'vertical',
    legend.key = element_blank(),
    legend.margin = margin(b = -5),
  )
p1

ggsave("Figures/Fig_Xc_Niche_prev_abun.pdf", p1)
