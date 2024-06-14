library(vegan)
library(ape)
library(janitor)
# library(ggvegan)
library(tidyverse)
library(RColorBrewer)

set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

map_clean <- clean_names(map)


dat_m <- dat %>%
  mutate_if(is.double, function(x, na.rm = FALSE)
    (x * 100)) %>%
  mutate_if(is.double, as.integer) %>%
  select(MAGs, map_clean$sample)%>%
  as.data.frame()%>%
  column_to_rownames(var = "MAGs")


dat_dist <- t(dat_m) %>%
  vegdist(method = "bray")

all(labels(dat_dist) == map_clean$sample)

dat.ano <- with(map_clean, anosim(dat_dist, site_c, distance = "bray"))
summary(dat.ano)

pcoa <- pcoa(dat_dist)

datfort <-  fortify(pcoa$vectors) %>%
  select(Axis.1, Axis.2) %>%
  rownames_to_column(var = "label") %>%
  rename(PCOA1 = Axis.1) %>%
  rename(PCOA2 = Axis.2)

datfort_sites <- datfort %>%
  filter(!(label == "GL140_5"))

datfinal <- datfort_sites %>%
  left_join(map_clean, join_by(label == sample))%>%
  select(PCOA1, PCOA2, label, site_c)%>%
  mutate(region = case_when(
    site_c == "New_Zealand" ~ "Southern Alps, New Zealand",
    site_c == "Nepal" ~ "Himalayas, Nepal",
    site_c == "Caucasus" ~ "Caucasus Mountains, Russia",
    site_c == "Kyrgyzstan" ~ "Pamir and Tian Shan, Kyrgyzstan",
    site_c == "Alps" ~ "European Alps",
    site_c == "Norway" ~ "Scandinavian Mountains, Norway",
    site_c == "Greenland" ~ "Southwest Greenland",
    site_c == "Uganda" ~ "Rwenzori Mountains, Uganda",
    site_c == "Ecuador" ~ "Ecuadorian Andes",
    .default = site_c
  ))

p <- ggplot(data = datfinal, aes(
  x = PCOA1,
  y = PCOA2,
  label = label,
  color = region
)) +
  geom_point(size = 3) + # add the point markers
  coord_equal() +
  theme(
    axis.title = element_text(
      size = 14,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 10, colour = "grey30"),
    legend.position = 'top',
    legend.box = 'vertical',
    legend.key = element_blank()
  )+
  scale_color_brewer(palette = "Paired")+
  guides(color=guide_legend(nrow=3,byrow=TRUE))+
  labs(color = "Region")
p


ggsave_fitmax("Figures/Fig_4b_PCOA_all.pdf",
              p,
              maxwidth = 20)
