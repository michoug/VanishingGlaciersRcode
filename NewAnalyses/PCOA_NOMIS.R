library(vegan)
library(ape)
library(janitor)
library(ggvegan)
library(tidyverse)

set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

map_clean <- clean_names(map)


dat_m <- dat %>%
  mutate_if(is.double, function(x, na.rm = FALSE)
    (x * 100)) %>%
  mutate_if(is.double, as.integer) %>%
  select(MAGs, map_clean$sample)

dat_m <- as.data.frame(dat_m)
rownames(dat_m) <- dat_m$MAGs
dat_m$MAGs <- NULL

dat_dist <- t(dat_m) %>%
  vegdist(method = "bray")

all(labels(dat_dist) == map_clean$sample)

dat.ano <- with(map_clean, anosim(dat_dist, site_c, distance = "bray"))
summary(dat.ano)

pcoa <- pcoa(dat_dist)

ordi <- ordisurf(pcoa ~ map_clean$site_c, plot = FALSE, scaling = 3)

datfort <-  fortify(pcoa$vectors) %>%
  select(Axis.1, Axis.2) %>%
  rownames_to_column(var = "label") %>%
  rename(PCOA1 = Axis.1) %>%
  rename(PCOA2 = Axis.2)

datfort_sites <- datfort %>%
  filter(!(label == "GL140_5"))

datfinal <- datfort_sites %>%
  left_join(map_clean, join_by(label == sample))

p <- ggplot(data = datfinal, aes(
  x = PCOA1,
  y = PCOA2,
  label = label,
  color = site_c
)) +
  geom_point(size = 3) + # add the point markers
  # geom_text_repel(size=4,vjust=0) +  # add the site labels
  coord_equal() +
  # stat_ellipse()+
  # theme_bw()+
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
    legend.key = element_blank(),
    legend.margin = margin(b = -5),
  )
p


ggsave_fitmax("NewAnalyses/PCOA_MAGs_Abund.pdf",
              p,
              maxwidth = 15)
