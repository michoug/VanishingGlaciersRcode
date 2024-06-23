library(vegan)
library(janitor)
library(tidyverse)
# library(scales)
library(ggordiplots)

set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

map_clean <- clean_names(map)%>%
  filter(!(sample == "GL140_5"))%>%
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

dat.ano <- with(map_clean, anosim(dat_dist, region, distance = "bray"))
summary(dat.ano)

nmds <- metaMDS(dat_dist, trymax = 1000)

datfort <-  as.data.frame(scores(nmds))

datfort_sites <- datfort %>%
  rownames_to_column(var = "label")

# datfinal <- datfort_sites %>%
#   left_join(map, join_by(label == sample))

nmds$stress <- ifelse(nmds$stress < 1e-2,
                      scientific(nmds$stress, digits = 2),
                      round(nmds$stress, digits = 2))

stress_value <- paste("Stress = ", nmds$stress, sep = "")

yvalue = max(datfort_sites$NMDS2) - 0.05 *
  (max(datfort_sites$NMDS2) - min(datfort_sites$NMDS2))
xvalue = min(datfort_sites$NMDS1) + 0.1 *
  (max(datfort_sites$NMDS1) - min(datfort_sites$NMDS1))


plt <- gg_ordiplot(nmds, groups = map_clean$region, ellipse = F, spiders = T)

p1 <- plt$plot +
  geom_point(size = 3) + 
  annotate("text",
           x = xvalue,
           y = yvalue,
           label = stress_value)+
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
    legend.key = element_blank(),
    legend.margin = margin(b = -5),
  )

p1

ggsave_fitmax("Figures/Fig_Xa_NMDS_pMAGs_Abund.pdf",
              p1,
              maxwidth = 15)
