## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2024-01-10
##
## Copyright (c) Grégoire Michoud, 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(janitor)
library(RColorBrewer)

source("customFunctions/plot_functions.R")

samps  <- read.table(file = "data/metadata_NOMIS_sed_chem.txt", header = TRUE, sep = "\t")
samps <- clean_names(samps)
samps_select <- samps %>%
  select(code_gl, site_c, lat_sn_dd, lon_sn_dd)%>%
  distinct()%>%
  mutate(region = case_when(
    site_c == "New_Zealand" ~ "Southern Alps, New Zealand",
    site_c == "Nepal" ~ "Himalayas, Nepal",
    site_c == "Caucasus" ~ "Caucasus Mountains, Russia",
    site_c == "Kyrgyzstan" ~ "Pamir and Tian Shan, Kyrgyzstan",
    site_c == "Alps" ~ "European Alps",
    site_c == "Norway" ~ "Scandinavian Mountains, Norway",
    site_c == "Greenland" ~ "Southwest Greenland",
    site_c == "Uganda" ~ "Rwenzori Mountains, Uganda",
    site_c == "Ecuador" ~ "Ecuadorian Andes"
  ))


world <- ne_countries(scale = "medium", returnclass = "sf")

p1 <- ggplot(data = world) +
  geom_sf(fill = NA)+
  coord_sf(xlim = c(-85, 180), ylim = c(-50, 73), expand = FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  labs(color = "Region")+
  # geom_text_repel(data= samps_select,aes(x=lon_sn_dd, y=lat_sn_dd, label=code_gl,
  #           color = site_c), fontface = "bold", max.overlaps = 100)+
  geom_point(data= samps_select,aes(x=lon_sn_dd, y=lat_sn_dd,
                                         color = region), size = 3.5)+
  scale_color_brewer(palette = "Paired")+
  theme_light()
p1
ggsave_fitmax("Figures/Fig_S1_map.pdf",p1, maxwidth = 12)
