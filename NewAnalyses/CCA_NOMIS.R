library(vegan)
library(ggvegan)
library(janitor)
library(tidyverse)


set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")

map_clean <- clean_names(map) %>%
  filter(!(sample == "GL140_5")) %>%
  mutate(site_c = case_when(
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
  )) %>%
  select(
    "sample",
    "site_c",
    "water_temp_c",
    "ph_p_h",
    "do_mg_l_1",
    "turb_ntu",
    "conductivity_u_s_cm_1",
    "gl_cov_percent",
    "chla_ug_g_1",
    "i6_so4_ug_l_1",
    "n3_srp_ug_l_1",
    "n5_no3_ug_l_1",
  ) %>%
  column_to_rownames(var = "sample")

dat_m <- dat %>%
  mutate_if(is.double, function(x, na.rm = FALSE)
    (x * 100)) %>%
  mutate_if(is.double, as.integer) %>%
  select(MAGs, map$Sample) %>%
  column_to_rownames(var = "MAGs") %>%
  select(-GL140_5)

dat_m <- t(dat_m)
dat_m <- as.data.frame(dat_m)

cca <- capscale(
  dat_m ~ water_temp_c +
    site_c +
    ph_p_h +
    turb_ntu + 
    gl_cov_percent +
    chla_ug_g_1 + 
    # i6_so4_ug_l_1 + 
    n3_srp_ug_l_1,
    # n5_no3_ug_l_1,
  map_clean,
  dist = "bray",
  na.action = na.omit
)

datfort <-  fortify(cca)

datfort_sites <- datfort %>%
  filter(!(label == "GL140_5")) %>%
  filter(score == "sites")

map_final <- map_clean %>%
  rownames_to_column(var = "label")

datfinal <- datfort_sites %>%
  left_join(map_final, join_by(label))

dat_env <- datfort %>%
  filter(score  == "biplot") %>%
  filter(!str_detect(label, "site_c"))%>%
  mutate(label = case_when(
    label == "chla_ug_g_1" ~ "Chla μg.g-1",
    label == "conductivity_u_s_cm_1" ~ "Conductivity (μs.cm-1)",
    label == "do_mg_l_1" ~ "DO (mg l-1)",
    label == "gl_cov_percent" ~ "Galcier Coverage (%)",
    label == "i6_so4_ug_l_1" ~ "SO4 (μg.l-1)",
    label == "n3_srp_ug_l_1" ~ "SRP (μg.l-1)",
    label == "n5_no3_ug_l_1" ~ "NO3 (μg.l-1)",
    label == "ph_p_h" ~ "pH",
    label == "turb_ntu" ~ "Turbity (NTU)",
    label == "water_temp_c" ~ "Water temperature °C",
    .default = label
  ))
  

axes <- paste0('CAP', c(1,2))
exp_var <- summary(cca)$concont$importance[2, axes]
axes <- paste0(axes, ' (', round(100 * exp_var, 2), '%)')

p <- ggplot(data = datfinal, aes(
  x = CAP1,
  y = CAP2,
  label = label,
  color = site_c
)) +
  geom_point(size = 3) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = CAP1,
      yend = CAP2
    ),
    data = dat_env,
    linewidth = 0.65,
    alpha = 0.5,
    colour = "grey60",
    arrow = arrow(length = unit(0.5, 'cm'))
  ) +
  geom_text(
    data = dat_env,
    aes(x = CAP1, y = CAP2),
    colour = "grey30",
    fontface = "bold",
    label = dat_env$label
  ) +
  coord_equal() +
  labs(x = axes[1], y = axes[2], color = "Mountain Range")+
  guides(color=guide_legend(nrow=4, byrow=F))+
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

adonis2(dat_m ~ water_temp_c +
          site_c +
          ph_p_h +
          turb_ntu + 
          gl_cov_percent +
          chla_ug_g_1 + 
          i6_so4_ug_l_1 +
          n3_srp_ug_l_1 +
          n5_no3_ug_l_1,
        map_clean,na.action = na.omit)
