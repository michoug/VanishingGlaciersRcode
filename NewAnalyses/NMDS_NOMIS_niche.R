library(vegan)
library(ggrepel)
library(janitor)
library(ggvegan)
library(tidyverse)

set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
spec_gen <- read_tsv("NewAnalyses/pMAGS_spec_gen.tsv")

type = "GENERALIST"

map_clean <- clean_names(map)

spec_gen <- spec_gen %>%
  select(MAGs, sign)

dat_m <- dat %>%
  left_join(spec_gen, join_by(MAGs)) %>%
  filter(sign == !!type) %>%
  mutate_if(is.double, function(x, na.rm = FALSE) (x * 100)) %>%
  mutate_if(is.double, as.integer) %>%
  select(MAGs, map_clean$sample)

dat_m <- as.data.frame(dat_m)
rownames(dat_m) <- dat_m$MAGs
dat_m$MAGs <- NULL

dat_dist <- t(dat_m) %>%
  avgdist(sample = 1000)

dat_rda <- metaMDS(t(dat_m), trymax = 1000)
datfort <-  fortify(dat_rda)

datfort_sites <- datfort %>%
  filter(score == "sites") %>%
  filter(!(label == "GL140_5"))


datfinal <- datfort_sites %>%
  left_join(map_clean, join_by(label == sample))

p <- ggplot(data = datfinal, aes(
  x = NMDS1,
  y = NMDS2,
  label = label,
  color = site_c
)) +
  geom_point(size = 3) + # add the point markers
  # geom_text_repel(size=4,vjust=0) +  # add the site labels
  coord_equal() +
  # stat_ellipse()+
  theme_bw()
p
ggsave_fitmax("NewAnalyses/NMDS_MAGs_Abund_generalist.pdf",p, maxwidth = 15)

map_num <- map_clean %>%
  select(where(is.numeric)) %>%
  select(
    "water_temp_c",
    "ph_p_h",
    "conductivity_u_s_cm_1",
    "turb_ntu",
    "ele_sn_m",
    "sn_sp_dist_m",
    "gl_sa_km2",
    "gl_cov_percent",
    "chla_ug_g_1",
    "n3_srp_ug_l_1",
    "n4_nh4_ug_l_1",
    "n5_no3_ug_l_1",
    "n6_no2_ug_l_1"
  )

nmds_envFit = envfit(dat_rda, map_num, na.rm = TRUE, permutations = 10000)

nmds_envScores = as.data.frame(scores(nmds_envFit, display = 'vectors')) %>%
  mutate(env_variable = rownames(.), pvalue = nmds_envFit$vectors$pvals) |>
  filter(pvalue <= 0.05)

en_coord_cont = as.data.frame(scores(nmds_envFit, "vectors")) *
  ordiArrowMul(nmds_envFit)

en_coord_cont_filter <-  en_coord_cont %>%
  filter(row.names(en_coord_cont) %in% nmds_envScores$env_variable)

en_coord_cont_filter


p_nmds = datfinal |>
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site_c), size = 4) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    data = en_coord_cont_filter,
    linewidth = 0.65,
    alpha = 0.5,
    colour = "grey60",
    arrow = arrow(length = unit(0.25, 'cm'))
  ) +
  geom_text(
    data = en_coord_cont_filter,
    aes(x = NMDS1, y = NMDS2),
    colour = "grey30",
    fontface = "bold",
    label = rownames(en_coord_cont_filter)
  ) +
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
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

p_nmds


