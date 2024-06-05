library(vegan)
library(ggrepel)
library(janitor)
library(ggvegan)
library(tidyverse)
library(scales)
library(ggpubr)

set.seed(123)
source("customFunctions/plot_functions.R")

dat <- read_tsv("data/pMAGs_cov_norm.txt.gz")
map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
spec_gen <- read_tsv("data/pMAGS_spec_gen.tsv")

map_clean <- clean_names(map)
map <- map_clean

getDist <-  function(dat, map, type) {
  spec_gen <- spec_gen %>%
    select(MAGs, sign)
  
  dat_m <- dat %>%
    left_join(spec_gen, join_by(MAGs)) %>%
    filter(sign == !!type) %>%
    mutate_if(is.double, function(x, na.rm = FALSE)
      (x * 100)) %>%
    mutate_if(is.double, as.integer) %>%
    select(MAGs, map_clean$sample)
  
  dat_m <- as.data.frame(dat_m)
  rownames(dat_m) <- dat_m$MAGs
  dat_m$MAGs <- NULL
  
  dat_dist <- t(dat_m) %>%
    vegdist(method = "bray")
  
  dat.ano <- with(map, anosim(dat_dist, site_c, distance = "bray"))
  sum_ano <- summary(dat.ano)
  return(sum_ano)
}

getCCA <- function(dat, map, type) {
  spec_gen <- spec_gen %>%
    select(MAGs, sign)
  
  dat_m <- dat %>%
    left_join(spec_gen, join_by(MAGs)) %>%
    filter(sign == !!type) %>%
    mutate_if(is.double, function(x, na.rm = FALSE)
      (x * 100)) %>%
    mutate_if(is.double, as.integer) %>%
    select(MAGs, map_clean$sample)
  
  dat_m <- as.data.frame(dat_m)
  rownames(dat_m) <- dat_m$MAGs
  dat_m$MAGs <- NULL
  
  cca <- cca(t(dat_m))
}

plotSimpleCCA <- function(cca, map) {
  datfort <-  fortify(cca)
  
  datfort_sites <- datfort %>%
    filter(!(label == "GL140_5"))%>%
    filter(score == "sites")
  
  datfinal <- datfort_sites %>%
    left_join(map, join_by(label == sample))
  
  p <- ggplot(data = datfinal, aes(
    x = CA1,
    y = CA2,
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
}

plotEnvFitCCA <- function(cca, map) {
  datfort <-  fortify(cca)
  
  datfort_sites <- datfort %>%
    filter(!(label == "GL140_5"))%>%
    filter(score == "sites")
  
  datfinal <- datfort_sites %>%
    left_join(map, join_by(label == sample))
  
  map_num <- map %>%
    select(where(is.numeric)) %>%
    select(
      "gl_cov_percent",
      "sn_sp_dist_m",
      "chla_ug_g_1",
      "conductivity_u_s_cm_1",
      "ph_p_h",
      "turb_ntu",
      "water_temp_c"
    )
  
  cca_envFit = envfit(cca, map_num, na.rm = TRUE, permutations = 10000)
  
  cca_envScores = as.data.frame(scores(cca_envFit, display = 'vectors')) %>%
    mutate(env_variable = rownames(.),
           pvalue = cca_envFit$vectors$pvals) |>
    filter(pvalue <= 0.05)
  
  en_coord_cont = as.data.frame(scores(cca_envFit, "vectors")) *
    ordiArrowMul(cca_envFit)
  
  en_coord_cont_filter <-  en_coord_cont %>%
    filter(row.names(en_coord_cont) %in% cca_envScores$env_variable)
  
  
  p_cca = datfinal |>
    ggplot(aes(x = CA1, y = CA2)) +
    geom_point(aes(color = site_c), size = 4) +
    geom_segment(
      aes(
        x = 0,
        y = 0,
        xend = CA1,
        yend = CA2
      ),
      data = en_coord_cont_filter,
      linewidth = 0.65,
      alpha = 0.5,
      colour = "grey60",
      arrow = arrow(length = unit(0.25, 'cm'))
    ) +
    geom_text(
      data = en_coord_cont_filter,
      aes(x = CA1, y = CA2),
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
    )
  p_cca
}

anos_gen <- getDist(dat, map_clean, "GENERALIST")
anos_non <- getDist(dat, map_clean, "NON SIGNIFICANT")
anos_spec <- getDist(dat, map_clean, "SPECIALIST")

cca_spec <- getCCA(dat, map_clean, "SPECIALIST")
cca_gen <- getCCA(dat, map_clean, "GENERALIST")
cca_non <- getCCA(dat, map_clean, "NON SIGNIFICANT")



plotEnvFitCCA(cca_spec, map_clean)

p1 <- plotSimpleCCA(cca_gen, map_clean)
p1

p2 <- plotSimpleCCA(cca_spec, map_clean)
p2

p_cca <- ggarrange(
  p1,
  p2,
  common.legend = TRUE,
  legend = "top",
  nrow = 1,
  labels = "auto",
  align = "h"
)
p_cca

ggsave_fitmax("NewAnalyses/CCA_MAGs_Abund_gen_spec.pdf", p_cca, maxwidth = 15)
