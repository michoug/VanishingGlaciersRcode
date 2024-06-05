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
  summary(dat.ano)
}

getPCOA <- function(dat, map, type) {
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
  
  dat_rda <- metaMDS(t(dat_m), trymax = 1000)
}

plotSimpleNMDS <- function(nmds, map) {
  datfort <-  fortify(nmds)
  
  datfort_sites <- datfort %>%
    filter(score == "sites") %>%
    filter(!(label == "GL140_5"))
  
  datfinal <- datfort_sites %>%
    left_join(map, join_by(label == sample))
  
  nmds$stress <- ifelse(nmds$stress < 1e-2,
                         scientific(nmds$stress, digits = 2),
                         round(nmds$stress, digits=2))
  
  stress_value <- paste("Stress = ", nmds$stress, sep = "")
  yvalue = max(datfort_sites$NMDS2) - 0.05 *
    (max(datfort_sites$NMDS2) - min(datfort_sites$NMDS2))
  xvalue = min(datfort_sites$NMDS1) + 0.1 *
    (max(datfort_sites$NMDS1) - min(datfort_sites$NMDS1))
  
  p <- ggplot(data = datfinal, aes(
    x = NMDS1,
    y = NMDS2,
    label = label,
    color = site_c
  )) +
    annotate("text", x = xvalue, y=yvalue, label= stress_value)+
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



plotEnvFitNMDS <- function(nmds, map) {
  datfort <-  fortify(nmds)
  
  datfort_sites <- datfort %>%
    filter(score == "sites") %>%
    filter(!(label == "GL140_5"))
  
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
  
  nmds_envFit = envfit(nmds, map_num, na.rm = TRUE, permutations = 10000)
  
  nmds_envScores = as.data.frame(scores(nmds_envFit, display = 'vectors')) %>%
    mutate(env_variable = rownames(.),
           pvalue = nmds_envFit$vectors$pvals) |>
    filter(pvalue <= 0.05)
  
  en_coord_cont = as.data.frame(scores(nmds_envFit, "vectors")) *
    ordiArrowMul(nmds_envFit)
  
  en_coord_cont_filter <-  en_coord_cont %>%
    filter(row.names(en_coord_cont) %in% nmds_envScores$env_variable)
  

  nmds$stress <- ifelse(nmds$stress < 1e-2,
                        scientific(nmds$stress, digits = 2),
                        round(nmds$stress, digits=3))
  
  stress_value <- paste("Stress = ", nmds$stress, sep = "")
  xmin_all <- min(c(datfort_sites$NMDS1, en_coord_cont_filter$NMDS1))
  xmax_all <- max(c(datfort_sites$NMDS1, en_coord_cont_filter$NMDS1))
  ymin_all <- min(c(datfort_sites$NMDS2, en_coord_cont_filter$NMDS2))
  ymax_all <- max(c(datfort_sites$NMDS2, en_coord_cont_filter$NMDS2))
  
  yvalue = ymax_all - 0.05 *
    (ymax_all - ymin_all)
  xvalue = xmin_all + 0.1 *
    (xmax_all - xmin_all)
  
  
  p_nmds = datfinal |>
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = site_c), size = 4) +
    annotate("text", x = xvalue, y=yvalue, label= stress_value)+
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
    )
  p_nmds
}


anos_gen <- getDist(dat, map_clean, "GENERALIST")
anos_non <- getDist(dat, map_clean, "NON SIGNIFICANT")
anos_spec <- getDist(dat, map_clean, "SPECIALIST")

nmds_spec <- getNMDS(dat, map_clean, "SPECIALIST")
nmds_gen <- getNMDS(dat, map_clean, "GENERALIST")
nmds_non <- getNMDS(dat, map_clean, "NON SIGNIFICANT")

plotSimpleNMDS(nmds_non, map_clean)

p1 <- plotSimpleNMDS(nmds_gen, map_clean)
p1


p2 <- plotEnvFitNMDS(nmds_spec, map_clean)
p2
p_nmds <- ggarrange(
  p1,
  p2,
  common.legend = TRUE,
  legend = "top",
  nrow = 1,
  labels = "auto",
  align = "h"
)
p_nmds

ggsave_fitmax("NewAnalyses/NMDS_MAGs_Abund_gen_spec.pdf", p_nmds, maxwidth = 15)
