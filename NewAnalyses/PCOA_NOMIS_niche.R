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
  
  pcoa <- pcoa(dat_dist)
}

plotSimplePCOA <- function(pcoa, map) {
  datfort <-  fortify(pcoa$vectors)%>%
    select(Axis.1, Axis.2)%>%
    rownames_to_column(var = "label")%>%
    rename(PCOA1 = Axis.1)%>%
    rename(PCOA2 = Axis.2)
  
  datfort_sites <- datfort %>%
    filter(!(label == "GL140_5"))
  
  datfinal <- datfort_sites %>%
    left_join(map, join_by(label == sample))
  
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
}

anos_gen <- getDist(dat, map_clean, "GENERALIST")
anos_non <- getDist(dat, map_clean, "NON SIGNIFICANT")
anos_spec <- getDist(dat, map_clean, "SPECIALIST")

pcoa_spec <- getPCOA(dat, map_clean, "SPECIALIST")
pcoa_gen <- getPCOA(dat, map_clean, "GENERALIST")
pcoa_non <- getPCOA(dat, map_clean, "NON SIGNIFICANT")

p1 <- plotSimplePCOA(pcoa_gen, map_clean)
p1

p2 <- plotSimplePCOA(pcoa_spec, map_clean)
p2

p_pcoa <- ggarrange(
  p1,
  p2,
  common.legend = TRUE,
  legend = "top",
  nrow = 1,
  labels = "auto",
  align = "h"
)
p_pcoa

ggsave_fitmax("NewAnalyses/PCOA_MAGs_Abund_gen_spec.pdf", p_pcoa, maxwidth = 15)
