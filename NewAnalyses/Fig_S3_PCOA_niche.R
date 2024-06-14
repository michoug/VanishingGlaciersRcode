library(vegan)
library(janitor)
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
  return(dat.ano)
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

plotSimplePCOA <- function(pcoa, map, ano) {
  datfort <-  fortify(pcoa$vectors)%>%
    select(Axis.1, Axis.2)%>%
    rownames_to_column(var = "label")%>%
    rename(PCOA1 = Axis.1)%>%
    rename(PCOA2 = Axis.2)
  
  datfort_sites <- datfort %>%
    filter(!(label == "GL140_5"))
  
  datfinal <- datfort_sites %>%
    left_join(map, join_by(label == sample))%>%
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
  
  xmin_all <- min(datfinal$PCOA1)
  xmax_all <- max(datfinal$PCOA1)
  ymin_all <- min(datfinal$PCOA2)
  ymax_all <- max(datfinal$PCOA2)
  
  yvalue = (ymax_all - 0.05 * (ymax_all - ymin_all))
  xvalue = xmin_all + 0.1 *
    (xmax_all - xmin_all)
  
  text <- paste("R = ", round(ano$statistic, digits = 2), "\npvalue = ", ano$signif, sep = "")
  dat_text <- data.frame(text, xvalue, yvalue, region = "black")
  
  p <- ggplot() +
    geom_point(data = datfinal, aes(
      x = PCOA1,
      y = PCOA2,
      color = region),size = 3) + # add the point markers
    coord_equal() +
    geom_text(data = dat_text,
              x = xvalue,
              y = yvalue,
              aes(label = text))+
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
    )+
    guides(color=guide_legend(nrow=3,byrow=TRUE))+
    labs(color = "Region")
  p
}

anos_gen <- getDist(dat, map_clean, "GENERALIST")
anos_spec <- getDist(dat, map_clean, "SPECIALIST")

pcoa_spec <- getPCOA(dat, map_clean, "SPECIALIST")
pcoa_gen <- getPCOA(dat, map_clean, "GENERALIST")

p1 <- plotSimplePCOA(pcoa_gen, map_clean, anos_gen)
p1

p2 <- plotSimplePCOA(pcoa_spec, map_clean, anos_spec)
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

ggsave_fitmax("Figures/Fig_S4_PCOA_niche.pdf", p_pcoa, maxwidth = 15)
