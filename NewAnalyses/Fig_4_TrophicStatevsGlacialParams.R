## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2024-01-22
##
## Copyright (c) Grégoire Michoud, 2024
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(tidyverse)
library(ggpubr)
library(janitor)

source("customFunctions/plot_functions.R")

dat <- read_tsv("../Correlations/datResponseExplNOMISMAGs.txt")

dat_4plot <- dat %>%
  clean_names()%>%
  rename(site = "name") %>%
  select(site,contains("trophy"), "none", contains("gl_"), contains("dist"),chla_ug_g_1) %>%
  select(-c(gl_a_km2, phototrophy))%>%
  rename(
    "Chemolithoautotrophy (%)\n" = "autotrophy",
    "Chemolitoheterotrophy\nOther" = "heterotrophy",
    "Chemoorganotrophy\nAerobic respiration\n(%)" = "heterotrophy_aerobic_respiration",
    "Chemoorganotrophy\nAnaerobic respiration\n(%)" = "heterotrophy_anaerobic_respiration",
    "Chemoorganotrophy\nFermentation\n(%)" = "heterotrophy_fermentation",
    "Mixotrophy\nautotroph/phototroph\n(%)" = "mixotrophy_autotrophy_phototrophy",
    "Mixotrophy\nHeterotroph/autotroph\n(%)" = "mixotrophy_heterotrophy_autotrophy",
    "Mixotrophy\nheterotroph/phototroph\n(%)" = "mixotrophy_heterotrophy_phototrophy",
    "Mixotrophy\nOther (%)" = "mixotrophy_other",
    "No Trophic State %\n" = "none",
    "Chla (ln ug.g-1)\n" = "chla_ug_g_1",
    "Glacier Coverage (ln proportion)" = "gl_cov_percent",
    "Glacier Surface\nArea (km^2)" = "gl_sa_km2",
    "Distance to\nthe snout (m)" = "sn_sp_dist_m"
  )

createplot <- function(data,glacialVar, TrophicVar){
  
  plot_list <- list()
  glacialVar <- sym(glacialVar)
  
  for (var in TrophicVar) {
    
    p1 <- ggplot(data, aes(x = !!glacialVar, y = .data[[var]]*100))+
      geom_point()+
      ylab(var)+
      geom_smooth(method = "gam")+
      theme_classic()+
      theme(
        strip.text.y = element_text(angle = 0),
        text = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank()
      )
    # return(p1)
    plot_list[[var]] <- p1
    
  }
  
  return(plot_list)
}

trophicParam <- c("Chemoorganotrophy\nAerobic respiration\n(%)", "Mixotrophy\nHeterotroph/autotroph\n(%)")
pCov <- createplot(dat_4plot, "Glacier Coverage (ln proportion)", trophicParam)

trophicParam <- c("Chemoorganotrophy\nAerobic respiration\n(%)")
pChla <- createplot(dat_4plot, "Chla (ln ug.g-1)\n", trophicParam)

plotList <- c(pCov, pChla)
p <- ggarrange(plotlist = plotList, labels = "auto", ncol = 1)
p
ggsave_fitmax("Figures/Fig_4_trophicvsParameters.pdf",p, maxwidth = 4)
