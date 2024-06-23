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

dat <- read_tsv("data/pMAGs_response_vs_explanatory.txt")

dat_4plot <- dat %>%
  clean_names()%>%
  rename(site = "name") %>%
  select(site,contains("trophy"), "none", contains("gl_"), contains("dist"),chla_ug_g_1) %>%
  select(-c(gl_a_km2, phototrophy))%>%
  rename(
    "Chemoorganotrophy\nAerobic respiration\n(%)" = "heterotrophy_aerobic_respiration",
    "Mixotrophy\nHeterotroph/autotroph\n(%)" = "mixotrophy_heterotrophy_autotrophy",
    "Chla (ln ug.g-1)\n" = "chla_ug_g_1"
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


trophicParam <- c("Mixotrophy\nHeterotroph/autotroph\n(%)","Chemoorganotrophy\nAerobic respiration\n(%)")
pChla <- createplot(dat_4plot, "Chla (ln ug.g-1)\n", trophicParam)


p <- ggarrange(plotlist = pChla, labels = "auto", ncol = 2)
p
ggsave_fitmax("Figures/Fig_6ab_trophicvsParameters.pdf",p, maxwidth = 8)
