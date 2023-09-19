## ---------------------------
##
## Script name: plotLengthComplCont.R
##
## Purpose of script: Plot different parameters describing MAGs
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-14
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(ggpubr)
library(tidyverse)

source("customFunctions/plot_functions.R")

dat_mags	<- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGs_checkm_stats.txt")
tax <- read_tsv("../Prokaryotes/TaxQual/NOMIS_MAGS_tax.tsv")
stats <- read_tsv("../Prokaryotes/Stats_rRNA//NOMIS_MAGs_statsTable.txt")

dat_plot <- dat_mags %>%
  left_join(tax)%>%
  left_join(stats)%>%
  mutate(Quality =if_else(Completeness >= 90 & Contamination <= 5,"High","Medium"))

dat_plot$N50 <- dat_plot$N50 / 1000

p1 <- ggboxplot(dat_plot, x = "Quality", y = c("N50","Completeness","Count", "Contamination"),
          color = "Quality", palette = "jco", title = "",
          xlab = "",
          ylab = c("N50 (kbp)","Completeness (%)","Contigs number", "Contamination (%)"))

p1$N50 <- p1$N50 + rotate()
p1$Completeness <- p1$Completeness + rotate()
p1$Contamination <- p1$Contamination + rotate()
p1$Count <- p1$Count + rotate()
  
p2 <- ggarrange(plotlist = p1, common.legend = TRUE)

ggsave_fitmax("Figures/Fig_1a_NOMIS_MAGs_stats.pdf", p2, maxwidth = 10)
