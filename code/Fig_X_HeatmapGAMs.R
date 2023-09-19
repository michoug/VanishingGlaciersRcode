## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Grégoire Michoud
##
## Date Created: 2023-09-18
##
## Copyright (c) Grégoire Michoud, 2023
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


library(tidyverse)
library(janitor)
library(ggpubr)

source("customFunctions/plot_functions.R")

dat_mod <- read_tsv("../Correlations/Models_good.txt")


dat_heat <- dat_mod%>%
  filter(!(Estimate == "Estimate"))%>%
  filter(!(variables == "(Intercept)"))%>%  
  clean_names()%>%
  select(resp, variables, t_value)%>%
  mutate(t_value = as.numeric(t_value))%>%
  complete(resp,variables,fill=list(t_value=0))%>%
  filter(!(resp %in% c("pili", "flagella")))%>%
  mutate(resp = case_when(
    resp == "alg_mag_diversity" ~ "Algal Diversity",
    resp == "alg_mag_rich" ~ "Algal Richness",
    resp == "ammonia_oxidation" ~ "Ammonia oxidation",
    resp == "ass_diss_sulfate_reduction" ~ "Assimilatory/dissimilatory\nsulfate reduction",
    resp == "assimilatory_nitrate_reduction" ~ "Assimilatory nitrate reduction",
    resp == "autotrophy" ~ "Chemolithoautotrophy",
    resp == "biotin_biosynthesis" ~ "Biotin biosynthesis",
    resp == "biotin_tranport" ~ "Biotin transport",
    resp == "calvin_cycle" ~ "Calvin cycle",
    resp == "chemotaxis" ~ "Chemotaxis",
    resp == "co_oxidation_c1_metabolism" ~ "CO oxidation", 
    resp == "cobalamin_biosynthesis" ~ "Cobalamin biosynthesis",
    resp == "denitrification" ~ "Denitrification",
    resp == "dissimilatory_nitrate_reduction" ~ "Dissimilatory nitrate reduction",
    resp == "dissimilatory_nitrite_reduction" ~ "Dissimilatory nitrite reduction",
    # resp == "dissimilatory_sulfate_reduction" ~ "Dissimilatory sulfate reduction",
    resp == "dmsp" ~ "DMSP",
    # resp == "flagella" ~ "Flagella",
    resp == "formaldehyde_assimilation" ~ "Formaldehyde assimilation",
    resp == "heterotrophy" ~ "Chemolitoheterotrophy\nMix",
    resp == "heterotrophy_aerobic_respiration" ~ "Chemoorganotrophy\naerobic respiration",
    resp == "heterotrophy_anaerobic_respiration" ~ "Chemoorganotrophy\nanaerobic respiration",
    resp == "heterotrophy_fermentation" ~ "Chemoorganotrophy\nfermentation",
    resp == "mixotrophy_autotrophy_phototrophy" ~ "Mixotrophy\nautotroph/phototroph",
    resp == "mixotrophy_heterotrophy_autotrophy" ~ "Mixotrophy\nheterotroph/autotroph",
    resp == "mixotrophy_heterotrophy_phototrophy" ~ "Mixotrophy\nheterotroph/phototroph",
    resp == "mixotrophy_other" ~ "Mixotrophy\nother",
    resp == "none" ~ "None",
    resp == "p_mag_diversity" ~ "Prokaryotic Diversity",
    resp == "p_mag_rich" ~ "Prokaryotic Richness",
    # resp == "pili" ~ "Pili",
    resp == "quorum_sensing" ~ "Quorum sensing",
    resp == "reductive_citrate_cycle" ~ "Reductive citrate cycle",
    resp == "sulfite_oxidation" ~ "Sulfite oxidation",
    resp == "thiamine_biosynthesis_prokaryotes" ~ "Thiamine biosynthesis",
    resp == "thiamine_transporter" ~ "Thiamine transporter",
    resp == "thiosulfate_oxidation_sox" ~ "Thiosulfate oxidation",
    resp == "urea_metabolism" ~ "Urea metabolism",
    resp == "v_mag_diversity" ~ "Viral Diversity",
    resp == "v_mag_rich" ~ "Viral Richness",
    resp == "vitamin_b12_transport" ~ "Vitamin B12 transport",
    # resp == "wood_ljungdahl_pathway" ~ "Wood-Ljungdahl pathway",
    resp == "x3_hp_4_hb_pathway" ~ "3-HP 4-HB pathway",
    .default = "Other"
  ))%>%
  filter(!(resp == "Other"))%>%
  mutate(variables = case_when(
    variables == "chla_ug_g_1" ~ "Chla",
    variables == "conductivity_u_s_cm_1" ~ "Conductivity",
    variables == "ele_sn_m" ~ "Elevation",
    variables == "gl_cov_percent" ~ "Glacier Coverage",
    variables == "gl_sa_km2" ~ "Glacier Surface\nArea",
    variables == "n3_srp_ug_l_1" ~ "SRP",
    variables == "n4_nh4_ug_l_1" ~ "NH4",
    variables == "n5_no3_ug_l_1" ~ "NO3",
    variables == "n6_no2_ug_l_1" ~ "NO2",
    variables == "ph_p_h" ~ "pH",
    variables == "sn_sp_dist_m" ~ "Distance to\nthe snout", 
    variables == "turb_ntu" ~ "Turbidity",
    variables == "water_temp_c" ~ "Water\ntemperature"
  ))

levelsResp <- c("Vitamin B12 transport", "Thiamine transporter", "Thiamine biosynthesis",
                "Quorum sensing", "DMSP", "Cobalamin biosynthesis", "Chemotaxis",
                "Biotin transport", "Biotin biosynthesis", "CO oxidation",
                "Thiosulfate oxidation", "Sulfite oxidation",
                "Assimilatory/dissimilatory\nsulfate reduction", "Urea metabolism",
                "Denitrification", "Dissimilatory nitrite reduction","Dissimilatory nitrate reduction",
                "Assimilatory nitrate reduction", "Ammonia oxidation",
                "Reductive citrate cycle", "Formaldehyde assimilation", "Calvin cycle",
                "3-HP 4-HB pathway", "None", "Mixotrophy\nother", "Mixotrophy\nheterotroph/phototroph",
                "Mixotrophy\nheterotroph/autotroph", "Mixotrophy\nautotroph/phototroph",
                "Chemoorganotrophy\nfermentation", "Chemoorganotrophy\nanaerobic respiration",
                "Chemoorganotrophy\naerobic respiration", "Chemolitoheterotrophy\nMix",
                "Chemolithoautotrophy", "Viral Richness", "Viral Diversity",
                "Prokaryotic Richness", "Prokaryotic Diversity", "Algal Richness", "Algal Diversity")

algal <- levelsResp[1:9]
sulfur <- levelsResp[11:13]
nitrogen <- levelsResp[14:19]
carbon <- levelsResp[20:23]
trophic <- levelsResp[24:33]
div <- levelsResp[34:length(levelsResp)]

levelsExp <- c("Elevation","Glacier Coverage","Glacier Surface\nArea","Distance to\nthe snout",
               "Chla","Conductivity","SRP","NH4","NO3","NO2","pH",
               "Turbidity","Water\ntemperature")

glacial <- levelsExp[1:4]
  
dat_heat$resp <- factor(dat_heat$resp, levels = levelsResp )
dat_heat$variables <- factor(dat_heat$variables, levels = levelsExp)

dat_heat_final <- dat_heat %>%
  mutate(cat_elem = case_when(
    variables %in% glacial~ "Glacial variables",
    .default = "Physicochemical variables"
  ))%>%
  mutate(cat_resp = case_when(
    resp %in% algal ~ "Algal\nGenes",
    resp == "CO oxidation" ~ "C1 metabolism",
    resp %in% sulfur ~ "Sulfur\nmetabolism",
    resp %in% nitrogen ~ "Nitrogen\nmetabolism",
    resp %in% carbon ~ "Carbon\nmetabolism",
    resp %in% trophic ~ "Trophic\nstate",
    resp %in% div ~ "Diversity\nRichness"
  ))

p1 <- ggplot(dat_heat_final, aes(variables,resp, fill = t_value))+
  geom_tile(color = "black")+
  scale_fill_gradient2()+
  theme_minimal()+
  labs(x = NULL, y = NULL, fill = "T value")+
  facet_grid(rows = vars(cat_resp), 
             cols = vars(cat_elem), scales = "free", space="free_y")+
  theme(strip.text.y = element_text(angle = 0),
        text = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank())

p1

ggsave_fitmax("Figures/Fig_X_heatmapGams.pdf", p1, maxwidth = 15, maxheight = 13)
