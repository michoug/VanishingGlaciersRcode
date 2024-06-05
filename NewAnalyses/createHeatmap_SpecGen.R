library(tidyverse)
library(janitor)
library(fs)

source("customFunctions/plot_functions.R")

directory = "temp"

if (!dir.exists(directory)) {dir.create(directory)}

unzip("NewAnalyses/ModelSpec.zip", exdir = directory, overwrite = T)
unzip("NewAnalyses/ModelGen.zip", exdir = directory, overwrite = T)

modSpec <- dir_ls(paste(directory, "ModelSpec/", sep = "/"))%>%
  map_dfr(read_tsv, progress = F, show_col_types = F)%>%
  mutate(type = "specialist")


modGen <- dir_ls(paste(directory, "ModelGen/", sep = "/"))%>%
  map_dfr(read_tsv, progress = F, show_col_types = F)%>%
  mutate(type = "generalist")

models <- rbind(modSpec, modGen)

dat_heat <- models%>%
  filter(!(Estimate == "Estimate"))%>%
  filter(!(variables == "(Intercept)"))%>%  
  filter(`Pr(>|t|)`<= 0.05)%>%
  clean_names()%>%
  select(resp, variables, t_value, type)%>%
  mutate(t_value = as.numeric(t_value))%>%
  complete(resp,variables,type,fill=list(t_value=0))%>%
  filter(!(variables %in% c("pili", "flagella")))


levelsResp <- c("vitamin_b12_transport","thiamine_transporter",
                "thiamine_biosynthesis_prokaryotes", "quorum_sensing",
                "dmsp","cobalamin_biosynthesis","chemotaxis","biotin_tranport",
                "biotin_biosynthesis","co_oxidation_c1_metabolism",
                "thiosulfate_oxidation_sox","sulfite_oxidation","ass_diss_sulfate_reduction",
                "urea_metabolism", "denitrification", "dissimilatory_nitrate_reduction",
                "dissimilatory_nitrite_reduction", "ammonia_oxidation","reductive_citrate_cycle",
                "formaldehyde_assimilation", "calvin_cycle", "x3_hp_4_hb_pathway",
                "phototrophy", "none","mixotrophy_other","mixotrophy_heterotrophy_phototrophy",
                "mixotrophy_heterotrophy_autotrophy","mixotrophy_autotrophy_phototrophy",
                "heterotrophy_fermentation", "heterotrophy_anaerobic_respiration", 
                "heterotrophy_aerobic_respiration", "heterotrophy", "autotrophy",
                "p_mag_rich","p_mag_diversity")

levelsExp <- c("gl_cov_percent","chla_ug_g_1")#,"sn_sp_dist_m",
               # "conductivity_u_s_cm_1","ph_p_h","turb_ntu","water_temp_c")

  
dat_heat$resp <- factor(dat_heat$resp, levels = levelsResp )
dat_heat$variables <- factor(dat_heat$variables, levels = levelsExp)

dat_heat_filter <- dat_heat %>%
  na.omit() %>%
  group_by(resp) %>%
  mutate(sum = sum(t_value))%>%
  ungroup()%>%
  filter(sum != 0)

p1 <- ggplot(dat_heat_filter, aes(variables,resp, fill = t_value))+
  geom_tile()+
  scale_fill_gradient2()+
  facet_wrap(~type)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p1

ggsave_fitmax("NewAnalyses/heatmapGamsSpecGen.pdf", p1, maxwidth = 10)
