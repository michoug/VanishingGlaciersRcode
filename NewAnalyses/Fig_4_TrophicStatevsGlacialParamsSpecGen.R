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
library(fs)
library(janitor)
library(scales)
library(patchwork)

source("customFunctions/plot_functions.R")

params <- c(
  "autotrophy",
  "calvin_cycle",
  "chemotaxis",
  "co_oxidation_c1_metabolism",
  "cobalamin_biosynthesis",
  "dissimilatory_nitrite_reduction",
  "formaldehyde_assimilation",
  "heterotrophy",
  "heterotrophy_aerobic_respiration",
  "mixotrophy_autotrophy_phototrophy",
  "none",
  "p_mag_rich",
  "quorum_sensing",
  "reductive_citrate_cycle",
  "thiosulfate_oxidation_sox",
  "urea_metabolism",
  "vitamin_b12_transport",
  "x3_hp_4_hb_pathway",
  "heterotrophy_anaerobic_respiration",
  "mixotrophy_heterotrophy_autotrophy",
  "mixotrophy_heterotrophy_phototrophy",
  "mixotrophy_other"
)

datSpec <- read_tsv("NewAnalyses/datResponseExplAMGsNOMISMAGsSpec.txt")
datGen <- read_tsv("NewAnalyses/datResponseExplAMGsNOMISMAGsGen.txt")

directory = "temp"

if (!dir.exists(directory)) {
  dir.create(directory)
}

unzip("NewAnalyses/ModelSpec.zip",
      exdir = directory,
      overwrite = T)
unzip("NewAnalyses/ModelGen.zip",
      exdir = directory,
      overwrite = T)

modSpec <- dir_ls(paste(directory, "ModelSpec/", sep = "/")) %>%
  map_dfr(read_tsv, progress = F, show_col_types = F) %>%
  mutate(type = "specialist")


modGen <- dir_ls(paste(directory, "ModelGen/", sep = "/")) %>%
  map_dfr(read_tsv, progress = F, show_col_types = F) %>%
  mutate(type = "generalist")

models <- rbind(modSpec, modGen) %>%
  filter(variables != "(Intercept)") %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  rename(p_value = `Pr(>|t|)`)

models_clean <- models %>%
  mutate(
    variables = case_when(
      variables == "chla_ug_g_1" ~ "Chla (ln ug.g-1)\n",
      variables == "gl_cov_percent" ~ "Glacier Coverage (ln proportion)",
      .default =  variables
    )
  )

dat_4plotSpec <- datSpec %>%
  clean_names() %>%
  rename(site = "name") %>%
  select(site, gl_cov_percent, chla_ug_g_1, all_of(params)) %>%
  rename("Chla (ln ug.g-1)\n" = "chla_ug_g_1",
         "Glacier Coverage (ln proportion)" = "gl_cov_percent") %>%
  mutate(type = "specialist")

dat_4plotGen <- datGen %>%
  clean_names() %>%
  rename(site = "name") %>%
  select(site, gl_cov_percent, chla_ug_g_1, all_of(params)) %>%
  rename(# "Chemoorganotrophy\nAerobic respiration\n(%)" = "heterotrophy_aerobic_respiration",
    # "Mixotrophy\nautotroph/phototroph\n(%)" = "mixotrophy_autotrophy_phototrophy",
    # "Calvin Cycle" = "calvin_cycle",
    # "3HP/4HB Pathway" = "x3_hp_4_hb_pathway",
    # "Chemotaxis" = "chemotaxis",
    "Chla (ln ug.g-1)\n" = "chla_ug_g_1",
    "Glacier Coverage (ln proportion)" = "gl_cov_percent") %>%
  mutate(type = "generalist")

dat_all <- rbind(dat_4plotGen, dat_4plotSpec)

glaParams <- sort(
  c(
    "heterotrophy",
    "chemotaxis",
    "thiosulfate_oxidation_sox",
    "mixotrophy_autotrophy_phototrophy",
    "quorum_sensing"
  )
)

model_sel <- models_clean %>%
  filter(resp %in% glaParams) %>%
  filter(variables == "Chla (ln ug.g-1)\n")

data_sel <- dat_all %>%
  select(all_of(glaParams), "Chla (ln ug.g-1)\n", type) %>%
  pivot_longer(cols = chemotaxis:thiosulfate_oxidation_sox)%>%
  left_join(model_sel, join_by(type, name == resp))%>%
  na.omit()

xmin_all <- min(data_sel[[glacialVar]])
xmax_all <- max(data_sel[[glacialVar]])
ymin_all <- min(data_sel[[var]])
ymax_all <- max(data_sel[[var]])

yvalue = (ymax_all - 0.05 * (ymax_all - ymin_all)) * 100
xvalue = xmin_all + 0.1 *
  (xmax_all - xmin_all)

dat_text <- data_sel %>%
  select(type, p_value) %>%
  distinct() %>%
  mutate(p_value = case_when(
    is.na(p_value) ~  ">0.05",
    p_value < 1e-2 ~ scientific(p_value, digits = 2),
    .default = as.character(round(p_value, digits = 3))
  )) %>%
  mutate(p_value = paste("pvalue =" , p_value, sep = " "))


p1 <- ggplot(data_sel, aes(x = `Chla (ln ug.g-1)\n`, y = value * 100)) +
  geom_point() +
  facet_grid(name ~ type, scales = "free")+
  # ylab(var) +
  geom_smooth(method = "gam", aes(color = p_value)) +
  theme_classic() +
  theme(
    strip.text.y = element_text(angle = 0),
    text = element_text(colour = "black", size = 14),
    panel.grid.major = element_blank()
  )
  facet_grid( ~ type) +
  geom_text(data = dat_text,
            x = xvalue,
            y = yvalue,
            aes(label = p_value))
# return(p1)
plot_list[[var]] <- p1


ggsave_fitmax(
  "NewAnalyses/Fig_4_trophicvsParameters_pChla.pdf",
  p,
  maxheight = 45,
  maxwidth = 10
)
