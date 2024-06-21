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
library(janitor)
library(scales)
library(patchwork)

source("customFunctions/plot_functions.R")

models <- read_tsv("data/pMAGs_gam_models_spec_gen.txt")
datSpec <- read_tsv("data/pMAGs_response_vs_explanatory_specialists.txt")
datGen <- read_tsv("data/pMAGs_response_vs_explanatory_generalists.txt")

glaParams <- sort(
  c(
    "heterotrophy",
    "chemotaxis",
    "thiosulfate_oxidation_sox",
    "mixotrophy_autotrophy_phototrophy",
    "quorum_sensing"
  )
)


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
  select(site, gl_cov_percent, chla_ug_g_1, all_of(glaParams)) %>%
  rename("Chla (ln ug.g-1)\n" = "chla_ug_g_1",
         "Glacier Coverage (ln proportion)" = "gl_cov_percent") %>%
  mutate(type = "specialist")

dat_4plotGen <- datGen %>%
  clean_names() %>%
  rename(site = "name") %>%
  select(site, gl_cov_percent, chla_ug_g_1, all_of(glaParams)) %>%
  rename(
    "Chla (ln ug.g-1)\n" = "chla_ug_g_1",
    "Glacier Coverage (ln proportion)" = "gl_cov_percent") %>%
  mutate(type = "generalist")

dat_all <- rbind(dat_4plotGen, dat_4plotSpec)

model_sel <- models_clean %>%
  filter(resp %in% glaParams) %>%
  filter(variables == "Chla (ln ug.g-1)\n")

data_sel <- dat_all %>%
  select(all_of(glaParams), "Chla (ln ug.g-1)\n", type) %>%
  pivot_longer(cols = chemotaxis:thiosulfate_oxidation_sox)%>%
  left_join(model_sel, join_by(type, name == resp))%>%
  mutate(color_pvalue = case_when(
    is.na(p_value) ~  "blue",
    p_value > 5e-2 ~ "blue",
    .default = "red"
  ))%>%
  mutate(good_name = case_when(
    name == "chemotaxis" ~ "Chemotaxis",
    name == "heterotrophy" ~ "Chemolitoheterotrophy\nOther",
    name == "mixotrophy_autotrophy_phototrophy" ~ "Mixotrophy\nautotroph/phototroph",
    name == "quorum_sensing" ~ "Quorum Sensing",
    name == "thiosulfate_oxidation_sox" ~ "Thiosulfate Oxidation",
    .default = name
  ))

# dat_text <- data_sel %>%
#   select(type, p_value, name) %>%
#   distinct() %>%
#   mutate(p_value = case_when(
#     is.na(p_value) ~  ">0.05",
#     p_value < 1e-2 ~ scientific(p_value, digits = 2),
#     .default = as.character(round(p_value, digits = 3))
#   )) %>%
#   mutate(p_value = paste("pvalue =" , p_value, sep = " "))


p1 <- ggplot(data_sel, aes(x = `Chla (ln ug.g-1)\n`, y = value * 100)) +
  geom_point() +
  facet_grid(good_name ~ type, scales = "free")+
  # ylab(var) +
  geom_smooth(method = "gam", aes(color = color_pvalue)) +
  scale_color_identity(guide = "legend")+
  theme_classic() +
  labs(y = "Relative abundance of MAGs (%)")+
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position="none",
    text = element_text(colour = "black", size = 14),
    panel.grid.major = element_blank()
  )

p1


ggsave_fitmax(
  "Figures/Fig_4b_trophicvsParameters_SpecGen.pdf",
  p1,
  maxwidth = 10
)
