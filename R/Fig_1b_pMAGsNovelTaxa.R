## ---------------------------
##
## Script name: barplotNovelTaxa.R
##
## Purpose of script: Create a barplot showing for each taxonomical level, the number of known and unknowns taxa
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


library(tidyverse)
library(ggplot2)
library(scales)

dat	<- read_tsv("data/pMAGs_tax.tsv", col_names = T)

dat2Plot <- dat %>%
  pivot_longer(cols = !MAGs) %>%
  mutate(value = gsub(".__", "", value)) %>%
  mutate(status = if_else(value == "", "UnKnown", "Known")) %>%
  select(-MAGs) %>%
  group_by(name, status) %>%
  summarise(
    sum = n(),
    perc = sum / 2868,
    label = percent(perc, accuracy = 0.1)
  ) %>%
  filter(!(name %in% c("d", "p", "s"))) %>% # Remove taxa levels if there is no unknown
  mutate(label = if_else(status == "UnKnown", paste(sum, " MAGs (", label, ")", sep = ""), "")) %>%
  mutate(
    name = case_when(
      name == "c" ~ "Class",
      name == "o" ~ "Order",
      name == "g" ~ "Genus",
      name == "f" ~ "Family"
      # .default = "other".default = "other"
    )
  ) %>%
  mutate(name = factor(name, levels = c("Class", "Order", "Family", "Genus")))


p1 <-
  ggplot(dat2Plot, aes(
    x = perc,
    y = reorder(name, desc(name)),
    fill = status,
    label = label
  )) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(nudge_x = 0.25) +
  scale_fill_manual(
    values = c("Known" = "white", "UnKnown" = "black"),
    labels = c("Known taxa", "Novel taxa")
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlim(0, 1) +
  labs(fill = NULL, x = "Percentage of MAGs (%)", y = NULL) +
  theme_classic() +
  theme(legend.position = "top")

ggsave_fitmax("Figures/Fig_1b_barplotNovelTaxa.pdf", p1)
