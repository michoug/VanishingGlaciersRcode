library(vegan)
library(tidyverse)
library(janitor)
library(linkET)

set.seed(123)

source("customFunctions/plot_functions.R")

map <- read_tsv("data/metadata_NOMIS_sed_chem.txt")
covMAGs <- read_tsv("data/pMAGs_cov_norm.txt.gz")
KO4Cor <- read_tsv("data/pMAGs_KO_modules.tsv.gz")

goodNames <- c("Chla μg.g-1"="chla_ug_g_1",
               "Conductivity (μs.cm-1)" = "conductivity_u_s_cm_1",
               "DO (mg l-1)" = "do_mg_l_1",
               "Galcier Coverage (%)" = "gl_cov_percent",
               "SO4 (μg.l-1)" = "i6_so4_ug_l_1",
               "SRP (μg.l-1)" = "n3_srp_ug_l_1",
               "NO3 (μg.l-1)" = "n5_no3_ug_l_1",
               "pH" = "ph_p_h",
               "Turbity (NTU)" = "turb_ntu",
               "Water temperature °C" = "water_temp_c"
               )

map_clean <- map %>%
  clean_names()%>%
  select("sample",
         "code_gl",
         "site_c",
         "water_temp_c",
         "ph_p_h",
         "do_mg_l_1",
         "conductivity_u_s_cm_1",
         "turb_ntu",
         "gl_cov_percent",
         "chla_ug_g_1",
         "i6_so4_ug_l_1",
         "n3_srp_ug_l_1",
         "n5_no3_ug_l_1",
         )%>%
  column_to_rownames(var = "sample")%>%
  select(code_gl, site_c, sort(colnames(.)))%>%
  rename(all_of(goodNames))

KO4Cor <- KO4Cor %>%
  select(rownames(map_clean))

KO4Cor_t <- as.data.frame(t(KO4Cor))

correlationMap <- cor(map_clean[3:12], method = "spearman", use = "complete.obs")

covMAGsPrep <- covMAGs %>%
  select(MAGs,rownames(map_clean))%>%
  column_to_rownames(var = "MAGs")

covMAGsPrep <- as.data.frame(t(covMAGsPrep))

correlationMAGs <- cor(covMAGsPrep, method = "spearman")

fullsepc <- cbind(covMAGsPrep,KO4Cor_t)

mantel <- mantel_test(fullsepc,map_clean[3:12],
                      spec_select = list(tax = 1:2855,
                                          ko_modules = 2856:3264
                      #                    specialist = 902:2855)
                      ))%>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


p1 <- qcorrplot(correlationMap,type = "lower",diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))


ggsave_fitmax("Figures/Fig_Xb_mantel.pdf",p1, maxheight = 10)
