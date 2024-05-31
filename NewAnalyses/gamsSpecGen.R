library(mgcv)
library(tidyverse)
library(performance)
library(janitor)
# library(vegan)
library(janitor)

setwd("~/../switchdrive/Institution/NOMIS_MAGs/Correlations/")

meta <- read_tsv("../metadata_NOMIS_sed_chem_noNAs.txt")
dat_prok <- read_tsv("Prok4CorrelGeneralist.txt")
dat_euk <- read_tsv("Euk4Correl.txt")
dat_vir <- read_tsv("Virus4Correl.txt")
dat_AMG <- read_tsv("AMGs4Correl.txt")

meta_correl <- meta %>%
  clean_names()%>%
  # select_if(~ !any(is.na(.)))%>%
  select(c(sample, site_c,water_temp_c,
           ph_p_h, do_mg_l_1,
           do_sat_saturation,
           conductivity_u_s_cm_1,
           turb_ntu, lat_sn_dd,
           lon_sn_dd, ele_sn_m,
           max_el_m, mean_el_m,
           sn_sp_dist_m, sn_sp_ele_m,
           gl_sa_km2, gl_cov_percent,
           gl_a_km2, chla_ug_g_1,
           n3_srp_ug_l_1, n4_nh4_ug_l_1, n5_no3_ug_l_1, n6_no2_ug_l_1,
           sba_cells_g_1))

log_const <- function(x){return(log(x + (min(x[which(x > 0)])/2)))}
chem <- c("water_temp_c","ph_p_h","do_mg_l_1","do_sat_saturation","conductivity_u_s_cm_1",
          "turb_ntu","ele_sn_m","max_el_m","mean_el_m","sn_sp_dist_m",
          "sn_sp_ele_m","gl_sa_km2","gl_cov_percent",
          "gl_a_km2","chla_ug_g_1","sba_cells_g_1","n3_srp_ug_l_1", "n4_nh4_ug_l_1", "n5_no3_ug_l_1", "n6_no2_ug_l_1" )

for(i in chem){
  meta_correl[[i]] <- log_const(meta_correl[[i]]) 
}


dat_prok_sel <- dat_prok %>%
  select(-c(`Dicarboxylate-hydroxybutyrate cycle`,
            `Hydroxypropionate-hydroxybutylate cycle`))

dat <- dat_prok_sel %>%
  # select(name,contains("trop"), contains("MAG"))%>%
  left_join(dat_euk)%>%
  left_join(dat_vir)%>%
  left_join(meta_correl, by = join_by(name == sample))
  # inner_join(meta_correl, by = join_by(variable == sample)) %>%
  # rename(name = variable)


dat$glacier <- gsub("_.*", "", dat$name)
dat_sel <- dat %>%
  clean_names()%>%
  select(-sba_cells_g_1)%>%
  select(-do_mg_l_1)%>%
  # select(-matches("n\\d"))%>%
  na.omit()

response <- colnames(dat_sel)[2:60]

options(warn = 0)

# response <- response[1:5]
for(resp in response){
  expl <- c("water_temp_c","ph_p_h","conductivity_u_s_cm_1",
            "turb_ntu","ele_sn_m","sn_sp_dist_m","gl_sa_km2",
            "gl_cov_percent","chla_ug_g_1", 
            "n3_srp_ug_l_1","n4_nh4_ug_l_1","n5_no3_ug_l_1","n6_no2_ug_l_1")
  
  
  for(j in 1:(length(expl)+1)){
    
    formula <- paste(resp, "~ s(lat_sn_dd, lon_sn_dd, bs=\'sos\', k=-1, m=1)")
    
    if(length(expl) > 0){
      for(i in expl){
        # formula <-  paste(formula, " + s(",i,",bs = 'cr', k = 3)")
        formula <-  paste(formula, "+", i)
      }
    }
    
    # ctrl <- lmeControl(opt='optim')
    
    gam <- gamm(data = dat_sel, formula = as.formula(formula), 
                # correlation = corAR1(form = ~ 1 | glacier), 
                method = 'REML')
    name <- paste("gam_",j, sep = "")
    assign(name,gam)
    summary <- summary(gam$gam)
    summary_f <- as.data.frame(summary$p.table)
    
    summary_f <- summary_f[!(row.names(summary_f)) == "(Intercept)",]
    
    if(length(row(summary_f))== 0){break}
    
    minRow <- summary_f %>%
      filter(`t value` == min(`t value`))
    
    val <- rownames(minRow)
    # val <- gsub("s\\((.*)\\)", "\\1", val, perl = T)
    # print(val)
    expl <- expl[expl != val]
  }
  
  
  
  BFdat <- test_bf(gam_14$gam, gam_13$gam, gam_12$gam,gam_11$gam,
                   gam_10$gam, gam_9$gam,
                   gam_8$gam, gam_7$gam, gam_6$gam, gam_5$gam,
                   gam_4$gam, gam_3$gam, gam_2$gam, gam_1$gam)
  
  BFdf <- as.data.frame(BFdat)
  
  BFdf$gam <- c("gam_14$gam", "gam_13$gam","gam_12$gam","gam_11$gam",
                "gam_10$gam","gam_9$gam","gam_8$gam",
                "gam_7$gam","gam_6$gam","gam_5$gam",
                "gam_4$gam","gam_3$gam","gam_2$gam","gam_1$gam")
  
  BFdf <- BFdf %>%
    na.omit()%>%
    filter(BF == max(BF))
  
  goodgam <- paste(BFdf$gam)
  goodgam <- gsub("\\$gam","", goodgam)
  
  final_summary <- as.data.frame(summary(get(goodgam)$gam)$p.table)
  final_summary$resp <- resp
  final_summary$variables <- rownames(final_summary)
  
  name_resp <- paste("ModelGen/",resp, "_model.txt", sep = "")
  write_tsv(final_summary, name_resp)
  print(resp)
}



