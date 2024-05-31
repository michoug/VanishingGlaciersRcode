library(tidyverse)
library(vegan)
library(ggvegan)

tad90 <- read_tsv("NewAnalyses/MAGs_TAD90.txt")
tad80 <- read_tsv("NewAnalyses/MAGs_TAD80.txt")

tad_clean <- function(tad, number){
  dat_clean <- tad %>%
  pivot_longer(-Genome) %>%
  mutate(name = gsub("$", paste("_",number, sep = ""), name))
}

tad90_clean <- tad_clean(tad90, 90)
tad80_clean <- tad_clean(tad80, 80)

tad_all <- rbind(tad90_clean, tad80_clean)

tad_all_clean <- tad_all %>%
  mutate(name = gsub("MAGs.fa.(.*)_mg.r1.preprocessed.fq.gz Trimmed Mean", "\\1", name)) %>%
  pivot_wider(names_from = name, values_from = value)%>%
  select(-Genome) %>%
  mutate_if(is.double, function(x, na.rm = FALSE)
    (x * 100)) %>%
  mutate_if(is.double, as.integer)

nmds_dist  <- vegdist(t(tad_all_clean), method = "bray")

dat_nmds <- metaMDS(t(tad_all_clean))

dat_fort <- as.data.frame(dat_nmds$points)
dat_fort$type <- gsub(".*_(\\d\\d)$", "\\1", rownames(dat_fort))

p <- ggplot(dat_fort, aes(x = MDS1, y = MDS2, colour = type)) +
  geom_point() +
  theme_classic() +
  # theme(legend.position = "none") +
  labs(x = "NMDS1", y = "NMDS2") +
  ggtitle("Comparisons of coverage methods")

p

ggsave("NewAnalyses/nmds_tads.pdf", p, width = 10, height = 10)

all(rownames(dat_fort) == labels(nmds_dist))

dat.ano <- with(dat_fort, anosim(nmds_dist, type, distance = "bray"))
dat.ano

# anosim(x = nmds_dist, grouping = type, distance = "bray")
# Dissimilarity: bray
# 
# ANOSIM statistic R: -0.002012
# Significance: 0.794
# 
# Permutation: free
# Number of permutations: 999