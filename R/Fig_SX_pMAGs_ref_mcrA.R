library(tidyverse)
library(gggenomes)
library(ggpubr)

source("customFunctions/plot_functions.R")

genome_11 <- "GFS_11005_182"

ids_genome_11 <- c(
  "GFS_11005_182",
  "NC_015416.1_subset",
  "NC_008553.1_subset",
  "NC_003552.1_subset",
  "NC_017527.1_subset"
)

genome_32 <- "GFS_3223_91"

ids_genome_32 <- c(
  "NZ_JQLR01000001.1_subset",
  "NZ_CP009516.1_subset",
  "NZ_CP032683.1_subset",
  "NZ_CP009515.1_subset",
  "NC_013665.1_subset"
)

genes_to_plot <- c(
  "mcrA","mcrB","mcrC","mcrD","mcrG",
  "topA", "pelA"
)

## Functions

rename_genes <- function(table, id){
  table_rename <- table %>%
    mutate({{id}} := case_when(
      get(id) == "NC_003901.1_subset" ~ "Methanosarcina_mazei_Go1",
      get(id) == "NC_003552.1_subset" ~ "Methanosarcina_acetivorans_C2A",
      get(id) == "NC_013665.1_subset" ~ "Methanocella_paludicola_SANAE",
      get(id) == "NC_008553.1_subset" ~ "Methanothrix_thermoacetophila_PT",
      get(id) == "NC_009464.1_subset" ~ "Methanocella_arvoryzae_MRE50",
      get(id) == "NC_015416.1_subset" ~ "Methanothrix_soehngenii_GP6",
      get(id) == "NC_017527.1_subset" ~ "Methanothrix_harundinacea_6Ac",
      get(id) == "NC_017034.1_subset" ~ "Methanocella_conradii_HZ254",
      get(id) == "NZ_JQLR01000001.1_subset" ~ "Methanosarcina_soligelidi_SMA_21",
      get(id) == "NZ_CP009501.1_subset" ~ "Methanosarcina_thermophila_TM_1",
      get(id) == "NZ_CP009520.1_subset" ~ "Methanosarcina_vacuolata_Z_761",
      get(id) == "NZ_CP009530.1_subset" ~ "Methanosarcina_barkeri_227",
      get(id) == "NZ_CP009506.1_subset" ~ "Methanosarcina_siciliae_T4_M",
      get(id) == "NZ_CP009515.1_subset" ~ "Methanosarcina_lacustris_Z_7289",
      get(id) == "NZ_CP009516.1_subset" ~ "Methanosarcina_horonobensis_HB_1",
      get(id) == "NZ_CP032683.1_subset" ~ "Methanosarcina_flavescens_E03.2",
      get(id) == "NZ_LMVP01000517.1_subset" ~ "Methanosarcina_spelaei_MC_15",
      .default = get(id)
    ))
}

plot_synteny <- function(genome, ids_sel, genes_to_plot) {
  
  ids_sel <- c(genome, ids_sel)
  
  meth_links <- read_paf("data/pMAGs_ref_mcrA_contigs.paf") %>%
    filter(seq_id %in% ids_sel) %>%
    filter(seq_id2 %in% ids_sel)
  
  meth_links <- rename_genes(meth_links, "seq_id")
  meth_links <- rename_genes(meth_links, "seq_id2")
  
  order <- meth_links %>%
    filter(seq_id == genome) %>%
    filter(seq_id2 != genome) %>%
    arrange(desc(map_length))
  
  order_final <- factor(c(genome, order$seq_id2),
                        levels = c(genome, order$seq_id2))
  
  meth_seq <- read_seq_len("data/pMAGs_ref_mcrA_contigs.fasta") %>%
    filter(seq_id %in% c(genome,ids_sel)) 
  
  meth_seq <- rename_genes(meth_seq,"seq_id")%>%
    mutate(seq_id = factor(seq_id, order_final))
  
  meth_seq$seq_id <- meth_seq$seq_id[order(meth_seq$seq_id)]
  
  meth_egg <- read_tsv("data/pMAGs_ref_mcrA_contigs_annot.tsv", comment = "##") %>%
    select("#query", Preferred_name) %>%
    rename("Genes" = "#query") %>%
    filter(Preferred_name %in% genes_to_plot)
  
  
  meth_genes <- read_gff3("data/pMAGs_ref_mcrA_contigs.gff") %>%
    filter(seq_id %in% ids_sel) %>%
    mutate(gc_cont = as.numeric(gc_cont)) %>%
    mutate(Genes = gsub("\\d+_(\\d+)", "\\1", feat_id)) %>%
    mutate(Genes = paste(seq_id, Genes, sep = "_")) %>%
    left_join(meth_egg, join_by(Genes))
  
  meth_genes <- rename_genes(meth_genes,"seq_id")
  
  p1 <- gggenomes(seqs = meth_seq,
                  genes = meth_genes,
                  links = meth_links) %>%
    # add_sublinks(meth_links) %>%
    sync() +
    geom_seq() +
    geom_bin_label() +
    geom_gene(aes(fill = Preferred_name)) +
    geom_link() +
    geom_gene_tag(aes(label = Preferred_name),
                  nudge_y = 0.1,
                  check_overlap = TRUE) +
    scale_fill_brewer("Genes", palette = "Dark2", na.value = "cornsilk3")
  
  p1
}

## Plotting


p1 <- plot_synteny(genome_11, ids_genome_11, genes_to_plot)
p2 <- plot_synteny(genome_32, ids_genome_32, genes_to_plot)

p <- ggarrange(p1, p2, common.legend = T, ncol = 1, labels = "auto")

ggsave_fitmax("Figures/Fig_sX_Methanogens.pdf", p, maxwidth = 14)




