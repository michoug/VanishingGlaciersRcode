library(targets)
library(tarchetypes)

options(tidyverse.quiet = TRUE)

source("R/functions.R")
source("R/packages.R")



tar_plan(
  tar_file_read(
    checkm,
    "data/NOMIS_MAGs_checkm_stats.txt",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    tax,
    "data/NOMIS_MAGS_tax.tsv",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    stats,
    "data/NOMIS_MAGs_statsTable.txt",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    bacttree,
    "data/gtdbtk.bac120.classify_user_midpoint.tree",
    read.tree(!!.x)
  ),
  
  tar_file_read(
    bactdataset,
    "data/MAGS_Presence_Datasets.txt",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    covsum,
    "data/MAGs_cov_sum.txt",
    read_tsv(!!.x)
  ),

  tar_file_read(
    euktree,
    "data/tree_fast.tree",
    read.tree(!!.x)
  ),
  
  tar_file_read(
    euktaxmags,
    "data/metadata_input.tsv",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    euktaxref,
    "data/metadata.tsv",
    read_tsv(!!.x)
  ),
  
  tar_file_read(
    eukcov,
    "data/Eukaryotes_MAGs_cov_norm.txt",
    read_tsv(!!.x)
  ),

  tar_target(
      fig1a_data,
      plot_length_compl_cont(checkm, tax, stats)
  ),
  
  tar_target(
    fig1_plot,
    ggsave_fitmax("Figures/Fig_1a_NOMIS_MAGs_stats.pdf", fig1a_data, maxwidth = 10),
    format = "file_fast"
  ),
  
  tar_target(
    fig1b_data,
    plot_barplot_novel_tax(tax)
  ),
  
  tar_target(
    fig1b_plot,
    ggsave_fitmax("Figures/Fig_1b_barplotNovelTaxa.pdf", fig1b_data),
    format = "file_fast"
  ),
  
  tar_target(
    fig1c_data,
    plot_bacterial_tree(tax, bactdataset, covsum,bacttree)
  ),
  
  tar_target(
    fig1c_plot,
    ggsave_fitmax("Figures/Fig_1c_bacterialTree.pdf", fig1c_data,maxheight =  15),
    format = "file_fast"
  ),
  
  tar_target(
    fig2a_data,
    plot_euk_tree(euktaxmags, euktaxref, eukcov, euktree)
  ),
  
  tar_target(
    fig2a_plot,
    ggsave_fitmax("Figures/Fig_2a_EukaryoticTreeAbove30.pdf", fig2a_data,maxwidth = 10),
    format = "file_fast"
  )
)
