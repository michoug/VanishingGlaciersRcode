plot_length_compl_cont <- function(dat_mags, tax, stats){
  
  dat_plot <- dat_mags %>%
    left_join(tax)%>%
    left_join(stats)%>%
    mutate(Quality =if_else(Completeness >= 90 & Contamination <= 5,"High","Medium"))%>%
    mutate(LengthNorm2 = LengthNorm2/1e6)
  
  p <- ggboxplot(dat_plot, x = "Quality", y = c("GC","Completeness","LengthNorm2", "Contamination"),
                  color = "Quality", palette = "jco", title = "",
                  xlab = "",
                  ylab = c("GC content (%)","Completeness (%)","Genome Length (Mbp)", "Contamination (%)"))
  
  p$GC <- p$GC + ggpubr::rotate()
  p$Completeness <- p$Completeness + ggpubr::rotate()
  p$Contamination <- p$Contamination + ggpubr::rotate()
  p$LengthNorm2 <- p$LengthNorm2 + ggpubr::rotate()
  
  p1 <- ggarrange(plotlist = p, common.legend = TRUE)
  p1
}

plot_euk_tree <- function(tax_mags, tax_ref, cov,tree){
  
  tax_mags <- tax_mags %>%
    select(`Unique ID`, `Lower Taxonomy`)%>%
    mutate(type = if_else(grepl("^GC",`Unique ID`),"ref","mags"))
  
  taxToPlot <- tax_mags%>%
    filter(type == "mags")%>%
    filter(`Unique ID`%in% tree$tip.label)
  
  tax_ref <- tax_ref%>%
    select(`Unique ID`, `Lower Taxonomy`)%>%
    mutate(type = "ref")
  
  tax <- rbind(tax_mags, tax_ref)
  colnames(tax) <- c("MAGs","Tax", "type")
  
  tax <- tax %>%
    mutate(taxa_good = if_else(Tax %in% taxToPlot$`Lower Taxonomy`, Tax, "Other"))
  
  tax_temp <- tax %>%
    filter(MAGs %in% tree$tip.label)
  
  cov$MAGs <- gsub("_","-",cov$MAGs)
  
  cov_max <- cov %>%
    rowwise()%>%
    summarise(total = sum(c_across(where(is.numeric))))
  
  cov_max$MAGs <- cov$MAGs
  
  cov_max_all <- cov_max %>%
    filter(MAGs %in% tree$tip.label)%>%
    right_join(tax, multiple = "all")%>%
    select(MAGs, total)%>%
    # mutate_all(~replace_na(.,1))%>%
    distinct()
  
  cov_max_all <- as.data.frame(cov_max_all)
  
  rownames(cov_max_all) <- cov_max_all$MAGs
  cov_max_all$MAGs <- NULL
  cov_max_all$count <- log10(cov_max_all$total)
  cov_max_all$total <- NULL
  
  meta	<- split(tax$MAGs, tax$taxa_good)
  
  tree_meta	<- groupOTU(tree, meta)
  
  getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))
  
  color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
             "#FF7F00", "#FFFF33", "#A65628", "black", "#999999")
  
  color[1]	<- "black"
  
  p1	<-
    ggtree(tree_meta,
           layout = 'circular',
           aes(color = group),
           #branch.length = "none"
    ) + #
    geom_tree() +
    theme_tree() +
    # geom_tiplab()+
    geom_treescale(width	= 0.1) +
    scale_color_manual(values	= color,
                       na.value = "transparent",
                       guide = "none") +
    theme(legend.position = "right")+
    new_scale_colour()+
    new_scale_fill()
  
  p2 <- p1 +
    new_scale_colour()+
    new_scale_fill()+
    geom_fruit(
      data=tax,
      pwidth	= 0.05,
      geom=geom_bar,
      mapping=aes(y=MAGs, x = 4, fill = taxa_good),
      orientation="y",
      stat="identity"
    )+
    scale_fill_manual(values	= color[-1])+
    labs(fill = "Taxa")
  
  # p2
  
  p3 <- p2 %<+% tax_temp +
    geom_tippoint(aes(color = type))+
    new_scale_colour()+
    new_scale_fill()
  
  p4 <- gheatmap(p3, cov_max_all, offset=0.2, width=0.05,
                 colnames=FALSE,color = NULL)+
    scale_fill_viridis_c()+
    labs(fill = "Normalized log10\nabundance")
  
  p4
}

plot_barplot_novel_tax <- function(dat){
  dat2Plot <- dat %>%
    pivot_longer(cols = !MAGs)%>%
    mutate(value = gsub(".__","", value))%>%
    mutate(status =if_else(value == "", "UnKnown", "Known")) %>%
    select(-MAGs)%>%
    group_by(name, status)%>%
    summarise(sum = n(), perc = sum/2868, label = percent(perc, accuracy = 0.1))%>%
    filter(!(name %in% c("d","p","s")))%>% # Remove taxa levels if there is no unknown
    mutate(label = if_else(status == "UnKnown",paste(sum, " MAGs (", label, ")", sep = ""),""))%>%
    mutate(name = case_when(
      name == "c" ~ "Class",
      name == "o" ~ "Order",
      name == "g" ~ "Genus",
      name == "f" ~ "Family"
      # .default = "other".default = "other"
    ))%>%
    mutate(name = factor(name, levels = c("Class", "Order", "Family", "Genus")))
  
  
  p1 <- ggplot(dat2Plot, aes(x = perc, y = reorder(name, desc(name)), fill = status, label = label))+
    geom_bar(stat = "identity", color = "black")+
    geom_text(nudge_x = 0.25)+
    scale_fill_manual(values = c("Known" = "white", "UnKnown" = "black"), labels = c("Known taxa", "Novel taxa"))+
    scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
    xlim(0,1)+
    labs(fill = NULL, x = "Percentage of MAGs (%)", y = NULL)+
    theme_classic()+
    theme(legend.position="top")
  p1
}

plot_bacterial_tree <- function(dat,dataset,covmax,BacTree){
  
  bacDat <- dat %>%
    mutate(p_c = if_else(p == "p__Proteobacteria", c, p))%>%
    mutate(p_c = gsub(".__","",p_c, perl = T))%>%
    mutate(p_c = gsub("_.$","",p_c, perl = T))%>%
    filter(d == "d__Bacteria")%>%
    mutate(Abundance	= 0.2)
  
  mostAbundantTax	<- bacDat%>%
    group_by(p_c)%>%
    summarise(total	= n())%>%
    arrange(desc(total))%>%
    slice(1:16)
  
  list <- mostAbundantTax$p_c
  list <- c(list, "Nitrospirota")
  
  
  bacDat	<- bacDat%>%
    mutate(p_c	= if_else(p_c %in% list,  p_c, "Other"))
  
  bacDatset <- dataset%>%
    filter(MAGs %in% bacDat$MAGs)
  
  bacDatset$Busi <- gsub("Busi", "Busi et al., Nat Com, 2022", bacDatset$Busi)
  bacDatset$ENSEMBLE <- gsub("ENSEMBLE", "Michoud et al, L&O, 2023", bacDatset$ENSEMBLE)
  bacDatset$Tibet <- gsub("Tibet", "Tibetan Glacier Genome and Gene", bacDatset$Tibet)
  bacDatset$Tara <- gsub("Tara", "Tara Oceans", bacDatset$Tara)
  
  bacDatset <- as.data.frame(bacDatset)
  rownames(bacDatset) <- bacDatset$MAGs
  bacDatset$MAGs <- NULL
  
  bactcov <- covmax%>%
    filter(MAGs %in% bacDat$MAGs)%>%
    select(MAGs, count)
  
  bactcov <- as.data.frame(bactcov)
  rownames(bactcov) <- bactcov$MAGs
  bactcov$MAGs <- NULL
  
  bactcov$count <- log10(bactcov$count)
  
  a	<- split(bacDat$MAGs, bacDat$p_c)
  
  tree	<- groupOTU(BacTree, a)
  
  getPaletteBact = colorRampPalette(brewer.pal(9, "Set1"))
  
  bactColor	<- getPaletteBact(length(unique(bacDat$p_c)) + 1)
  
  bactColor[1]	<- "black"
  
  
  p1	<-
    ggtree(tree,
           layout = 'circular',
           aes(color = group),
           #branch.length = "none"
    ) + #
    geom_tree() +
    theme_tree() +
    geom_treescale(width	= 0.1) +
    scale_color_manual(values	= bactColor,
                       na.value = "transparent",
                       guide = "none") +
    #geom_text2(aes(subset=!isTip,	label=node), hjust=-.3)+
    theme(legend.position = "right")
  # p1

  p2 <- p1 +
    new_scale_colour()+
    new_scale_fill()+
    geom_fruit(
      data=bacDat,
      pwidth	= 0.01,
      geom=geom_bar,
      mapping=aes(y=MAGs, fill = p_c,x = 1),
      # orientation="y",
      stat="identity",
    )+
    scale_fill_manual(values	= bactColor[-1])+
    labs(fill = "Taxa")+
    new_scale_colour()+
    new_scale_fill()

  # p2
  
  p3 <- gheatmap(p2, bacDatset,width = 0.2,offset = 0.1,# offset=8, width=0.6,
                 colnames=FALSE, color = NULL) +
    scale_fill_discrete(na.translate = F)+
    labs(fill = "Datasets")+
    new_scale_colour()+
    new_scale_fill()
  # p3
  # 
  p4 <- gheatmap(p3, bactcov, offset=0.6, width=0.05,
                 colnames=FALSE, color = NULL)+
    scale_fill_viridis_c()+
    labs(fill = "Normalized log10\nabundance")

  p4
}

ggsave_fitmax <- function(filename,
                          plot,
                          maxheight = 7,
                          maxwidth  = maxheight,
                          units     = "in", ...) {
  if(is.null(plot)) return(FALSE)
  dims = get_dims(ggobj     = plot,
                  maxheight = maxheight,
                  maxwidth  = maxwidth,
                  units     = units)
  ggplot2::ggsave(filename = filename,
                  plot   = plot,
                  height = dims$height,
                  width  = dims$width,
                  units  = units, ...)
}

get_dims <- function(ggobj,
                     maxheight,
                     maxwidth = maxheight,
                     units    = "in", ...) {
  
  # Internal helper function:
  # Treat all null units in a unit object as if they were inches.
  # This is a bad idea in gneral, but I use it here as a workaround.
  # Extracting unit names from non-atomic unit objects is a pain,
  # so questions like "which rows of this table layout have null heights?"
  # are hard to answer. To work around it, I exploit an (undocumented!)
  # quirk: When calculating the size of a table layout inside a Grid plot,
  # convertUnit(...) treats null units as zero.
  # Therefore
  #	(convertHeight(grob_height, "in", valueOnly=TRUE)
  #	- convertHeight(null_as_if_inch(grob_height), "in", valueOnly=TRUE))
  # does the inverse of convertUnit: It gives the sum of all *null* heights
  # in the object, treating *fixed* units as zero.
  #
  # Warning: I repeat, this approach ONLY makes any sense if
  #	convertUnit(unit(1, "null"), "in", "x", valueOnly=T) == 0
  # is true. Please check that it is before calling this code.
  .null_as_if_inch = function(u){
    stopifnot(packageVersion("grid") < "4.0")
    if(!grid::is.unit(u)) return(u)
    if(is.atomic(u)){
      if("null" %in% attr(u, "unit")){
        d = attr(u, "data")
        u = unit(
          x=as.vector(u),
          units=gsub("null", "in", attr(u, "unit")),
          data=d)
      }
      return(u)
    }
    if(inherits(u, "unit.arithmetic")){
      l = .null_as_if_inch(u$arg1)
      r = .null_as_if_inch(u$arg2)
      if(is.null(r)){
        args=list(l)
      }else{
        args=list(l,r)
      }
      return(do.call(u$fname, args))
    }
    if(inherits(u, "unit.list")){
      return(do.call(grid::unit.c, lapply(u, .null_as_if_inch)))
    }
    return(u)
  }
  
  if(inherits(ggobj, "ggplot") && !isTRUE(ggobj$respect) &&
     is.null(ggobj$theme$aspect.ratio) && is.null(ggobj$coordinates$ratio) &&
     is.null(ggplot2::theme_get()$aspect.ratio)) {
    return(list(height = maxheight, width = maxwidth))
  }
  
  tmpf = tempfile(pattern = "dispos-a-plot", fileext = ".png")
  png(filename = tmpf,
      height   = maxheight,
      width    = maxwidth,
      units    = units,
      res      = 120, ...)
  
  on.exit({
    dev.off()
    unlink(tmpf)
  })
  
  if (inherits(ggobj, "ggplot")) {
    g = ggplot2::ggplotGrob(ggobj)
  } else if (inherits(ggobj, "gtable")) {
    g = ggobj
  } else {
    stop("Don't know how to get sizes for object of class ", deparse(class(ggobj)))
  }
  
  stopifnot(grid::convertUnit(grid::unit(1, "null"), "in", "x", valueOnly = TRUE) == 0)
  known_ht = sum(grid::convertHeight(g$heights, units, valueOnly = TRUE))
  known_wd = sum(grid::convertWidth(g$widths,   units, valueOnly = TRUE))
  free_ht  = maxheight - known_ht
  free_wd  = maxwidth  - known_wd
  
  if (packageVersion("grid") >= "4.0.0") {
    null_rowhts <- as.numeric(g$heights[grid::unitType(g$heights) == "null"])
    null_colwds <- as.numeric(g$widths[grid::unitType(g$widths) == "null"])
    panel_asps <- (
      matrix(null_rowhts, ncol = 1)
      %*% matrix(1 / null_colwds, nrow = 1))
  } else {
    all_null_rowhts <- (
      grid::convertHeight(.null_as_if_inch(g$heights), "in", valueOnly = TRUE)
      - grid::convertHeight(g$heights, "in", valueOnly = TRUE))
    all_null_colwds <- (
      grid::convertWidth(.null_as_if_inch(g$widths), "in", valueOnly = TRUE)
      - grid::convertWidth(g$widths, "in", valueOnly = TRUE))
    null_rowhts <- all_null_rowhts[all_null_rowhts > 0]
    null_colwds <- all_null_colwds[all_null_colwds > 0]
    panel_asps <- (matrix(null_rowhts, ncol = 1) %*% matrix(1 / null_colwds, nrow = 1))
  }
  
  panel_asps = matrix(null_rowhts, ncol = 1) %*% matrix(1 / null_colwds, nrow = 1)
  max_rowhts = free_ht/sum(null_rowhts) * null_rowhts
  max_colwds = free_wd/sum(null_colwds) * null_colwds
  rowhts_if_maxwd = max_colwds[1] * panel_asps[, 1]
  colwds_if_maxht = max_rowhts[1]/panel_asps[1, ]
  height = min(maxheight, known_ht + sum(rowhts_if_maxwd))
  width  = min(maxwidth,  known_wd + sum(colwds_if_maxht))
  return(list(height = height, width = width))
}
