## Modified from https://github.com/IndrajeetPatil/ggstatsplot/issues/654

library("QsRutils")
library(statsExpressions)
library(ggstatsplot)

########### Get the letters to display:
CalcMultCompLetters <- function(data){
  
  #data is the output from pairwise_comparison function. It will extract the contrasts and the p-vals and make a dataframe with mult comp letters.
  
  pvals = data$p.value
  names(pvals) = paste0(data$group1, "-",data$group2)
  
  results <- multcompView:::multcompLetters(pvals)
  return(results)
}

########## Add compact letter display to plot:
AddLetters <- function(data, plot, x, y,ytext,
                          type = "parametric", tr = .2, paired = FALSE, var.equal = FALSE, p.adjust.method = "holm", k=2L){
  
  # data is the original data.frame
  # plot is the object plot called from ggbetweenstats and saved in the env. remember to set `pairwise.comparisons = FALSE`
  # Additional options are for the pairwise comparisons. See `pairwiseComparisons::pairwise_comparisons` for details
  
  #---------------- Run pairwise comparisons:
  
  # Run multiple comparisons or not? depends on the no. of levels of the independent variable
  if (nlevels(factor(data %>% dplyr::pull({{ x }}))) == 1){
    cat("\nInput df has only one level, no pairwise comparisons can be made. Returning the original plot...\n")
    return(plot)
  }
  
  # Run pairwise_comparisons:
  plot_stats <- pairwise_comparisons(
    data = data,
    x = {{ x }},
    y = {{ y }},
    type = type,
    tr = tr,
    paired = paired,
    var.equal = var.equal,
    p.adjust.method = p.adjust.method,
    k = k
  )
  
  # Get the letters to display:
  letters <- CalcMultCompLetters(plot_stats)

  
  #---------------- Plot:
  
  #1. Get the top of each box (in terms of y-coordinates):
  df = as.data.frame(plot$data[,c(1:2)])
  box.rslt <- graphics::boxplot(df[,2] ~ factor(df[,1]), plot = FALSE, data = df)
  
  #2. Make a df with letters and the top of the boxes:
  x <- c(1:length(letters$Letters))
  y <- box.rslt$stats[5, ]
  cbd <- letters$Letters
  ltr_df <- data.frame(x, y, cbd)
  
  #3. If we plot the CLDs at the coordinates in ltr_df, they will over plot the tops of the whiskers . We need to nudge the CLDs upward to avoid the overlap. To determine how much to nudge, I will get the range of the Y-axis and nudge upward 5% of this range.
  lmts <- QsRutils::get_plot_limits(plot)
  y.range <- lmts$ymax - lmts$ymin
  y.nudge <- 1 * y.range
  
  #4. Update the plot:
  plot = plot +
    geom_text(data = ltr_df, aes(x=x, y=ytext, label=cbd))
  
  return(plot)
  
}


### Example:
#library(ggstatsplot)
#library(pairwiseComparisons)
#library(ggplot2)
#plot = ggbetweenstats(mtcars, cyl, mpg, pairwise.comparisons = FALSE)
#plot_stat = pairwise_comparisons(mtcars, cyl, mpg)
#AddLetters(mtcars, plot = plot, x= cyl, y=mpg)




############### Add compact letter display - grouped function:
grouped_AddLetters <- function(data_list, plot_list, x, y, excluded_plots = NULL,
                               type = "parametric", tr = .2, paired = FALSE, var.equal = FALSE, p.adjust.method = "holm", k=2L){
  
  #data_list is a list of dataframes created as per example here: https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html#ggbetweenstats-1
  #plot_list is a list of ggstatsplot created as per example at the same link as above.
  #excluded_plots is a vector of df names (i.e. the names of each df in `data_list`) that you do NOT want to be changed by adding the compact letter display.
  # Additional options are for the pairwise comparisons. See `pairwiseComparisons::pairwise_comparisons` for details
  
  # Get names of plots:
  nomi = names(data_list)
  
  # Do not add letters to excluded plots:
  if (!(is.null(excluded_plots))){
    `%!in%` <- purrr::compose(`!`, `%in%`) #Make a "not in" function
    keep = nomi %!in% excluded_plots
    nomi = nomi[keep]
  }
  
  # Select each one plot at the time and update it
  for (plot_number in nomi){
    data = data_list[[plot_number]]
    plot = plot_list[[plot_number]]
    
    updatedplot = AddLetters(data = data, plot = plot, x={{x}}, y={{y}},type = type, tr = tr, paired = paired, var.equal = var.equal, p.adjust.method = p.adjust.method, k=k) #To understand the usage of `{{ }}` see: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html#fn2
    plot_list[[plot_number]] = updatedplot
  }
  return(plot_list)
}

### Example:
#library(ggstatsplot)
#library(pairwiseComparisons)
#library(ggplot2)
#library(gapminder)
#year_list <-
#  gapminder::gapminder %>%
#  dplyr::filter(year %in% c(1967, 1987, 2007), continent != "Oceania") %>%
#  split(f = .$year, drop = TRUE)

#plot_list <-
#  purrr::pmap(
#    .l = list(
#      data = year_list,
#      x = "continent",
#      y = "lifeExp",
#      title = list(
#        "Year: 1967",
#        "Year: 1987",
#        "Year: 2007"
#      ),
#      type = list("r", "bf", "np"),
#      pairwise.comparisons = list("FALSE","FALSE","FALSE"),
#      plot.type = list("box", "boxviolin", "violin")
#    ),
#    .f = ggstatsplot::ggbetweenstats
#  )

#updated_plot = grouped_AddLetters(year_list, plot_list = plot_list,  x=continent, y=lifeExp, excluded_plots = c("1987","1967"))

#ggstatsplot::combine_plots(
#  plotlist = updated_plot,
#  annotation.args = list(title = "Changes in life expectancy across continents (1967-2007)"),
#  plotgrid.args = list(ncol = 3)
#)
