## Copied from https://github.com/kevinwolz/hisafer

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
