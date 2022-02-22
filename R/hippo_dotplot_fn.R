## input should have 2 slots, input$Cluster_mean and input$Cluster_zero_proportion
## input: a list object with 2 items inside as selected.genes, one slot name is 'Cluster_mean', another slot name is 'Cluster_zero_proportion'
## features or topN: indicate either features to be selected, or first 'topN' column features to be selected
## cols: color options.
## col.min: Minimum scaled average expression threshold (everything smaller will be set to this), default = -2.5
## col.max:	Maximum scaled average expression threshold (everything larger will be set to this), default = 2.5
## dot.min: The fraction of zero proportion at which to draw the smallest dot (default is 0), ranging from 0 to 1. Genes with less than this value will have no dot drawn.
## dot.scale: Scale the size of the points,
## scale: Determine whether the data is scaled, TRUE for default
## scale.by: Scale the size of the points by 'size' or by 'radius', by default points are sclaed by 'size'
## scale.min: Set lower limit for points scaling to indicate the percentage of zero proportion, by default decide by data automatically
## scale.max: Set upper limit for points scaling to indicate the percentage of zero proportion, by default decide by data automatically
## fontsize.*: font size options for legend, x-, y-axis respectively
## fontangle.*: the x-, y-axis label font align angels, 0 -horizontal align, 90-vertical align,

makeDotplot <- function(input, features = NULL, topN = NULL, cols = c("lightgrey", "blue"),
                        col.min = -2.5, col.max = 2.5, dot.min = 0.01, dot.scale = 6,
                        scale = TRUE, scale.by = "size", scale.min = NA, scale.max = NA,
                        fontsize.legend = 20, fontsize.x = 20, fontsize.y = 20, fontangle.x = 90, fontangle.y = 0) {
  require(ggplot2)
  require(cowplot)
  require(Seurat)

  scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

  input$Cluster_zero_proportion <- 1 - input$Cluster_zero_proportion
  if (is.null(features) & !is.null(topN)) {
    topN <- min(topN, ncol(input$Cluster_mean))
    data.plot.prep <- lapply(input, function(x) x[,1:topN])
    features       <- colnames(input$Cluster_mean)[1:topN]
  } else if (!is.null(features) & is.null(topN)) {
    data.plot.prep <- lapply(input, function(x) x[, match(features,colnames(x))])
  } else {
    stop("please provide either option 'features' or 'topN'.")
  }

  data.plot.mean           <- reshape2::melt(data.plot.prep$Cluster_mean)
  data.plot.zeroProportion <- reshape2::melt(data.plot.prep$Cluster_zero_proportion)

  if (any(data.plot.mean$Var1 != data.plot.zeroProportion$Var1) | any(data.plot.mean$Var2 != data.plot.zeroProportion$Var2)) stop('input data.plot error!')

  data.plot.prep2   <- data.frame(avg.exp = data.plot.mean$value,
                                  pct.exp = data.plot.zeroProportion$value,
                                  features.plot = data.plot.mean$Var2,
                                  id = data.plot.mean$Var1)


  avg.exp.scaled  <- sapply(X = unique(x = data.plot.prep2$features.plot),
                            FUN = function(x) {
                              data.use <- data.plot.prep2[data.plot.prep2$features.plot == x, "avg.exp"]
                              if (scale) {
                                data.use <- scale(x = data.use)
                                data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                              } else {
                                data.use <- log1p(x = data.use)
                              }
                              return(data.use)
                            })

  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot.prep2$avg.exp.scaled <- avg.exp.scaled
  data.plot.prep2$features.plot <- factor(x = data.plot.prep2$features.plot, levels = features)
  data.plot.prep2$pct.exp[data.plot.prep2$pct.exp < dot.min] <- NA
  data.plot.prep2$pct.exp <- data.plot.prep2$pct.exp * 100
  if (!is.na(x = scale.min)) {
    data.plot.prep2[data.plot.prep2$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot.prep2[data.plot.prep2$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  color.by = 'avg.exp.scaled'
  if (fontangle.y < 90) {
    hjustVal = 1 ##right align
    vjustVal.y = 0.5
  } else if (fontangle.y > 90 & fontangle.y < 180) {
    hjustVal = 0 ##left align
    vjustVal.y = 1
  } else {
    hjustVal = 0 ##left align
    vjustVal.y = 0.5
  }

  plot <- ggplot(data = data.plot.prep2, mapping = aes_string(x = "features.plot", y = "id"))
  plot <- plot + geom_point(mapping = aes_string(size = "pct.exp", color = color.by))
  plot <- plot + scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max))
  plot <- plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  plot <- plot + guides(size = guide_legend(title = "Percent Expressed"))
  plot <- plot + labs(x = "", y = "HIPPO Cluster" ) + theme_cowplot()
  if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  plot <- plot + theme(legend.text = element_text(color = "black", size = fontsize.legend),
                       axis.text.x = element_text(color = "black", size = fontsize.x, angle = fontangle.x, vjust = 0.5),
                       axis.text.y = element_text(color = "black", size = fontsize.y, angle = fontangle.y, hjust = hjustVal, vjust = vjustVal.y))
  # , axis.title = element_blank() )

  return(plot)
}

