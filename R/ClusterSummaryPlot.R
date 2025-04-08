#' Visualize summaries of variables by cluster
#'
#' \code{ClusterSummaryPlot} is a function to create a heatmap plot that
#' summarizes continuous and categorical variable data between clusters. This
#' function is meant to be used with the residual k-means method. The plot
#' consists of a heatmap separated vertically between the continuous variables
#' on top and the categorical variables below.
#'
#' @param data This is a data frame containing columns for the continuous
#' variables (\code{contvars}) and categorical variables (\code{catvars}) to be
#' summarized by this function.
#' @param clustvec This is a character vector of cluster membership. This
#' should be of the same length as the number of rows (i.e., observations) in
#' \code{data}.
#' @param contvars This is a character vector of column names for continuous
#' variables in the data frame \code{data}.
#' @param catvars This is a character vector of column names for categorical
#' variables in the data frame \code{data}. These should be binary variables of
#' class \code{character}, \code{factor}, or \code{integer} with positive values
#' being equal to one of the following values:
#' \code{c(1, "y", "Y", "yes", "Yes", "YES")}. For categorical variables with >2
#' unique values, these shoudl be dummy encoded into binary variables.
#' @param contgroups An optional character vector of group names for the
#' \code{contvars} of the same length as \code{contvars}. This grouping is
#' shown in the output plot. If \code{contgroups} is provided, \code{catgroups}
#' must also be provided. This defaults to \code{NULL}, meaning that no
#' grouping of features is applied.
#' @param catgroups An optional character vector of group names for the
#' \code{catvars} of the same length as \code{catvars}. This grouping is
#' shown in the output plot. If \code{catgroups} is provided, \code{contgroups}
#' must also be provided. This defaults to \code{NULL}, meaning that no
#' grouping of features is applied.
#' @param contwhitecutoff An optional value on either the original scale of the
#' \code{contvars} or on a z-score scale (if \code{scale = TRUE}) above which
#' the font of the values displayed on the heatmap will be black and below which
#' or equal to which the font will be white. This is useful for seeing values
#' against a dark background. This defaults to \code{NULL}, meaning that all
#' text font will be black.
#' @param catwhitecutoff A similar optional value to \code{contwhitecutoff}, but
#' on the scale of 0 - 100 for the percentages of the binary categorical
#' variables.
#' @param contpalette The name of a viridis palette to apply to the continuous
#' variable summary heatmap (see \code{\link[ggplot2]{scale_colour_viridis_d}}
#' for palette names and
#' \href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{
#' the introductory vignette to viridis} to see the palettes). This defaults to
#' \code{"viridis"}.
#' @param catpalette Similar to \code{contpalette} but for the categorical
#' variable summary heatmap. This defaults to \code{"plasma"}.
#' @param palettealpha A value between 0 and 1 to scale the transparency of the
#' heatmap background fill colors. This defaults to 0.7.
#' @param groupdiscpalette If feature groups are defined, the heatmap will be
#' filled based on a discrete palette and both \code{contpalette} and
#' \code{catpalette} will be ignored. This can be the name of a palette
#' available through the package \code{RColorBrewer} (see
#' \code{\link[RColorBrewer]{display.brewer.all}}) or a custom palette function.
#' This defaults to \code{"Set1"}.
#' @param title An optional title to add to the plot. This defaults to
#' \code{NULL}, meaning that no title is applied.
#' @param titlesize The size of the title font if \code{title} is not
#' \code{NULL}. This defaults to \code{12}.
#' @param xlabelsize The font size of the x-axis title, which will be the word
#' "Cluster". This defaults to \code{12}.
#' @param yaxistextsize The font size of the text along the y-axis, which will
#' be the variable/feature names. This defaults to \code{10}.
#' @param xaxistextsize The font size of the x-axis text, which will be the
#' cluster names. This defaults to \code{10}.
#' @param na.rm This is a logical (boolean) value determining whether to remove
#' missing values from the total sample size (i.e., the denominator) when
#' calculating the percentages for each binary categorical variable. For
#' example, if there are 2 "yes" responses out of 10 total observations with 4
#' missing values, the percentage for that variable will be 20% if
#' \code{na.rm = FALSE} and 33.3% if \code{na.rm = TRUE}. This defaults to
#' \code{TRUE}.
#'
#' @returns \code{ClusterSummaryPlot} returns a compound plot of class
#' \code{ggplot} consisting of two separate heatmaps for the continuous and
#' categorical variables.
#'
#' @import cowplot ggplot2 reshape2 ggh4x
#' @export ClusterSummaryPlot
#'
#' @examples
#' # Generate simulated data, cluster it, and plot the results
#'
#' #Single data frame of "Hub condition, medium separation" from Day et al. 2024
#'
#' hubhighsep <- GRCsim(nDisClust = 4, nCohortClust = 4, ncontvars = 8,
#' ncatvars = 7, DisSepVal = 0.6, CohortSepVal = 0.2, catq = 0.8, nSignal = 15,
#' nNoise = 0, nOutliers = 0, nrep = 1, DisClustSizes = c(600, 200, 100, 100),
#' CohortClustSizes = c(400, 200, 300, 100), CDS = F, nContCDSrootvars = 2,
#' nCatCDSrootvars = 2, nContCDSgenvars = 3, nCatCDSgenvars = 3,
#' DisClustseed = 200, CohortClustseed = 300, CDSrho = 0.7)$DataList
#'
#' clusters <- residkm(hubhighsep[[1]][, c(paste0("x", 1:15), "Cohort")],
#' groupcolumn = "Cohort", altfeatnames = paste0("Feature ", 1:15))
#'
#' # residkm already creates one of these plots as part of its output
#' clusters$SummaryPlot
#'
#' # We can revise the plot to group outcomes and add a custom palette
#'
#' custompalette <- function(n) {
#' cols <- c("orange2", "deeppink3", "steelblue3", "forestgreen", "grey50",
#' "medium turquoise")
#' cols[1:n]
#' }
#'
#' ClusterSummaryPlot(hubhighsep[[1]][, c(paste0("x", 1:15))],
#' clusters$Kmeans$cluster, paste0("x", 1:8), paste0("x", 9:15),
#' contgroups = rep(c("Resp", "Psych"), c(3, 5)),
#' catgroups = rep("History", 7), groupdiscpalette = custompalette)
#'
ClusterSummaryPlot <- function(data, clustvec, contvars, catvars, contgroups=NULL,
                             catgroups= NULL, contwhitecutoff = NULL,
                             catwhitecutoff = NULL, scale = T,
                             contpalette = "viridis", catpalette = "plasma",
                             palettealpha = 0.7, groupdiscpalette = "Set1",
                             title = NULL, titlesize = 12, xlabelsize = 12,
                             yaxistextsize = 10, xaxistextsize = 10, na.rm = T){
  posresponses <- c("Y", 1, "y", "yes", "Yes", "YES")
  checkpositivesymbols <- function(x) any(
    grepl(paste(posresponses, collapse="|"), x))
  checkpositives <- sapply(data[, catvars], checkpositivesymbols)
  if(any(!checkpositives)){
    stop(paste0("It is unclear how positive categorical responses are coded \n",
                "in the provided data frame. Please make sure the categorical \n",
                "columns are binary with positive responses being one of the \n",
                "following values: 1, 'Y', 'y', 'yes', 'Yes', or 'YES'. \n",
                "Categorical variables with >2 unique values should be \n",
                "converted into dummy-coded binary variables."))
  }
  grpbool1 <- !is.null(contgroups) & is.null(catgroups)
  grpbool2 <- is.null(contgroups) & !is.null(catgroups)
  if(grpbool1|grpbool2) stop(paste0("contgroups and catgroups should either",
                                    " both be NULL or both be provided."))
  data$Cluster <- as.factor(clustvec)
  countsumm <- summary(data$Cluster)
  clustcountdf <- data.frame(Cluster = names(countsumm), Count = countsumm)
  data$Count <- clustcountdf[match(data$Cluster,
                                   clustcountdf$Cluster), "Count"]
  data$Cluster <- as.factor(paste0(as.character(data$Cluster),
                                   " (N=", data$Count, ")"))
  data$Count <- NULL
  datacat_long <- reshape2::melt(data[, which(
    names(data) %in% c("Cluster", catvars))], id.vars = "Cluster")
  if(na.rm){
    catpercfunc <- function(x) (length(which(grepl(
      paste(posresponses, collapse = "|"), x)))/length(x[which(!is.na(x))]))*100
  } else {
    catpercfunc <- function(x) (length(which(grepl(paste(
      posresponses, collapse = "|"), x)))/length(x))*100
  }
  heatdfcat <- aggregate(
    value ~ Cluster + variable, datacat_long,
    catpercfunc, na.action = NULL)
  if(!is.null(catgroups)){
    heatdfcat$Type <- catgroups[match(heatdfcat$variable, catvars)]
    heatdfcat$Type <- factor(heatdfcat$Type, levels = unique(catgroups))
    allkeydf <- data.frame(Type = c(unique(contgroups), unique(catgroups)))
    allkeydf$Type <- factor(allkeydf$Type, levels = allkeydf$Type)
    if(is.function(groupdiscpalette)){
      allkeydf$Fill = groupdiscpalette(nrow(allkeydf))
    } else {
      allkeydf$Fill = RColorBrewer::brewer.pal(nrow(allkeydf), groupdiscpalette)
    }
    allkeydf$Value <- 1:nrow(allkeydf)
    fillcols <- allkeydf$Fill
    allkeyplot <- ggplot(allkeydf, aes(x = 1, y = 1, fill = Type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = fillcols) +
      theme(legend.position = "bottom", legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = 'bold'))
    l1 <- cowplot::get_plot_component(allkeyplot,
                                      'guide-box-bottom', return_all = T)
  }

  datacont_orig <- datacont <- data[, c(contvars)]
  if(any(scale)){
    if(length(scale) == 1){scale <- rep(scale, length(contvars))}
    for(i in 1:length(contvars)){
      if(scale[i]){datacont[, contvars[i]] <-
        as.numeric(scale(datacont[, contvars[i]]))}
    }
  }
  datacont$Cluster <- datacont_orig$Cluster <- data$Cluster
  datacont_long <- reshape2::melt(datacont, id.vars = "Cluster")
  heatdfcont <- aggregate(value ~ Cluster + variable, datacont_long, mean)
  heatdfcont2 <- aggregate(value ~ Cluster + variable, datacont_long, sd)
  heatdfcont$SD <- heatdfcont2$value
  rm(heatdfcont2)
  datacont_orig_long <- reshape2::melt(datacont_orig, id.vars = "Cluster")
  heatdfcont_orig <- aggregate(value ~ Cluster + variable,
                               datacont_orig_long, mean)
  heatdfcont_orig2 <- aggregate(value ~ Cluster + variable,
                                datacont_orig_long, sd)
  heatdfcont_orig$SD <- heatdfcont_orig2$value
  rm(heatdfcont_orig2)
  heatdfcont$CenterLab <- heatdfcont_orig$value
  heatdfcont$SDLab <- heatdfcont_orig$SD

  for(i in 1:nrow(heatdfcont)){
    dcp_c <- c(decimalplaces(heatdfcont[i, "CenterLab"]),
               decimalplaces(heatdfcont[i, "CenterLab"], T))
    dcp_s <- c(decimalplaces(heatdfcont[i, "SDLab"]),
               decimalplaces(heatdfcont[i, "SDLab"], T))
    if(dcp_c[2] > 2) {heatdfcont[i, "CenterLab"] <-
      signif(heatdfcont[i, "CenterLab"], 3)
    } else if(dcp_c[1] > 0){heatdfcont[i, "CenterLab"] <-
      round(heatdfcont[i, "CenterLab"], 1)}
    if(dcp_s[2] > 1) {heatdfcont[i, "SDLab"] <-
      signif(heatdfcont[i, "SDLab"], 2)
    } else if(dcp_s[1] > 0){heatdfcont[i, "SDLab"] <-
      round(heatdfcont[i, "SDLab"], 1)}
  }
  heatdfcont$TextLabel <- paste0(heatdfcont$CenterLab,
                                 " (", heatdfcont$SDLab, ")")

  heatdfcat$variable <- factor(heatdfcat$variable,
                               levels = rev(levels(heatdfcat$variable)))
  heatdfcont$variable <- factor(heatdfcont$variable,
                                levels = rev(levels(heatdfcont$variable)))
  if(!is.null(contwhitecutoff)){
    heatdfcont$TextColor <- ifelse(heatdfcont$value > contwhitecutoff,
                                   "black", "white")
  } else {
    heatdfcont$TextColor <- "black"
  }
  if(!is.null(catwhitecutoff)){
    heatdfcat$TextColor <- ifelse(heatdfcat$value > catwhitecutoff,
                                   "black", "white")
  } else {
    heatdfcat$TextColor <- "black"
  }

  heatdfcont$StripTitle <- "Mean (SD)"
  heatdfcat$StripTitle <- "Percent"

  if(na.rm){
    heatdfcont <- heatdfcont[heatdfcont$Cluster != "NA", ]
    heatdfcat <- heatdfcat[heatdfcat$Cluster != "NA", ]
  }

  if(!is.null(contgroups)){
    heatdfcont$Type <- contgroups[match(heatdfcont$variable, contvars)]
    heatdfcont$Type <- factor(heatdfcont$Type, levels = unique(contgroups))
    heatdfcont$FillCol <- allkeydf[match(heatdfcont$Type, allkeydf$Type), "Fill"]
    g1 <- ggplot(heatdfcont, aes(y = variable, x = Cluster, fill = FillCol)) +
      cowplot::theme_cowplot() + geom_tile(aes(alpha = value)) +
      geom_text(aes(label = TextLabel, color = TextColor), size = 3.5) +
      ggh4x::facet_nested(StripTitle + Type ~ ., space = "free", scales = "free") +
      scale_alpha_continuous(range = c(0.05, palettealpha), guide = "none") +
      scale_fill_identity("Type", labels = levels(heatdfcont$Type)) +
      scale_color_identity(guide = "none") +
      theme(legend.position = "none")
  } else {
    g1 <- ggplot(heatdfcont, aes(y = variable, x = Cluster, fill = value)) +
      facet_grid(StripTitle ~ ., space = "free", scales = "free") +
      geom_tile(alpha = palettealpha) + cowplot::theme_cowplot() +
      geom_text(aes(label = TextLabel, color = TextColor), size = 3.5) +
      scale_fill_viridis_c(option = contpalette) +
      scale_color_identity(guide = "none") + labs(fill = "Z") +
      theme(legend.position = "right",
            legend.key.height = unit(length(
              unique(heatdfcont$variable))/5, 'line'),
            legend.title.align = 0.5)
  }

  g1 <- g1 +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = yaxistextsize))

  if(!is.null(title)){g1 <- g1 + ggtitle(title) +
    theme(plot.title = element_text(size = titlesize, hjust = 0.5))}

  if(!is.null(catgroups)){
    heatdfcat$FillCol <- allkeydf[match(heatdfcat$Type, allkeydf$Type), "Fill"]
    g2 <- ggplot(heatdfcat, aes(y = variable, x = Cluster, fill = FillCol)) +
      ggh4x::facet_nested(StripTitle + Type ~ ., space = "free", scales = "free") +
      geom_tile(aes(alpha = value)) + cowplot::theme_cowplot() +
      geom_text(aes(label = round(value, 1), color = TextColor), size = 4) +
      scale_alpha_continuous(range = c(0.05, palettealpha), guide = "none") +
      scale_fill_identity("Type", labels = levels(heatdfcat$Type)) +
      scale_color_identity() + theme(legend.position = "none")
  } else {
    g2 <- ggplot(heatdfcat, aes(y = variable, x = Cluster, fill = value)) +
      facet_grid(StripTitle ~ ., space = "free", scales = "free") +
      geom_tile(alpha = 0.6) +
      geom_text(aes(label = round(value, 1), color = TextColor), size = 4) +
      scale_fill_viridis_c(option = catpalette) + cowplot::theme_cowplot() +
      scale_color_identity(guide = "none") + labs(fill = "%") +
      theme(legend.position = "right",
            legend.key.height = unit(length(
              unique(heatdfcont$variable))/5, 'line'),
            legend.title.align = 0.5)
  }

  g2 <- g2 +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = xlabelsize),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(size = xaxistextsize),
          axis.text.y = element_text(size = yaxistextsize))

  catlength <- length(unique(heatdfcat$variable))
  contlength <- length(unique(heatdfcont$variable))
  catplotht <- (catlength)/(catlength + contlength)
  contplotht <- contlength/(catlength + contlength)
  fig <- cowplot::plot_grid(g1, g2, nrow = 2, align = "v",
                            rel_heights = c(contplotht, catplotht))
  if(!is.null(catgroups)){
    fig <- cowplot::plot_grid(fig, cowplot::ggdraw(l1), ncol = 1,
                              rel_heights = c(1, 0.1))
  }
  return(fig)
}
