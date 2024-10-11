#' Bar plot of cluster center values
#'
#' \code{PlotClusterCenters} is a function that plots the centers of continuous
#' features by cluster. This function is meant to be used with
#' \code{\link[ResidualKMeans]{residkm}} and is part of its default output to
#' cluster the continuous residuals of mixed type (i.e., continuous and
#' categorical) features regressed on some grouping variable, such as study
#' cohort. However, this function can be applied to any clustering results with
#' continuous features.
#'
#' @param centers This is a data frame or matrix of cluster centers as typically
#' output by the \code{\link[stats]{kmeans}} function. The format is one row for
#' each cluster with the cluster names being saved as \code{rownames} and a
#' column for each feature with the appropriate \code{colnames} naming those
#' features. This defaults to NULL, in which case one could instead provide a
#' data frame \code{Data} containing only the features used in a cluster
#' analysis and a vector of cluster centers \code{clustvec}.
#' @param Data This is an optional data frame containing continuous features.
#' This should only include columns of features that were part of the
#' clustering. This can be provided along with the vector of cluster membership
#' \code{clustvec} as an alternative to providing \code{centers}. This defaults
#' to \code{NULL}, meaning that \code{centers} must be provided instead.
#' @param clustvec This is an optional vector of cluster membership. This should
#' have the same length as \code{nrow(Data)}. This defaults to \code{NULL},
#' meaning that \code{centers} must be provided instead. If \code{clustvec} is
#' provided with either \code{centers} or \code{Data}, cluster-specific counts
#' are added to the labels of clusters (e.g., "Cluster 1 (N=200)").
#' @param altfeatnames A character vector of alternative feature names to
#' replace the feature \code{colnames} in \code{centers} or \code{Data}.
#' @param featgroups An optional character vector of group names for the
#' features of the same length as the features. This grouping is shown in the
#' output plot. This defaults to \code{NULL}, meaning that no grouping of
#' features is applied.
#' @param scalefeatures A logical value determining whether feature values
#' should be centered and scaled prior to determining cluster centers if
#' \code{Data} is provided but \code{centers} is not. This should be done if
#' the clustering done to obtain \code{clustvec} was done on centered and
#' scaled data, but \code{Data} has features on the original scale instead.
#' This defaults to \code{FALSE}, meaning that centers are calculated for the
#' features in \code{Data} on their provided scale.
#' @param clusterorder A vector of the preferred order of the cluster names. This
#' defaults to \code{NULL}, meaning that the order of clusters is taken to be
#' the row names of \code{centers} if that's provided or the order of unique
#' values of \code{clustvec} if \code{centers} is not provided.
#' @param fillpalette If \code{featgroups} is not provided, the bar plot is
#' filled based on the value of the feature coordinates at the cluster centers.
#' This is done using \code{\link[ggplot2]{scale_fill_viridis_c}}, and
#' \code{fillpalette} determines which viridis palette to use. This defaults to
#' \code{"viridis"}.
#' @param palettealpha This is a value between 0 and 1 determining the
#' transparency of the barplot. This defaults to 1, meaning it is entirely
#' opaque.
#' @param groupdiscpalette If feature groups are defined, the heatmap will be
#' filled based on a discrete palette and \code{fillpalette} will be ignored.
#' This can be the name of a palette available through the package
#' \code{RColorBrewer} (see \code{\link[RColorBrewer]{display.brewer.all}}) or
#' a custom palette function. This defaults to \code{"Set1"}.
#'
#' @returns This returns a vertical bar plot of features faceted vertically
#' by cluster and, if provided, groupings of features. This is a plot of class
#' \code{"ggplot"}.
#'
#' @import ggplot2 cowplot reshape2
#' @export PlotClusterCenters
#'
#' @examples
#' hubhighsep <- GRCsim(nDisClust = 4, nCohortClust = 4, ncontvars = 8,
#' ncatvars = 7, DisSepVal = 0.6, CohortSepVal = 0.2, catq = 0.8, nSignal = 15,
#' nNoise = 0, nOutliers = 0, nrep = 1, DisClustSizes = c(600, 200, 100, 100),
#' CohortClustSizes = c(400, 200, 300, 100), CDS = F,
#' DisClustseed = 200, CohortClustseed = 300, CDSrho = 0.7)$DataList
#'
#' clusters <- residkm(hubhighsep[[1]][, c(paste0("x", 1:15), "Cohort")],
#' groupcolumn = "Cohort", altfeatnames = paste0("Feature ", 1:15))
#'
#' # residkm automatically creates a center plot as part of its output
#' clusters$CenterPlot
#'
#' # We can revise the plot to group outcomes and add a custom palette
#'
#' custompalette <- function(n) {
#' cols <- c("orange2", "deeppink3", "steelblue3", "forestgreen", "grey50",
#' "medium turquoise")
#' cols[1:n]
#' }
#'
#' PlotClusterCenters(centers = clusters$Kmeans$centers,
#' clustvec = clusters$Kmeans$cluster, altfeatnames = paste0("Outcome ", 1:15),
#' featgroups = rep(c("Resp", "Psych", "History"), c(3, 5, 7)),
#' groupdiscpalette = custompalette)
#'
PlotClusterCenters <- function(centers = NULL, Data = NULL, clustvec = NULL,
                               altfeatnames = NULL, featgroups = NULL,
                               scalefeatures = F, clusterorder = NULL,
                               fillpalette = "viridis", palettealpha = 1,
                               groupdiscpalette = "Set1"){
  if(is.null(centers)){
    if(is.null(Data)|is.null(clustvec)){
      stop(paste0("If centers is NULL, you must provide data and clustvec."))
    }
    kmc <- as.data.frame(getcenters(Data, clustvec, scalefeatures))
    if(!is.null(clusterorder)) kmc <- kmc[clusterorder, ]
  } else {
    if(!is.null(Data)) warning(paste0(
      "centers was provided, so Data is being ignored."))
    if(!is.data.frame(centers)) centers <- as.data.frame(centers)
    kmc <- centers
    allintnames <- all(!is.na(as.integer(rownames(kmc))))
    if(allintnames){
      rownames(kmc) <- paste0("Cluster ", rownames(kmc))
    }
  }
  if(!is.null(altfeatnames)){
    featnames <- names(kmc) <- altfeatnames
  } else featnames <- names(kmc)
  kmc$Cluster <- factor(rownames(kmc), levels = rownames(kmc))
  if(!is.null(clustvec)){
    clusternobs <- summary(as.factor(clustvec))
    clusternobs <- clusternobs[match(rownames(centers), names(clusternobs))]
    levels(kmc$Cluster) <- paste0(levels(kmc$Cluster), " (N=", clusternobs, ")")
  }
  kmc <- reshape2::melt(kmc, id.vars = "Cluster")
  kmc$variable <- factor(kmc$variable, levels = rev(featnames))
  if(!is.null(featgroups)){
    ftd <- data.frame(Types = featgroups, Feats = featnames)
    kmc$Type <- factor(ftd[match(as.character(kmc$variable), ftd$Feats), "Types"],
                       levels = unique(featgroups))
    centerplot <- ggplot(kmc, aes(y = variable, x = value)) + theme_bw() +
      geom_bar(aes(fill = Type), stat = "identity", alpha = palettealpha) +
      geom_vline(xintercept = 0) + geom_text(aes(label = round(value, 2))) +
      ggh4x::facet_nested(Cluster + Type ~ ., scales = "free_y",
                          space = "free_y") +
      xlab("Center Z-score")
    if(!is.function(groupdiscpalette)){
      centerplot <- centerplot + scale_fill_brewer(palette = groupdiscpalette)
    } else {
      centerplot <- centerplot + scale_fill_gradientn(
        groupdiscpalette(length(unique(featgroups))))
    }
    centerplot <- centerplot +
      theme(legend.position = "bottom", axis.title.y = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = 'bold'))
  } else {
    centerplot <- ggplot(kmc, aes(y = variable, x = value)) + theme_bw() +
      geom_bar(aes(fill = value), stat = "identity", alpha = palettealpha) +
      geom_vline(xintercept = 0) + geom_text(aes(label = round(value, 2))) +
      scale_fill_viridis_c(option = fillpalette) +
      facet_grid(Cluster ~ .) + xlab("Center Z-score") +
      theme(legend.position = "none", axis.title.y = element_blank())
  }
  return(centerplot)
}


