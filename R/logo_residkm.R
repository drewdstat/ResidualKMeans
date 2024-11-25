#' Leave-One-Group-Out (LOGO) Sensitivity Analysis for Residual K-Means
#'
#' \code{logo_residkm} is a function to perform a sensitivity analysis for the
#' residual k-means algorithm. Since residual k-means was designed to perform
#' clustering independent of a grouping variable (e.g., cohort, location, group,
#' etc.), a useful sensitivity analysis is to evaluate whether one particular
#' group is driving the clustering despite the algorithm focusing on clustering
#' independent of that grouping variable. For example, if our residual k-means
#' analysis is focused on clustering variables independent of study cohort in a
#' multi-cohort study, we would want to make sure the clusters identified in the
#' full sample aren't primarily reliant on the inclusion of one particular cohort.
#' The leave-one-group-out (LOGO) analysis omits one group at a time and repeats
#' the entire residual k-means procedure. Then a list of summary results is
#' returned along with a measure of cluster stability for each LOGO iteration so
#' that the consistency and stability of the clustering results can be evaluated.
#'
#' @param data This is a data frame including only of the variables to be
#' clustered and a column with group identities for the residual k-means
#' algorithm.
#' @param groupvar This is a character value representing the column name of the
#' grouping variable in the data frame \code{data}. This defaults to
#' \code{"Cohort"}.
#' @param includefull This is a logical (boolean) value defining whether to
#' also run residual k-means on the full sample first and to include those
#' results in the function output, including in the plots if
#' \code{returnplots = TRUE}. Otherwise, only the LOGO results without the full
#' sample results will be returned. This defaults to \code{TRUE}.
#' @param krange This is a single integer or a vector of integers representing
#' k values (i.e., numbers of cluster partitions). If \code{includefull = TRUE}
#' and \code{length(krange) > 1}, an optimal number of clusters will be selected
#' by the \code{\link[NbClust]{NbClust}} algorithm as part of that full-sample
#' residual k-means run, and then that number of clusters will be applied to all
#' subsequent LOGO analyses. If \code{includefull = FALSE} and
#' \code{length(krange) > 1}, only the first value in that vector will be used
#' in all LOGO analyses, and otherwise if \code{length(krange) = 1}, then that
#' number of clusters will be applied to all LOGO iterations. This defaults to
#' \code{2:10} since \code{includefull} defaults to \code{TRUE}.
#' @param jaccardboot This is a logical value defining whether the LOGO analysis
#' (and also the full-sample analysis if \code{includefull = TRUE}) also
#' includes a bootstrapped Jaccard index analysis as implemented by
#' \code{\link[fpc]{clusterboot}}. The sample is bootstrapped many times, and
#' each time the similarity between the bootstrapped k-means clustering analysis
#' and the initial k-means clustering is evaluated with the Jaccard index. Given
#' the bootstrapping, this can take some time to run. This defaults to
#' \code{TRUE}.
#' @param jacb If \code{jaccardboot = TRUE}, this value is the number of
#' bootstrap iterations to implement in the bootstrap Jaccard index analysis.
#' This defaults to \code{1000}.
#' @param exclgroups This is an optional character vector to define which groups
#' within the \code{groupvar} variable to include in the LOGO analysis. For
#' example, if one has 5 groups in \code{groupvar} named Group 1, etc. and one
#' only wants to do a short LOGO analysis with just the first 3 groups, this
#' would be set to \code{exclgroups = paste0("Group ", 1:3)}. This value
#' defaults to \code{NULL}, meaning that all groups within \code{groupvar} will
#' be included in the LOGO analysis.
#' @param exclgroupnames This is an optional character vector of names for each
#' group within the \code{groupvar} variable. This must be of the same length as
#' the number of groups within \code{groupvar} or the length of
#' \code{exclgroups} if that is provided. This defaults to \code{NULL}.
#' @param returnplots This is a logical value defining whether to also output
#' three plots summarizing the LOGO results: another plotting the cluster
#' centers, another showing the variable summaries, and one visualizing the
#' Jaccard indices (if \code{jaccardboot = TRUE}). If \code{FALSE}, only the
#' list of LOGO \code{residkm} results will be returned. This defaults to
#' \code{TRUE}.
#' @param groupname This is an optional character value to use as the title
#' of the group variable in the output plots if \code{returnplots = TRUE}. If
#' this value is \code{NULL} (the default), then the group variable title will
#' be \code{"Group"}.
#' @param origcenters This is an optional matrix wherein each row represents
#' a cluster, each column represents a variable included in a residual k-means
#' clustering analysis, and the values in that matrix are the cluster center
#' residual z-scores values for those variables. This cluster center matrix is
#' output by every \code{residkm} run and stored as \code{$Kmeans$center} within
#' the results object. This matrix can be provided so that the cluster identity
#' assigned to every LOGO iteration can be matched to the original cluster
#' identities by finding the minimum Euclidean distance between each LOGO
#' cluster center and the original cluster centers. Otherwise, the \code{residkm}
#' output will assign Cluster 1, Cluster 2, etc. based on greatest to smallest
#' sample size within each cluster, and this may not always line up with the
#' original Cluster 1, Cluster 2, etc. It is therefore recommended to provide
#' a matrix of "reference" cluster centers for \code{origcenters} if
#' \code{includefull = FALSE} to ensure all LOGO clusters have the same
#' interpretation and identity. If \code{includefull = TRUE}, then the cluster
#' centers from that full sample residual k-means analysis will be used as the
#' "reference" cluster centers and \code{origcenters} will be ignored. This
#' defaults to \code{NULL}.
#' @param altfeatnames This is an optional character vector of names for each
#' variable (i.e., "feature") in the clustering analysis to be used in plots
#' that show variable-specific cluster centers or summary statistics. This
#' defaults to \code{NULL}.
#' @param origclusternames This is an optional character vector of names for
#' each of the original clusters to be used in plots. If this is set to
#' \code{NULL} (the default), these clusters will just by called Cluster 1, etc.
#' @param variable_pal This is an optional discrete color palette to be applied
#' to the LOGO center plot, which color codes center identities over LOGO
#' iterations by variable type. Therefore, if a palette is provided, it must
#' have enough colors to represent each of the included variables. This defaults
#' to \code{NULL}, and so the default ggplot discrete palette will be used.
#' @param ... Additional arguments to pass to \code{residkm}.
#'
#' @return \code{logo_residkm} returns a list that includes the following:
#' \item{Results}{A list of objects of class \code{"residkm"} for all the
#' LOGO results. If \code{jaccardboot = TRUE}, each object will include an
#' additional list object "Stability" for the \code{fpc::clusterboot} output.}
#' \item{CenterData}{If \code{returnplots = TRUE}, this data frame of cluster
#' center values for all clustered variables over LOGO iterations is output.}
#' \item{CenterPlot}{If \code{returnplots = TRUE}, this plot of cluster center
#' values for all clustered variables over LOGO iterations is output.}
#' \item{SummaryPlot}{If \code{returnplots = TRUE}, this multiplot of summary
#' statistics for each LOGO iteration is output. Note that this is drawn from
#' the summary plot output of \code{residkm} for each LOGO iteration, and
#' therefore the clusters are named based on decreasing sample size and
#' may change names between LOGO iterations (e.g., Cluster 3 is labelled as
#' Cluster 2).}
#' \item{JaccardData}{If \code{returnplots = TRUE} and
#' \code{jaccardboot = TRUE}, this data frame of Jaccard indices over LOGO
#' iterations is output.}
#' \item{JaccardPlot}{If \code{returnplots = TRUE} and
#' \code{jaccardboot = TRUE}, this plot of Jaccard indices over LOGO iterations
#' is output.}
#'
#' @import ggplot2 NbClust missMDA fpc ggh4x gplots reshape2
#' @export logo_residkm
#'
#' @examples
#' # Creating simulated data with GRCsim and analyze it with residkm
#'
#' #Single data frame of "Hub condition, medium separation" from Day et al. 2024
#'
#' hubhighsep <- GRCsim(nDisClust = 4, nCohortClust = 4, ncontvars = 8,
#' ncatvars = 7, DisSepVal = 0.6, CohortSepVal = 0.2, catq = 0.8, nSignal = 15,
#' nNoise = 0, nOutliers = 0, nrep = 1, DisClustSizes = c(600, 200, 100, 100),
#' CohortClustSizes = c(400, 200, 300, 100), CDS = F,
#' DisClustseed = 200, CohortClustseed = 300, CDSrho = 0.7)$DataList
#'
#' logo_ex1 <- logo_residkm(hubhighsep[[1]][, c(paste0("x", 1:15), "Cohort")],
#' krange = 4, jacb = 3)
#' # Note that we recommend the default jacb of 1000. jacb is low here for
#' # a speedy run through of the example, but it should be much higher.
#'
#' #Center Plot
#' logo_ex1$CenterPlot
#'
#' #Summary Plots
#' logo_ex1$MeanPlot
#'
#' #Jaccard Plot
#' logo_ex1$JaccardPlot
#'
logo_residkm <- function(data = NULL, groupvar = "Cohort", includefull = T,
                         krange = 2:10, jaccardboot = T, jacb = 1000,
                         exclgroups = NULL, exclgroupnames = NULL,
                         returnplots = T, groupname = NULL, origcenters = NULL,
                         altfeatnames = NULL, origclusternames = NULL,
                         variable_pal = NULL, ...){
  #Get group levels to omit
  if(!is.null(exclgroups) & !is.null(exclgroupnames) &
     length(exclgroups) != length(exclgroupnames)){
    stop(paste0("If exclgroups and exclgroupnames are both provided, they must",
                " have equal lengths."))
  }
  if(!is.factor(data[, groupvar])) data[, groupvar] <- as.factor(data[, groupvar])
  if(is.null(exclgroups)){
    grpquant <- summary(data[, groupvar])
    grpquant <- grpquant[order(-grpquant)]
    grplvls <- names(grpquant)
  } else {
    grplvls <- exclgroups
  }
  #Run full residual k-means if includefull
  if(includefull){
    fullres <- residkm(data, groupvar, krange, ksel = length(krange) > 1, ...)
    if(length(krange) > 1) krange <-
        length(unique(fullres$KChoice$Best.partition))
    if(is.null(origcenters)) origcenters <- fullres$Kmeans$centers
    if(jaccardboot){
      fullres$Stability <-
        fpc::clusterboot(scale(fullres$ResidualData), B = jacb,
                         k = krange, iter.max = 100, nstart = 1000,
                         clustermethod = fpc::kmeansCBI, count = F,
                         bootmethod = "boot")
    }
    if(is.null(origclusternames)) origclusternames <-
      paste0("Cluster ", 1:length(unique(fullres$Kmeans$cluster)))
  }
  if(!includefull & length(krange) > 1){
    krange <- krange[1]
    warning(paste0("krange can only be of length 1 if a range of numbers of ",
                   "clusters is not being supplied to an initial full sample",
                   " run of residual k-means. The number of clusters for the ",
                   "LOGO analysis is being set to the first number in krange: ",
                   krange, "."))
  }

  #LOGO function to apply
  runlogo <- function(grp, ...){
    outres <- residkm(data[which(data[, groupvar] != grp), ],
                             groupvar, krange = krange, ksel = F, ...)
    if(jaccardboot){
      outres$Stability <-
        fpc::clusterboot(scale(outres$ResidualData), B = jacb,
                         k = krange, iter.max = 100, nstart = 1000,
                         clustermethod = fpc::kmeansCBI, count = F,
                         bootmethod = "boot")
    }
    return(outres)
  } # end runlogo
  #run runlogo
  pbapply::pboptions(type = "timer")
  logores <- pbapply::pblapply(grplvls, runlogo)
  names(logores) <- grplvls
  if(includefull){
    resnm <- names(logores)
    logores$Full <- fullres
    logores <- logores[c("Full", resnm)]
  } else if(is.null(origcenters)){
    origcenters <- logores[[1]]$Kmeans$centers
  }

  if(returnplots){
    #renaming
    if(!is.null(groupname)) grpnm <- groupname else grpnm <- "Group"
    if(!is.null(exclgroupnames)) names(logores) <- exclgroupnames
    if(is.null(origclusternames)) origclusternames <- paste0("Cluster ",
                                                             1:krange)

    #Stability (Jaccard) Plots
    if(jaccardboot){
      jacmat <- sapply(logores, function(x) x[["Stability"]]$bootmean[
        matchclustercenter(x$Kmeans$centers, origcenters)])
      jacmat <- as.data.frame(jacmat)
      jacmat <- suppressMessages(reshape2::melt(jacmat))
      names(jacmat) <- c("Omitted", "Jaccard")
      jacmat$Cluster <- as.factor(rep(origclusternames,
                                  nrow(jacmat)/length(origclusternames)))
      jacmat$Omitted <- factor(jacmat$Omitted, levels = names(logores))
      jacfig <- ggplot(jacmat, aes(x = Omitted, y = Jaccard, color = Cluster)) +
        theme_bw() + xlab(paste0("Omitted ", grpnm)) +
        ylab("Bootstrap Jaccard Index") +
        geom_line(aes(group = Cluster), linewidth = 1.2) +
        geom_point(show.legend = F, size = 2) +
        geom_text(aes(label = round(Jaccard, 3)), size = 3,
                  vjust = -0.75, color = "black") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    #Center Plots
    logocents <- as.data.frame(do.call(
      "rbind", lapply(logores, function(x) x[["Kmeans"]]$centers)))
    tmp <- suppressMessages(reshape2::melt(sapply(
      logores, function(x) matchclustercenter(x[["Kmeans"]]$centers,
                                              origcenters))))
    logocents$Omitted <- tmp$Var2
    logocents$Cluster <- as.factor(tmp$value)
    rownames(logocents) <- NULL
    logocents <- suppressMessages(reshape2::melt(logocents, id.vars = c(
      "Omitted", "Cluster"), value.name = "value", variable.name = "Variable"))
    if(is.null(altfeatnames)) varnames <-
      colnames(logores[[1]][["Kmeans"]]$centers) else varnames <- altfeatnames
    levels(logocents$Variable) <- varnames
    logocents$Omitted <- factor(logocents$Omitted, levels = names(logores))
    levels(logocents$Cluster) <- origclusternames
    centfig <- ggplot(logocents, aes(x = Omitted, y = value, color = Variable)) +
      theme_bw() + xlab("") + ylab("Center Residual Z-Score") +
      geom_line(aes(group = Variable), linewidth = 1.3) +
      geom_point(show.legend = F, size = 2.5) +
      facet_grid(~Cluster) +
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            legend.position = "bottom", legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = 'bold'),
            legend.key.width = unit(0.5, 'in'),
            legend.key.height = unit(0.2, 'in'),
            strip.text = element_text(size = 12))
    if(!is.null(variable_pal)) centfig <- centfig +
      scale_color_manual(values = variable_pal)

    #Variable Summary Plots
    meanplotlist <- list()
    for(i in names(logores)){
      meanplotlist[[i]] <- logores[[i]][["SummaryPlot"]]
    }; rm(i)
    meanfig <- cowplot::plot_grid(plotlist = meanplotlist, ncol = 3,
                                  labels = names(logores), hjust = 0)
    reslist <- list(Results = logores, CenterData = logocents,
                    CenterPlot = centfig, SummaryPlot = meanfig)
    if(jaccardboot){
      reslist$JaccardData <- jacmat
      reslist$JaccardPlot <- jacfig
    }
  } else reslist <- list(Results = logores)
  return(reslist)
}
