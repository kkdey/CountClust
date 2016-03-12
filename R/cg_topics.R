#' @title The center of gravity of the clusters
#'
#' @param theta The theta matrix obtained from the topic model fitting (it is a G x K matrix where G is number of features, K number of topics)
#' @param feature.comp the numeric matrix (G x L) comprising of values for each feature g and feature metadata l.

#' @description This function computes the center of gravity for each cluster by taking weighted mean
#'  of each component of features where the weights are determined from the theta matrix of the topic model fit.
#'
#' @return  Returns a matrix of cluster centers of gravity for the L feature metadata.
#' @author Kushal K Dey, Matthew Stephens
#'
#' @examples
#'
#' N=360;
#' M=560;
#' lat <- rep(1:N, M);
#' long <- rep(1:M, each=N)
#' comp <- cbind(lat, long);
#' Topic_clus <- readRDS("../data/topic_clus_2_maps_no_abundance.rds")
#' center_gravity <- cg_topics(Topic_clus$theta, comp);
#'
#' @export


cg_topics <- function(theta, feature.comp)
{
  gravity_center_list <-  lapply(1:dim(theta)[2], function(l) return(t(theta[,l]) %*% feature.comp))
  return(gravity_center_list)
}
