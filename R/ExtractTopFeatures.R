#' @title The CountsClust function- a function for clustering and visualization of counts!
#'
#' @param theta The theta matrix obtained from the topic model fitting (it is a G x K matrix where G is number of features, K number of topics)
#' @param top_features  The number of top features per cluster that drives away that cluster from others. Default value is 10
#' @param method  The underlying model assumed for KL divergence measurement. Two choices considered- "bernoulli" and "poisson"
#'
#' @description This function uses the relativ expression profile of the clusters across the features, and applies a KL divergence mechanism to obtain
#' a list of features that are significant in driving the clusters. Often annotating these features may give us an idea as to whether there is a common
#' property among the features for each cluster.
#'
#' @author Kushal K Dey, Matthew Stephens
#'
#' @export
#'
#'

ExtractTopFeatures <- function(theta, top_features=10, method=c("poisson","bernoulli"))
{
  if(method=="poisson") {
  KL_score <- lapply(1:dim(theta)[2], function(n) {
    out <- t(apply(theta, 1, function(x){
      y=x[n] *log(x[n]/x) + x - x[n];
      return(y)
    }));
    return(out)
  })
  }

  if(method=="bernoulli"){
    KL_score <- lapply(1:K, function(n) {
      out <- t(apply(theta, 1, function(x){
        y=x[n] *log(x[n]/x) + (1 - x[n])*log((1-x[n])/(1-x));
        return(y)
      }));
      return(out)
    })
  }

  if(dim(theta)[2]==2){
    indices_mat=matrix(0,dim(theta)[2],top_features);

    for(k in 1:dim(theta)[2])
    {
      temp_mat <- KL_score[[k]][,-k];
      vec <- as.vector(temp_mat);
      indices_mat[k,] = order(vec, decreasing = TRUE)[1:top_features]
    }

  } else{
    indices_mat=matrix(0,dim(theta)[2],top_features);

    for(k in 1:dim(theta)[2])
    {
      temp_mat <- KL_score[[k]][,-k];
      vec <- apply(temp_mat, 1, function(x) max(x))
      indices_mat[k,] = order(vec, decreasing = TRUE)[1:top_features]
    }
  }

  return(indices_mat);
}
