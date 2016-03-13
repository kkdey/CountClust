#' @title extractTopFeatures- extracting top driving genes driving GoM clusters
#'
#' @param theta The cluster probability distribution/theta matrix obtained from the GoM model fitting (it is a G x K matrix where G is number of features, K number of topics)
#' @param top_features  The number of top features per cluster that drives away that cluster from others. Default value is 10
#' @param method  The underlying model assumed for KL divergence measurement. Two choices considered- "bernoulli" and "poisson"
#' @param options if "min", for each cluster k, we select features that maximize the minimum KL divergence of
#'        cluster k against all other clusters for each feature. If "max", we select features  that maximize the
#'        maximum KL divergence of cluster k against all other clusters for each feature.
#'
#' @return A matrix (K x top_features) which tabulates in k th  row the top feature indices driving the cluster k.
#'
#' @description This function uses relative expression profile of the GoM clusters for each feature, and applies a KL divergence mechanism to obtain
#' a list of top features that drive the clusters.
#'
#' @examples
#' data("MouseDeng2014.FitGoM")
#' theta_mat <- MouseDeng2014.FitGoM$clust_6$theta;
#' top_features <- ExtractTopFeatures(theta_mat, top_features=100,
#'                                   method="poisson", options="min");
#'
#' @export

ExtractTopFeatures <- function(theta, top_features=10, method=c("poisson","bernoulli"), options=c("min", "max"))
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
    KL_score <- lapply(1:dim(theta)[2], function(n) {
      out <- t(apply(theta, 1, function(x){
        y=x[n] *log(x[n]/x) + (1 - x[n])*log((1-x[n])/(1-x));
        return(y)
      }));
      return(out)
    })
  }

  indices_mat=matrix(0,dim(theta)[2],top_features);

  if(dim(theta)[2]==2){
    for(k in 1:dim(theta)[2])
    {
      temp_mat <- KL_score[[k]][,-k];
      if(options=="min"){vec <- apply(as.matrix(temp_mat), 1, function(x) min(x))}
      if(options=="max"){vec <- apply(as.matrix(temp_mat), 1, function(x) max(x))}
      #vec <- temp_mat;
      ordered_kl <- order(vec, decreasing = TRUE);
      counter <- 1
      flag <- counter
      while(flag <= top_features)
      {
        if(max(theta[ordered_kl[counter],])==k){
          indices_mat[k, flag] <- ordered_kl[counter];
          flag <- flag + 1;
          counter <- counter + 1;}
        else{
          counter <- counter + 1;
        }
      }
    }

  } else{
    for(k in 1:dim(theta)[2])
    {
      temp_mat <- KL_score[[k]][,-k];
      if(options=="min"){vec <- apply(temp_mat, 1, function(x) min(x))}
      if(options=="max"){vec <- apply(temp_mat, 1, function(x) max(x))}

      ordered_kl <- order(vec, decreasing = TRUE);
      counter <- 1
      flag <- counter
      while(flag <= top_features)
      {
        if(which.max(theta[ordered_kl[counter],])==k){
          indices_mat[k, flag] <- ordered_kl[counter];
          flag <- flag + 1;
          counter <- counter + 1;}else{
          counter <- counter + 1;
        }
      }
    }
  }

  return(indices_mat);
}
