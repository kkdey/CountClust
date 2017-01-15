#' @title Extracting top driving genes of GoM clusters
#'
#' @description This function uses relative gene expression profile of the GoM
#'        clusters and applies a KL-divergence based method to
#'        obtain a list of top features that drive each of the clusters.
#'
#' @param theta \eqn{\boldsymbol{theta}} matrix, the relative gene expression profile of the GoM clusters
#'                (cluster probability distributions)
#'                from the GoM model fitting (a \eqn{G x K} matrix where \eqn{G} is
#'                number of features, \eqn{K} number of topics).
#' @param top_features  The top features in each cluster \eqn{k} that are selected based on the feature's
#'                          ability to distinguish cluster \eqn{k} from cluster \eqn{1, \dots, K}
#'                          for all cluster \eqn{k \ne l}. Default: \eqn{10}.
#' @param method  The underlying model assumed for KL divergence measurement.
#'                  Two choices considered are "bernoulli" and "poisson". Default: poisson.
#' @param options if "min", for each cluster k, we select features that
#'                  maximize the minimum KL divergence of cluster k against
#'                  all other clusters for each feature. If "max", we select
#'                  features  that maximize the maximum KL divergence of cluster
#'                  k against all other clusters for each feature.
#' @param shared if TRUE, then we report genes that can be highly expressed in more than one cluster. Else, we stick to
#'               only those genes that are highest expressed only in a specific cluster.
#'
#' @return A matrix (K x top_features) which tabulates in k-th row the
#'        top feature indices driving the cluster k.
#'
#'
#' @examples
#' data("MouseDeng2014.FitGoM")
#' theta_mat <- MouseDeng2014.FitGoM$clust_6$theta;
#' top_features <- ExtractTopFeatures(theta_mat, top_features=100,
#'                                   method="poisson", options="min");
#'
#' @export

ExtractTopFeatures <- function(theta,
                               top_features = 10,
                               method = c("poisson","bernoulli"),
                               options=c("min", "max"),
                               shared = FALSE)
{
    if (is.null(method)) {
        warning("method is not specified! Default method is Poisson distribution.")
        method <- "poisson"
    }
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
            if(options=="min"){
                vec <- apply(as.matrix(temp_mat), 1, function(x) min(x))}
            if(options=="max"){
                vec <- apply(as.matrix(temp_mat), 1, function(x) max(x))}
            #vec <- temp_mat;
            ordered_kl <- order(vec, decreasing = TRUE);
            counter <- 1
            flag <- counter
            while(flag <= top_features)
            {
                if(!shared){
                if(which.max(theta[ordered_kl[counter],])==k){
                    indices_mat[k, flag] <- ordered_kl[counter];
                    flag <- flag + 1;
                    counter <- counter + 1;}
                else{counter <- counter + 1;}
                }else{
                    indices_mat[k, flag] <- ordered_kl[counter];
                    flag <- flag + 1;
                    counter <- counter + 1;
                }

            }
        }

    } else{
        for(k in 1:dim(theta)[2])
        {
            temp_mat <- KL_score[[k]][,-k];
            if(options=="min"){
                vec <- apply(temp_mat, 1, function(x) min(x))}
            if(options=="max"){
                vec <- apply(temp_mat, 1, function(x) max(x))}

            ordered_kl <- order(vec, decreasing = TRUE);
            counter <- 1
            flag <- counter
            while(flag <= top_features)
            {
                if(counter > dim(theta)[1]){
                  indices_mat[k,(flag:top_features)]=NA;
                  break
                }
                if(!shared){
                    if(which.max(theta[ordered_kl[counter],])==k){
                        indices_mat[k, flag] <- ordered_kl[counter];
                        flag <- flag + 1;
                        counter <- counter + 1;
                    } else {
                        counter <- counter + 1;
                    }
                }else{
                    indices_mat[k, flag] <- ordered_kl[counter];
                    flag <- flag + 1;
                    counter <- counter + 1;
                }
            }
        }
    }

    return(indices_mat);
}
