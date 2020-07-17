#' @title Extracting most highly correlated genes with GoM topics/clusters
#'
#' @description This function compares grades of membership profile for each
#'  cluster in GoM model fit with the data expression profile to identify
#'  genes that are mostly strongly associated with each topic.
#'
#'
#' @param omega \eqn{\boldsymbol{omega}} matrix, the relative grades of memberships
#'                from the GoM model fitting (a \eqn{N x K} matrix where \eqn{N} is
#'                number of samples, \eqn{K} number of topics).
#' @param data \eqn{G x N} matrix of the expression profile of genes across samples,
#'              where \eqn{G} is the number of features and \eqn{N} number of samples
#' @param num_genes The number of top associated genes with each cluster. Defaults to 100
#'
#' @return A list containing two items - a \eqn{K x num_genes} matrix of the
#'          top strongly associated/correlated indices/features for K clusters,
#'          and another \eqn{K x num_genes} matrix of the absolute values of the
#'          correlations.
#'
#' @examples
#'data("MouseDeng2014.FitGoM")
#' omega_mat <- MouseDeng2014.FitGoM$clust_6$omega;
#' read.data1 = function() {
#'     x = tempfile()
#'    download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
#'    z = get(load((x)))
#'    return(z)
#' }
#' Deng2014MouseESC <- read.data1()
#' deng.counts <- Biobase::exprs(Deng2014MouseESC)
#' out <- ExtractHighCorFeatures(omega_mat, deng.counts, num_genes=10)
#'
#' @export
#'
ExtractHighCorFeatures <- function(omega,
                                   data,
                                   num_genes=100){
    if(num_genes > dim(data)[1]){
        stop("num genes must be less than the number of features in data")
    }

    if(dim(data)[2] < dim(omega)[2]){
        stop("number of samples must be greater than K")
    }

    K <- dim(omega)[2]

    indices <- matrix(0,K,num_genes)
    cor_values <- matrix(0,K,num_genes)

    for(k in 1:K){
        abs_cor_genes <- array(0, dim(data)[1])
        for(l in 1:dim(data)[1]){
            if(stats::sd(data[l,])==0){ abs_cor_genes[l] <- 0}
            else{abs_cor_genes[l] <- abs(cor(data[l,], omega[,k]))}
        }
        indices[k,] <- order(abs_cor_genes, decreasing=TRUE)[1:num_genes]
        cor_values[k,] <- abs_cor_genes[indices[k,]]
    }

    ll <- list("feature_indices"=indices,
               "abs_cor_features"=cor_values)
    return(ll)
}
