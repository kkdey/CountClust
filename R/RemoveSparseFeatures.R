#' @title Removes features with a lot of zero counts
#'
#' @description This function deals with zero counts in the counts
#'   dataset. If for a feature, the proportion of zeros across the
#'   samples is greater than filter_prop, then we remove the feature.
#'
#' @param data Count data in a sample by feature matrix.
#' 
#' @param filter_prop Threshold proportion. If the proportion of zeros
#'   for the feature exceeds this threshold then we remove the feature
#'   altogether. Default is 0.9.
#'
#' @keywords counts data, feature extraction
#'
#' @return Returns a list with
#'    \item{data}{filtered data with sparse features removed}
#'    \item{sparse_features}{the feature names of the features
#'                           found sparse and removed}
#' 
#' @export
#'
#' @examples
#' mat <- rbind(c(2,0,3,0,4),c(4,5,5,0,0),c(30,34,63,25,0),c(0,0,0,0,0))
#' RemoveSparseFeatures(mat, filter_prop = 0.5)
#' RemoveSparseFeatures(mat)
#'
RemoveSparseFeatures <- function(data, filter_prop=0.9)
{
    zero_prop_cols <- apply(data,2,function(x) length(which(x==0))/length(x))

    if(max(zero_prop_cols)<filter_prop){
        message("no features filtered due to sparsity",
                domain = NULL, appendLF = TRUE)
        ll <- list("data" = data,
                   "sparse_features" = as.numeric())
        return(ll)
    } else{
        features_to_remove <- as.numeric(which(zero_prop_cols > filter_prop))
        data_filtered <- as.matrix(data[,-features_to_remove])
        ll <- list("data" = data_filtered,
                   "sparse_features" = colnames(data)[features_to_remove])
        return(ll)
    }
}
