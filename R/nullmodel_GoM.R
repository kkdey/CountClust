#' @title Null models for Grade of Membership (GoM) cluster validation
#'
#' @description Use null models (popular in ecology) to generate randomized
#' matrix of counts given the observed data matrix, fit the GoM model to these
#' null matrices and compare the fit on null model data with that on the
#' observed data. Used for validating the GoM clusters
#'
#' @param counts The counts matrix (N x G): N- the number of samples, G- number
#' of features
#' @param K  The number of clusters to fit
#' @param tol The tolerance of the GoM model fitted
#' @param null.model The type of nullmodel used (similar to the
#' randomizeMatrix() function argument in picante package)
#' @param iter_fill The number of swaps/fills in each randomized matrix build
#' @param iter_randomized The number of randomization matrices generated
#' @param plot If TRUE, plots density of log Bayes factor
#'
#' @return  Returns a list with
#'        \item{GoMBF.obs}{log BF for the observed counts with K=2 against the
#'        null with no clusters}
#'        \item{GoMBF.rand}{a vector of log BF for each randomized count matrix
#'        with K=2 against the null with no clusters}
#'        \item{pval}{the p-value of the observed log Bayes factor against the
#'        ones from randomized matrices}
#'
#' @importFrom  picante randomizeMatrix
#' @import slam
#' @importFrom stats dmultinom density
#' @export
#'
#' @examples
#'
#' data("ex.counts")
#' nullmodel_GoM(ex.counts,
#'               K=2,
#'               tol=500,
#'               null.model="frequency",
#'               iter_randomized=3,
#'               plot=FALSE)


nullmodel_GoM <- function(counts,
                          K,
                          tol=0.1,
                          null.model=c("frequency", "richness",
                                       "independentswap", "trialswap"),
                          iter_fill=1000,
                          iter_randomized=100,
                          plot=TRUE)
{
    bf_gom_rand <-
        unlist(lapply(1:iter_randomized,
                      function(n)
                      {
                          rand_counts <-
                              picante::randomizeMatrix(counts,
                                                       null.model=null.model,
                                                       iterations=iter_fill);
                          suppressMessages(
                              topics_rand <- FitGoM(rand_counts, K=K, tol=tol))
                          X <- CheckCounts(rand_counts);
                          n <- nrow(X);
                          p <- ncol(X)

                          ## Null model log probability
                          sx <- sum(X)
                          qnull <- slam::col_sums(X)/sx
                          null_prob <- sum( X$v*log(qnull[X$j]) ) -
                              0.5*(n+p)*(log(sx) - log(2*pi))

                          ## GoM model probability

                          gom_model_prob <-
                              sum(unlist(lapply(1:dim(counts)[1],
                                                function(n)
                                                {
                                                    return(dmultinom(counts[n,],
                                                                     prob= topics_rand[[1]]$omega[n,]%*%t(topics_rand[[1]]$theta),
                                                                     log=TRUE))
                                                })))

                          ##

                          bf_gom <- gom_model_prob - null_prob;
                          return(bf_gom)

                      }))

    topics_obs <- FitGoM(counts, K=K, tol=tol);
    X <- CheckCounts(counts); n <- nrow(X); p <- ncol(X)

    ## Null model log probability
    sx <- sum(X)
    qnull <- slam::col_sums(X)/sx
    null_prob <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))

    ## GoM model probability

    gom_model_prob <- sum(unlist(lapply(1:dim(counts)[1], function(n)
    {
        return(dmultinom(counts[n,],
                         prob= topics_obs[[1]]$omega[n,]%*%t(topics_obs[[1]]$theta),
                         log=TRUE))
    })))

    ##

    bf_gom_obs <- gom_model_prob - null_prob;
    pval_bf_gom <- length(which(bf_gom_obs < bf_gom_rand))/iter_randomized;

    ll <- list("GoMBF.obs"=bf_gom_obs,
               "GoMBF.rand"=bf_gom_rand,
               "pval"=pval_bf_gom)
    if(plot){
        plot(density(bf_gom_rand), col="blue",
             main="Log BF density over randomized matrices",
             xlab="Log Bayes facto", ylab="density",
             xlim=c(min(bf_gom_obs,bf_gom_rand)-100, max(bf_gom_obs,bf_gom_rand)+100))
        abline(v=bf_gom_obs, col="red")
    }
    return(ll)
}




CheckCounts <- function(counts){
    if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
    if(is.null(dimnames(counts)[[1]]))
    {dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
    if(is.null(dimnames(counts)[[2]]))
    { dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
    empty <- slam::row_sums(counts) == 0
    if(sum(empty) != 0){
        counts <- counts[!empty,]
        cat(paste("Removed", sum(empty), "blank documents.\n")) }
    return(slam::as.simple_triplet_matrix(counts))
}
