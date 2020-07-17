#' @title Compare GoM model fits across K or across different runs
#'   through log-likelihood, BIC and null loglikelihood
#'
#' @description This function takes the \code{FitGoM}/maptpx fitted
#'   model and computes log likelihood, BIC and null model loglikelihood
#'   for the fitted GoM models.
#'
#' @param data Matrix on which GoM model is fitted (samples along
#'   rows, genes along columns)
#' 
#' @param model \code{FitGoM} or\code{maptpx::topics} function output
#'   (either a class \code{topics} or a \code{list} of class
#'   \code{topics}).
#'
#' @return compGoM_models a vector list that returns the BIC and
#'   loglikelihood values for each of the fitted models in \code{model}.
#'
#' @keywords GoM, model fit
#'
#' @examples
#' read.data <- function() {
#'   x <- tempfile()
#'   download.file(paste0("https://cdn.rawgit.com/kkdey/",
#'                          "singleCellRNASeqMouseDeng2014",
#'                          "/master/data/Deng2014MouseEsc.rda"),
#'                 destfile = x, quiet = TRUE)
#'   z <- get(load((x)))
#'   return(z)
#'   }
#' Deng2014MouseESC <-read.data()
#'
#' # Extract observed counts
#' deng.counts <- Biobase::exprs(Deng2014MouseESC)
#'
#' # Import GoM fitting results
#' data("MouseDeng2014.FitGoM")
#' names(MouseDeng2014.FitGoM)
#'
#' compGoM(data = t(deng.counts),
#'            model = MouseDeng2014.FitGoM)
#' compGoM(data = t(deng.counts),
#'            model = MouseDeng2014.FitGoM$clust_3)
#'
#' @importFrom slam col_sums
#' @importFrom slam row_sums
#' @importFrom slam as.simple_triplet_matrix
#' 
#' @export
#'
compGoM <- function(data, model) {
    X <- CheckCounts(data+1e-04)
    p <- ncol(X)
    n <- nrow(X)

    ## Null model log probability
    sx <- sum(X)
    qnull <- col_sums(X)/sx
    null_loglik <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))

    if(class(model) == "topics"){
        theta <- model$theta
        omega <- model$omega
        probs <- omega %*% t(theta)

        if (!all.equal(dim(data), dim(probs))) {
            stop("Observed data dimension does not match \n
                 the model fit results dimension")
        }

        loglik <- Reduce(sum, sapply(1:NROW(data), function(i) {
            dmultinom(x = data[i,], prob = probs[i,], log = TRUE)
        }) )
        BIC <- -2*loglik + NROW(theta)*(NCOL(theta)-1)*log(NROW(data))
        ll <- list("BIC"=BIC, "loglik"=loglik, "null_loglik"=null_loglik)
        return(ll)
    }
    if(class(model) == "list"){
        compgom_list <- list()
        num_models <- length(model)
        for(j in 1:num_models){
            theta <- model[[j]]$theta
            omega <- model[[j]]$omega
            probs <- omega %*% t(theta)

            if (!all.equal(dim(data), dim(probs))) {
                stop("Observed data dimension does not match \n
                     the model fit results dimension")
            }

            loglik <- Reduce(sum, sapply(1:NROW(data), function(i) {
                dmultinom(x = data[i,], prob = probs[i,], log = TRUE)
            }) )
            BIC <- -2*loglik + NROW(theta)*(NCOL(theta)-1)*log(NROW(data))
            compgom_list[[j]] <- list("BIC"=BIC, "loglik"=loglik, "null_loglik"=null_loglik)
        }
        if(!is.null(names(model))) names(compgom_list) <- names(model)
        return(compgom_list)
    }
    else{
        stop("The model argument is either of class topics or of class list- where each element of the list is of
             class topics")
    }
}
CheckCounts <- function(fcounts){
    if(class(fcounts)[1] == "TermDocumentMatrix"){ fcounts <- t(fcounts) }
    if(is.null(dimnames(fcounts)[[1]])){ dimnames(fcounts)[[1]] <- paste("doc",1:nrow(fcounts)) }
    if(is.null(dimnames(fcounts)[[2]])){ dimnames(fcounts)[[2]] <- paste("wrd",1:ncol(fcounts)) }
    empty <- row_sums(fcounts) == 0
    if(sum(empty) != 0){
        fcounts <- fcounts[!empty,]
        cat(paste("Removed", sum(empty), "blank documents.\n")) }
    return(as.simple_triplet_matrix(fcounts))
}
