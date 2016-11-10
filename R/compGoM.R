#' @title compGoM: compare GoM model fits through log-likelihood, BIC and null loglikelihood
#'
#' @description This function takes the \code{FitGoM} output model fits
#'              and compute log likelihood, BIC and null model loglikelihood for the GoM models.
#'
#' @param data matrix on which GoM model is fitted (samples along rows, genes along columns)
#' @param model_output \code{FitGoM} output (a \code{list}).
#'
#' @return compGoM_models a vector of GoM model fit BIC, loglikelihood and null model loglikelihood for each model in FiGoM model input.
#'
#' @keywords GoM, model fit
#'
#' @examples
#'
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
#'            model_output = MouseDeng2014.FitGoM)
#' @importFrom slam col_sums row_sums as.simple_triplet_matrix
#' @export

compGoM <- function(data, model_output)
{
    # Get the numer GoM models
    num_models <- length(model_output)

    X <- CheckCounts(data+1e-04);
    p <- ncol(X)
    n <- nrow(X)

    ## Null model log probability
    sx <- sum(X)
    qnull <- slam::col_sums(X)/sx
    null_loglik <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))


    # Compute logLik for each GoM
    compGoM_models <- sapply(1:num_models, function(i) {

        theta <- model_output[[i]]$theta
        omega <- model_output[[i]]$omega
        probs <- omega %*% t(theta)

        if (!all.equal(dim(data), dim(probs))) {
            stop("Observed data dimension does not match \n
                 the model fit results dimension")
        }

        loglik <- Reduce(sum, sapply(1:NROW(data), function(i) {
            dmultinom(x = data[i,], prob = probs[i,], log = TRUE)
        }) )
        BIC <- -2*loglik + NROW(theta)*(NCOL(theta)-1)*log(NROW(data))

        ll <- list("BIC"=BIC, "loglik"=loglik, "null_loglik"=null_loglik);

        return(ll)
    })
    names(compGoM_models) <- names(model_output)

    return(compGoM_models)
}

CheckCounts <- function(fcounts){
    if(class(fcounts)[1] == "TermDocumentMatrix"){ fcounts <- t(fcounts) }
    if(is.null(dimnames(fcounts)[[1]])){ dimnames(fcounts)[[1]] <- paste("doc",1:nrow(fcounts)) }
    if(is.null(dimnames(fcounts)[[2]])){ dimnames(fcounts)[[2]] <- paste("wrd",1:ncol(fcounts)) }
    empty <- slam::row_sums(fcounts) == 0
    if(sum(empty) != 0){
        fcounts <- fcounts[!empty,]
        cat(paste("Removed", sum(empty), "blank documents.\n")) }
    return(slam::as.simple_triplet_matrix(fcounts))
}
