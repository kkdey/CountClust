#' @title log likelihood of the GoM model fits
#'
#' @description This function takes the \code{FitGoM} output model fits
#'              and compute log likelihood for the GoM models.
#'
#' @param model_output \code{FitGoM} output (a \code{list}).
#'
#' @return logLik_vec a vector of GoM model log likelihoods.
#'
#' @keywords GoM, model fit
#'
#' @examples
#'
#' # We need both the observed count matrix and
#' # the GoM model fit output
#'
#' # Import Deng et al data
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
#' model_logLik(data = deng.counts,
#'              model_output = MouseDeng2014.FitGoM)
#' @export

model_logLik <- function(data, model_output)
{
    # Get the numer GoM models
    num_models <- length(model_output)

    # Compute logLik for each GoM
    loglik_models <- sapply(1:num_models, function(i) {

        theta <- model_output[[i]]$theta
        omega <- model_output[[i]]$omega
        probs <- theta %*% t(omega)

        if (!all.equal(dim(data), dim(probs))) {
            stop("Observed data dimension does not match \n
                 the model fit results dimension")
        }

        loglik <- Reduce(sum, sapply(1:NCOL(data), function(j) {
            dmultinom(x = data[,j], prob = probs[,j], log = TRUE)
        }) )

        return(loglik)
    })
    names(loglik_models) <- names(model_output)

    return(loglik_models)
}
