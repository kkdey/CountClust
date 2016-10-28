#' @title Run Grade of Membership (GoM) model with multiple starting points !
#'
#' @description Fits grade of membership model \code{FitGoM()} to count data
#'        with multiple starting points and choose the best fit using BIC (Bayesian
#'        Information Criterion). the multiple starting points ensure that
#'        the output is more reliable.
#'
#' @param data counts data \eqn{N x G}, with \eqn{N}, the number of samples
#'       along the rows and \eqn{G}, number of genes along columns.
#' @param K the vector of clusters or topics to be fitted. Must be an integer,
#'       unlike in ]\code{FitGom()}. So you need to apply this function separately
#'       for each K.
#' @param tol Tolerance value for GoM model absolute log posterior increase
#'            at successive iterations (set to 0.1 as default).
#' @param burn_trials The number of trials with different starting points used.
#' @param path_rda The directory path for saving the GoM model output.
#'                  If NULL, it will return the output to console.
#' @param control Control parameters. Same as topics() function of
#'                 maptpx package.
#'
#' @return Outputs the best GoM model fit output for cluster K and saves it
#'         at the directory path in path_rda if the latter is provided.
#'
#' @references Matt Taddy. On Estimation and Selection for Topic Models.
#'                AISTATS 2012, JMLR W\&CP 22.
#'
#'            Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly.
#'                Inference of population structure using multilocus genotype
#'                data. Genetics 155.2 (2000): 945-959.
#'
#' @keywords counts data, clustering, Structure plot
#'
#'
#' @examples
#'
#' data("ex.counts")
#' out <- FitGoMpool(ex.counts, K=2, tol=100, burn_trials=3,
#'                    control=list(tmax=100))
#'
#' @importFrom maptpx topics
#' @import slam
#' @importFrom utils modifyList
#' @export
#'

FitGoMpool <- function(data,
                   K,
                   tol=0.1,
                   burn_trials = 10,
                   path_rda = NULL,
                   control=list())
{
  if(length(K) > 1)
    stop("For FitGoMpool, K must be an integer, run for separate K")

  out <- list()

  for(num in 1:burn_trials){
    out[[num]] <- FitGoM(data, K=K, tol=tol)
  }

  BIC_val <- array(0, burn_trials)
  for(n in 1:length(BIC_val)){
    BIC_val [n] <- compGoM(data, out[[n]])[,1]$BIC
  }

  Topic_clus <- out[[which.min(BIC_val)]][[1]]

  if(!is.null(path_rda)){
    save(Topic_clus, file = path_rda);
  }else{
    return(Topic_clus)
  }
}

