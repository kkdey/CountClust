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
#' @param options the measure used to choose best fit, either "BF" or "BIC" measures can be used.
#'        BF is more trustworthy, but BIC can be used for better model comparison.
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
                   options = c("BF", "BIC"),
                   path_rda = NULL,
                   control=list())
{
  if(length(K) > 1)
    stop("For FitGoMpool, K must be an integer, run for separate K")

  out <- list()

  for(num in 1:burn_trials){
    out[[num]] <- FitGoM(data, K=K, tol=tol)
  }

  if(options=="BIC"){
        BIC_val <- array(0, burn_trials)
        for(n in 1:length(BIC_val)){
            BIC_val [n] <- compGoM(data, out[[n]])[,1]$BIC
        }

        Topic_clus <- out[[which.min(BIC_val)]][[1]]
        ll <- list("topic_fit"=Topic_clus,
                   "BIC"=BIC_val[which.min(BIC_val)])
  }

  if(options=="BF"){
      BF_val <- array(0, burn_trials);
      for(n in 1:length(BF_val)){
          BF_val [n] <- as.numeric(out[[n]][[1]]$BF);
      }

      Topic_clus <- out[[which.max(BF_val)]][[1]]
      ll <- list("topic_fit"=Topic_clus,
                 "BF"=BF_val[which.max(BF_val)])
  }


  if(!is.null(path_rda)){
    save(Topic_clus, file = path_rda);
     return(ll)
  }else{
     return(ll)
  }
}

