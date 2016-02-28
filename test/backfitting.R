#' backfitting for correction K factors model
#'
#' correction the factor and loading matrix estimation using backfitting algorithm
#'
#' @param Y is the data matrix (N by P)
#' @param Lest is estimation for l to correct
#' @param Fest is estimation for f to correct
#' @param tautol this the stop criterion for convergence, number for iteration
#' @param numtau number for iteration in flash rank one.
#'
#' @details greedy algorithm on the residual matrix to get a rank one matrix decomposition
#'
#' @export backfitting
#'
#' @importFrom ashr ash
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{l}} {is a K by N matrix for loadings}
#'   \item{\code{f}} {is a K by P matrix for factors}
#'  }
#' @examples
#' N = 100
#' P = 200
#' Y = matrix(rnorm(N*P,0,1),ncol=P)
#' g = flash_v_K(Y,10)
#'

# tautol is the number of iterations here
backfitting = function(Y,Lest,Fest,tautol=100,numtau = 500){
  # backfitting with initial values
  epsilon = 1
  tau = 1
  while(epsilon>1e-5 & tau < tautol){
    tau = tau + 1
    # if the l or f is a vector
    if(is.vector(Lest) || is.vector(Fest)){
      residual = Y - Lest %*% t(Fest)
      preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
      residual = residual + Lest %*% t(Fest)
      r_flash = flash_VEM(residual,numtau = numtau)
      Lest = r_flash$l
      Fest = r_flash$f
      residual = residual - Lest %*% t(Fest)
    }else{
      K = dim(Lest)[2]
      # this one can be put out of the while loop
      residual = Y - Lest %*% t(Fest)
      preRMSfl = sqrt(mean((Lest %*% t(Fest))^2))
      for(i in 1:K){
        residual = residual + Lest[,i] %*% t(Fest[,i])
        r_flash = flash_VEM(residual,numtau = numtau)
        Lest[,i] = r_flash$l
        Fest[,i] = r_flash$f
        residual = residual - Lest[,i] %*% t(Fest[,i])
      }
      # remove the zero in the l and f
      while(i <= dim(Lest)[2] ){
        if(sum((Lest[,i])^2)==0 || sum((Fest[,i])^2)==0){
          Lest = Lest[,-i]
          Fest = Fest[,-i]
        }
        numfactor = ifelse(is.vector(Lest),1,dim(Lest)[2])
        if(numfactor == 1){
          break
        }
        i = i+1
      }
    }

    RMSfl = sqrt(mean((Lest %*% t(Fest))^2))
    epsilon = abs(preRMSfl - RMSfl)
  }
  return(list(l = Lest, f = Fest))
}
