#' title ash type model for f
#'
#' description use ash type model to maxmization
#'
#' @return Ef is the mean of f Ef2 is mean of f^2
#'
#' @keywords internal

# El is expectation of l, and El2 is the second moment of l, sigmae2 is estimation of sigmae^2
ATM_f = function(Y,El,El2,sigmae2){
  sum_El2 = sum(El2)
  sebeta = sqrt( sigmae2/(sum_El2) )
  betahat = (t(El) %*% Y) / (sum_El2)
  # betahat=(sum(l^2))^(-1)*(t(l)%*%Y)
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="+uniform")
  Ef = ATM$PosteriorMean
  SDf = ATM$PosteriorSD
  Ef2 = SDf^2 + Ef^2
  return(list(Ef = Ef, Ef2 = Ef2))
}

#' title ash type model for l
#'
#' description use ash type model to maxmization
#'
#' @return El is the mean of l,  El2 is mean of l^2
#'
#' @keywords internal
#'
ATM_l = function(Y,Ef,Ef2,sigmae2){
  sum_Ef2 = sum(Ef2)
  sebeta = sqrt(sigmae2/(sum_Ef2))
  betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
  # betahat=(sum(f^2))^(-1)*(t(f)%*%t(Y))
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="+uniform")
  El = ATM$PosteriorMean
  SDl = ATM$PosteriorSD
  El2 = SDl^2 + El^2
  return(list(El = El, El2 = El2))
}


# Fval = function()







#' Factor Loading Adaptive Shrinkage (VEM version)
#'
#' flash provide rank one matrix decomposition
#'
#' @param Y is the data matrix (N by P)
#' @param tol which is the stop criterion for the convergence, default is 1e-5
#' @param numtau number of iteration, default is 500. for the backfitting case, the number of tau should be 5 or 10.
#'
#' @details flash_VEM privide rank one matrix decomposition with variational EM algorithm.
#'
#' @export flash_VEM
#'
#' @importFrom ashr ash
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{l}} {is a N vector for loadings}
#'   \item{\code{f}} {is a P vector for factors}
#'   \item{\code{sigmae2}} {is mean of sigma square which is estimation for the noise variance}
#'  }
#' @examples
#' N = 100
#' P = 200
#' Y = matrix(rnorm(N*P,0,1),ncol=P)
#' g = flash_VEM(Y)
#'


# set the number of iteration as numtau
flash_VEM = function(Y,tol=1e-6,numtau = 500){
  #dealing with missing value
  Y[is.na(Y)] = 0
  # get initial value for l and f and sigmae
  El = svd(Y)$u[,1]
  El2 = El^2
  Ef = as.vector(t(El)%*%Y)
  Ef2 = Ef^2

  #start iteration
  sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

  par_f = ATM_f(Y,El,El2,sigmae2_v)
  Ef = par_f$Ef
  Ef2 = par_f$Ef2

  sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
  #sigmae2
  par_l = ATM_l(Y,Ef,Ef2,sigmae2_v)
  El = par_l$El
  El2 = par_l$El2

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < numtau ){
    tau = tau + 1
    pre_sigmae2 = sigmae2_v

    sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )

    par_f = ATM_f(Y,El,El2,sigmae2_v)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    if(sum(Ef^2)==0){
      l = 0
      f = 0
      break
    }
    sigmae2_v = mean( Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2)) )
    #sigmae2
    par_l = ATM_l(Y,Ef,Ef2,sigmae2_v)
    El = par_l$El
    El2 = par_l$El2
    if(sum(El^2)==0){
      l = 0
      f = 0
      break
    }
    epsilon = abs(pre_sigmae2 - sigmae2_v )
  }
  return(list(l = El, f = Ef, sigmae2 = sigmae2_v))
}
