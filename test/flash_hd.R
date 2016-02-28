# flash for poison data (possible) and non constant variance

#' title ash type model for f
#'
#' description use ash type model to maxmization
#'
#' @return Ef is the mean of f,  Ef2 is mean of f^2
#'
#' @keywords internal
#'
# El is expectation of l, and El2 is the second moment of l, sigmae2 is estimation of sigmae^2
ATM_fhd = function(Y,El,El2,sigmae2){
  sum_El2 = t(El2) %*% (1/sigmae2)
  sum_El2 = as.vector(sum_El2)
  sebeta = sqrt( 1/(sum_El2) )
  betahat = as.vector( t(El) %*% (Y/sigmae2) ) / (sum_El2)
  # betahat=(sum(l^2))^(-1)*(t(l)%*%Y)
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="normal")
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
ATM_lhd = function(Y,Ef,Ef2,sigmae2){
  sum_Ef2 = (1/sigmae2) %*% Ef2
  sum_Ef2 = as.vector(sum_Ef2)
  sebeta = sqrt(1/(sum_Ef2))
  betahat = as.vector( (Y/sigmae2) %*% Ef ) / (sum_Ef2)
  # betahat = as.vector(t(Ef) %*% t(Y/sigmae2)) / (sum_Ef2)
  # betahat=(sum(f^2))^(-1)*(t(f)%*%t(Y))
  betahat=as.vector(betahat)
  ATM = ash(betahat, sebeta, method="fdr", mixcompdist="normal")
  El = ATM$PosteriorMean
  SDl = ATM$PosteriorSD
  El2 = SDl^2 + El^2
  return(list(El = El, El2 = El2))
}


# Fval = function()




#' title estimate the variance using anova
#'
#' sigma_ij = mu + a_i + b_j
#'
#' @return sigmae2_v fitted value from the anova model
#'
#' @keywords internal
#'
# function for nonconstant variance
sigmaest = function(residualsqr){
  N = dim(residualsqr)[1]
  P = dim(residualsqr)[2]
  # r_ij^2 = mu + a_i + b_j
  mu = mean(residualsqr)
  a = apply(residualsqr,1,mean) - mu
  b = apply(residualsqr,2,mean) - mu
  # add ash in to the a and b
  # a = ash(a,sd(a))$PosteriorMean
  # b = ash(b,sd(b))$PosteriorMean
  # a = a - mean(a)
  # b = b - mean(b)

  sigmae2_v = mu + matrix(rep(a,P),ncol = P) + matrix(rep(b,each = N),ncol = P)
  return(sigmae2_v)
}
# another contrast
#sigmaest = function(residualsqr){
#  N = dim(residualsqr)[1]
#  P = dim(residualsqr)[2]
#  # r_ij^2 = mu + a_i + b_j
#  a = rep(0,N)
#  b = rep(0,P)
#  a[2:N] = sapply(seq(2,N),function(x){mean(residualsqr[x,] - residualsqr[1,])})
#  b[2:P] = sapply(seq(2,P),function(x){mean(residualsqr[,x] - residualsqr[,1])})
#  mu = (residualsqr[1,1]+sum(residualsqr[2:N,1]-a[2:N]) + sum(residualsqr[1,2:P]-b[2:P]))/(N+P-1)
#  sigmae2_v = mu + matrix(rep(a,P),ncol = P) + matrix(rep(b,each = N),ncol = P)
#  return(sigmae2_v)
#}


#' Factor Loading Adaptive Shrinkage (heteroscedasticity version)
#'
#' flash provide rank one matrix decomposition
#'
#' @param Y is the data matrix (N by P)
#' @param tol which is the stop criterion for the convergence, default is 1e-5
#' @param numtau number of iteration, default is 500. for the backfitting case, the number of tau should be 5 or 10.
#' @param sigmae2 which is a matrix for the known noise variance matrix. the default is 1 which is not useful
#' @param partype type for the noise variance, which takes values from "constant": constant variance,
#'  "anova": variance is from anova model, "loganova": log variance is from anova model, "poisson":
#'  variance is equal to mean, else user should provide known value for sigmae2
#'
#' @details flash_hd privide rank one matrix decomposition with variational EM algorithm.
#'
#' @export flash_hd
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
flash_hd = function(Y, tol=1e-5, numtau = 500, partype = "constant", sigmae2 = 1){
  N = dim(Y)[1]
  P = dim(Y)[2]
  #dealing with missing value
  Y[is.na(Y)] = 0
  # get initial value for l and f and sigmae
  El = svd(Y)$u[,1]
  El2 = El^2
  Ef = as.vector(t(El)%*%Y)
  Ef2 = Ef^2

  #start iteration
  # residual matrix initialization
  sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  if(partype == "anova"){
    sigmae2 = sigmaest(sigmae2_v)
  }else if(partype == "loganova"){
    sigmae2 = exp(sigmaest(log(sigmae2_v)))
  }else if( partype == "poisson"){
    # where variance is equal to mean
    sigmae2 = El %*% t(Ef)
  } else if(partype == "constant") {
    # this is for constant variance which can be easier
    sigmae2 = mean(sigmae2_v) * matrix(rep(1,N*P),ncol = P)
  }else{
    sigmae2 = sigmae2
  }

  par_f = ATM_fhd(Y,El,El2,sigmae2)
  Ef = par_f$Ef
  Ef2 = par_f$Ef2

  sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  if(partype == "anova"){
    sigmae2 = sigmaest(sigmae2_v)
  }else if(partype == "loganova"){
    sigmae2 = exp(sigmaest(log(sigmae2_v)))
  }else if( partype == "poisson"){
    # where variance is equal to mean
    sigmae2 = El %*% t(Ef)
  } else if(partype == "constant") {
    # this is for constant variance which can be easier
    sigmae2 = mean(sigmae2_v) * matrix(rep(1,N*P),ncol = P)
  }else{
    sigmae2 = sigmae2
  }

  par_l = ATM_lhd(Y,Ef,Ef2,sigmae2)
  El = par_l$El
  El2 = par_l$El2

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < numtau){
    tau = tau + 1
    if(partype == "constant"){
      pre_clkhd = mean(sigmae2)
    }else {
      pre_clkhd = mean(sigmae2_v/sigmae2)
    }

    sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
    if(partype == "anova"){
      sigmae2 = sigmaest(sigmae2_v)
    }else if(partype == "loganova"){
      sigmae2 = exp(sigmaest(log(sigmae2_v)))
    }else if( partype == "poisson"){
      # where variance is equal to mean
      sigmae2 = El %*% t(Ef)
    } else if(partype == "constant") {
      # this is for constant variance which can be easier
      sigmae2 = mean(sigmae2_v) * matrix(rep(1,N*P),ncol = P)
    }else{
      sigmae2 = sigmae2
    }

    par_f = ATM_fhd(Y,El,El2,sigmae2)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    if(sum(Ef^2)==0){
      l = 0
      f = 0
      break
    }
    sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
    if(partype == "anova"){
      sigmae2 = sigmaest(sigmae2_v)
    }else if(partype == "loganova"){
      sigmae2 = exp(sigmaest(log(sigmae2_v)))
    }else if( partype == "poisson"){
      # where variance is equal to mean
      sigmae2 = El %*% t(Ef)
    } else if(partype == "constant") {
      # this is for constant variance which can be easier
      sigmae2 = mean(sigmae2_v) * matrix(rep(1,N*P),ncol = P)
    }else{
      sigmae2 = sigmae2
    }

    par_l = ATM_lhd(Y,Ef,Ef2,sigmae2)
    El = par_l$El
    El2 = par_l$El2
    if(partype == "constant"){
      clkhd = mean(sigmae2)
    }else {
      clkhd = mean(sigmae2_v/sigmae2)
    }
    epsilon = abs(pre_clkhd - clkhd )
  }
  return(list(l = El, f = Ef, sigmae2 = sigmae2))
}
