#' @title Re-ordering cluster membership proportion matrices and Information
#'      calculation
#'
#' @description This function computes a re-ordering of the clusters from GoM
#'      model fit in one model to make it comparable with that from another.
#'      The two models are applied on the same set of samples with same number
#'      of clusters, but the features may change from one model to another.
#'      The two models may not be of same type as well. One could be a DAPC
#'      model, the other a standard topic model. Aids in checking for
#'      consistency in topic proportion patterns across multiple GoM methods
#'      or across different types of feature sets.
#'
#' @param omega1 cluster membership proportion matrix (N x K) from model 1
#' @param omega2 cluster membership proportion matrix (N x K) from model 2
#'
#'
#' @return  Returns a list containing
#'            \item{kl.dist}{A symmetric KL divergence matrix across the
#'                            re-ordered clusters of two omega matrices}
#'            \item{kl.order_model2_topics}{re-ordering of the clusters
#'                                for omega2 to match the clusters for omega1
#'                                based on KL divergence}
#'            \item{kl.information_content}{A measure based on KL information
#'                                to record how much information in omega2
#'                                is explained by omega1. Varies from 0 to 1}
#'            \item{cor.dist}{A correlation matrix across the re-ordered
#'                                 clusters of two omega matrices}
#'            \item{cor.order_model2_topics}{re-ordering of the clusters
#'                                 for omega2 to match the clusters for
#'                                  omega1 based on correlation information}
#'            \item{cor.information_content}{A measure based on correlation
#'                                 information to record how much information
#'                                 in omega2 is explained by omega1. Varies from 0 to 1}
#'
#' @examples
#' tt=10;
#' omega1=matrix(rbind(gtools::rdirichlet(tt*10,c(3,4,2,6)),
#'                      gtools::rdirichlet(tt*10,c(1,4,6,3)),
#'                       gtools::rdirichlet(tt*10,c(4,1,2,2))), nrow=3*10*tt);
#' omega2=matrix(rbind(gtools::rdirichlet(tt*10,c(1,2,4,6)),
#'                        gtools::rdirichlet(tt*10,c(1,4,6,3)),
#'                       gtools::rdirichlet(tt*10,c(3,1,5,2))), nrow=3*10*tt);
#' out <- compare_omega(omega1, omega2)
#'
#' @import gtools
#' @export

compare_omega <- function(omega1, omega2)
{
    omega1[omega1==0] <- 1e-20;
    omega2[omega2==0] <- 1e-20;
    cor.out <- 1 - cor(omega1, omega2)
    kl.out <- matrix(0,dim(omega1)[2],dim(omega2)[2]);
    for(m in 1:dim(omega1)[2])
    {
      for(n in 1:dim(omega2)[2])
      {
        KLdiv.mat <- flexmix::KLdiv(as.matrix(cbind(omega1[,m],omega2[,n])));
        kl.out[m,n] <- 0.5* (KLdiv.mat[1,2]+ KLdiv.mat[2,1]);
      }
    }

    kl.order_model2_topics <- apply(kl.out, 1, function(x) which.min(x))
    kl.information_content <- exp(-(mean(apply(kl.out, 1, min))/ mean(kl.out)));

    cor.order_model2_topics <- apply(cor.out, 1, function(x) which.min(x))
    cor.information_content <- exp(-(mean(apply(cor.out, 1, min))/ mean(cor.out)));

    ll <- list("kl.dist"=kl.out,
               "kl.order_model2_topics" = kl.order_model2_topics,
               "kl.information_content"=kl.information_content,
               "cor.dist"=cor.out,
               "cor.order_model2_topics"=cor.order_model2_topics,
               "cor.information_content"=cor.information_content);
    return(ll)
}
