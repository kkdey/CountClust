#' @title Re-ordering topic proportion matrices and Information calculation
#'
#' @param omega1 topic proportion fitted matrix (N x K) from model 1:
#' @param omega2 topic proportion fitted matrix (N x K) from model 2:

#' @description This function computes a re-ordering of the clusters from GoM model fit in
#' in one model to make it comparable with that from another. The two models are applie on the
#' same set of samples with same number of clusters, but the features may change from one model to another.
#' The two models may not be of same type as well. One could be a DAPC model, the other a standard topic model.
#' The idea is to check for consistency in topic proportion patterns across multiple GoM methods or across
#' different types of feature sets.
#'
#' @author Kushal K Dey, Matthew Stephens
#'
#' @export

compare_omega <- function(omega1, omega2)
{
    cor.out <- 1 - cor(omega1, omega2)
    kl.out <- matrix(0,dim(omega1)[2],dim(omega2)[2]);
    for(m in 1:dim(omega1)[2])
    {
      for(n in 1:dim(omega2)[2])
      {
        kl.out[m,n] <- 0.5* philentropy::distance(t(cbind(docweights1[,m], docweights2[,n])), method="kullback-leibler")
        + 0.5* philentropy::distance(t(cbind(docweights2[,n], docweights1[,m])), method="kullback-leibler") ;
      }
    }

    kl.order_model2_topics <- apply(kl.out, 1, function(x) which.min(x))
    kl.information_content <- exp(-(mean(apply(kl.out, 1, min))/ mean(kl.out)));

    cor.order_model2_topics <- apply(cor.out, 1, function(x) which.min(x))
    cor.information_content <- exp(-(mean(apply(cor.out, 1, min))/ mean(cor.out)));

    ll <- list("kl.dist"=kl.out, "kl.order_model2_topics"=kl.order_model2_topics,
                "kl.information_content"=kl.information_content,
                "cor.dist"=cor.out, "cor.order_model2_topics"=cor.order_model2_topics,
               "cor.information_content"=cor.information_content);
    return(ll)
}
