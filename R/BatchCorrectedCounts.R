#' @title Obtain Batch effect Corrected counts !
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param batch_lab the batch label vector
#' @param use_parallel if TRUE, we do a parallel analysis over featres, else serial application.
#'
#' @description This function first converts counts data to log CPM data , then apply a linear model with the batch effect as a factor. We take the sum of
#'              intercept, residuals and mean batch effect across all the batches and then inverse transform it back to counts to get rid of batch effects.
#'
#' @return Returns a counts data. with same dimension as the input data, but which is corrected for batch_lab.
#' @author Kushal K Dey, Joyce HSiao
#' @keywords counts data, batch effect
#' @examples
#' N=500;
#' K=4;
#' G=100;
#' Label.Batch=c(rep(1,N/4),rep(2,N/4),rep(3,N/4),rep(4,N/4));
#' alpha_true=matrix(rnorm((K)*G,0.5,1),nrow=(K));
#' library(gtools)
#' T=10;
#' omega_true=matrix(rbind(gtools::rdirichlet(T*10,c(3,4,2,6)),gtools::rdirichlet(T*10,c(1,4,6,3)),
#'                       gtools::rdirichlet(T*10,c(4,1,2,2)),gtools::rdirichlet(T*10,c(2,6,3,2)),
#'                       gtools::rdirichlet(T*10,c(3,3,5,4))), nrow=N);
#' B=max(Label.Batch);
#' sigmab_true=2;
#' beta_true=matrix(0,B,G);
#' for(g in 1:G)
#' {
#'    beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
#' }
#' read_counts=matrix(0,N,G);
#' for(n in 1:N){
#'    for(g in 1:G)
#'      {
#'         read_counts[n,g]=rpois(1, omega_true[n,]%*%exp(alpha_true[,g]
#'                                    + beta_true[Label.Batch[n],g]));
#'      }
#'  }
#'
#'  batchcorrect.counts <- BatchCorrectedCounts(read_counts, Label.Batch, use_parallel=FALSE);
#'
#' @export
#'


BatchCorrectedCounts <- function(data, batch_lab,use_parallel=TRUE)
{
  out_voom <- limma::voom(data);
  trans_data <- out_voom$E;
  lib_size <- rowSums(data);
  if(use_parallel){
    batch_removed_counts_mean <- do.call(cbind, parallel::mclapply(1:dim(trans_data)[2],function(g)
                                        {
                                            out <- lm(trans_data[,g] ~  as.factor(batch_lab),
                                                      contrasts = list(batch_lab = "contr.sum") )
                                            return(round(exp((out$coefficients[1] + out$residuals)/6)*(lib_size+1)-0.4))
                                         }, mc.cores=parallel::detectCores()));
  }

  if(!use_parallel){
    batch_removed_counts_mean <- do.call(cbind, lapply(1:dim(trans_data)[2],function(g)
                                        {
                                          out <- lm(trans_data[,g] ~ as.factor(batch_lab),
                                                    contrasts = list(batch_lab = "contr.sum"))
                                          return(round(exp((out$coefficients[1] + out$residuals)/6)*(lib_size+1)-0.4));
                                        }));
  }

  if (dim(batch_removed_counts_mean)[2]!=dim(data)[2])
    stop("The batch corrected data is not of same dimension as the counts data : try changing use_parallel")
  batch_corrected_counts <- round(batch_removed_counts_mean);
  rownames(batch_corrected_counts) = rownames(data);
  colnames(batch_corrected_counts) = colnames(data);
  return(batch_corrected_counts)
}


