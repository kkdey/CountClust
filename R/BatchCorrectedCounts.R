#' @title Obtaining Batch effect Corrected counts !
#'
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param batch_lab the batch label vector
#' @param use_parallel if TRUE, we do a parallel analysis over featres, else serial application.
#'
#' @description This function first converts counts data to CPM data/ sqrt data, then apply a linear model with the batch effect as a factor. We take the sum of
#'              intercept, residuals and mean batch effect across all the batches and then exponentiate it to bring the data back almost to the original
#'              space of the counts data, except for the batch adjustment.Generate a Poisson random number for each of the means. The obtained counts
#'              are then devoid of batch effects.
#'
#' @author Kushal K Dey, Joyce HSiao
#' @keywords counts data, batch effect
#' @export
#'


BatchCorrectedCounts <- function(data, batch_lab,use_parallel=TRUE, use_cpm=TRUE)
{
  if(use_cpm){
  out_voom <- limma::voom(data);
  trans_data <- out_voom$E;
  }
  if(!use_cpm){
    trans_data <- sqrt(data);
  }
  if(use_parallel){
    batch_removed_counts_mean <- do.call(cbind, parallel::mclapply(1:dim(trans_data)[2],function(g)
                                        {
                                            out <- lm(trans_data[,g] ~  as.factor(batch_lab),
                                                      contrasts = list(batch_lab = "contr.sum") )
                                            return(round(exp(out$coefficients[1] + out$residuals)-0.4))
                                         }, mc.cores=parallel::detectCores()));
  }

  if(!use_parallel){
    batch_removed_counts_mean <- lapply(1:dim(trans_data)[2],function(g)
                                        {
                                          out <- lm(trans_data[,g] ~ batch_lab,
                                                    contrasts = list(batch_lab = "contr.sum"))
                                          return(round(exp(out$coefficients[1] + mean(out$coefficients[-1])+out$residuals)-0.4));
                                          return(ll)
    })
  };

  if (dim(batch_removed_counts_mean)[2]!=dim(data)[2])
    stop("The batch corrected data is not of same dimension as the counts data : try changing use_parallel")
  batch_corrected_counts <- round(batch_removed_counts_mean);
  rownames(batch_corrected_counts) = rownames(data);
  colnames(batch_corrected_counts) = colnames(data);
  return(batch_corrected_counts)
}


