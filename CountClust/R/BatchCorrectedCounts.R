#' @title Obtaining Batch effect Corrected counts!
#'
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param expand_factor a expanding factor required to convert batch corrected CPM values to batch corrected counts (default=100)
#' @param batch_labs the batch label vector
#' 
#' 
#' @description This function first converts counts data to CPM data, then apply a linear model with the batch effect as a factor. The batch corrected 
#'              CPM residuals are then exponentiated and multiplied by expand_factor. These values for each sample and each feature is taken as the 
#'              mean for a Poisson random number generation. The obtained counts are then devoid of batch effects.
#'
#' @author Kushal K Dey, Joyce HSiao
#' @keywords counts data, batch effect
#' @export
#' 


BatchCorrectedCounts <- function(data, expand_factor, batch_labs)
{
  cpm_data <- voom(data)$E;
  batch_removed_cpm <- matrix(0,dim(cpm_data)[1], dim(cpm_data)[2]);
  if(use_parallel){
    batch_removed_cpm <- do.call(cbind,mclapply(1:dim(cpm_data)[2],function(g) as.numeric(lm(cpm_data[,g] ~ individual_id+batch_id)$residuals),mc.cores=detectCores()));
  }
  if(!use_parallel){
    batch_removed_cpm <- do.call(cbind,lapply(1:dim(cpm_data)[2],function(g) as.numeric(lm(cpm_data[,g] ~ individual_id+batch_id)$residuals)));
  }
  
  if (dim(batch_removed_cpm)[2]!=dim(cpm_data)[2])
    stop("The batch corrected data is not of same dimension as the counts data : try changing use_parallel")
  
  exp_batch_removed_cpm_data <- expand_factor*exp(batch_removed_cpm);
  batch_corrected_counts <- apply(exp_batch_removed_cpm_data,c(1,2),function(x) rpois(1,x));
  return(batch_corrected_counts)
}


