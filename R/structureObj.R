#' @title Topic model fitting and Structure Plot!
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param nclus_vec the vector of clusters or topics to be fitted.
#' @param tol the tolerance value for topic model fit (set to 0.001 as default)
#' @param path The directory path where we want to save the data and Structure plots.

#' @description This function takes the counts data (no. of samples x no. of features) and the value of K, the number of topics or
#' cluster to fit, along with sample metadata information and fits the topic model (due to Matt Taddy, check package
#' maptpx) and outputs the Structure plot with the column bars in
#' Structure plot arranged as per the sample metadata information.
#'
#' @author Kushal K Dey, Matt Taddy, Matthew Stephens, Ida Moltke
#'
#' @references Matt Taddy.On Estimation and Selection for Topic Models. AISTATS 2012, JMLR W\&CP 22.
#'
#'            Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly. Inference of population structure using multilocus genotype data.
#'            Genetics 155.2 (2000): 945-959.
#' @keywords counts data, clustering, Structure plot
#'
#' @export
#'



StructureObj <- function(data, nclus_vec, tol, path_rda)
{
  ## dealing with blank rows: we first remove them

  indices_blank <- as.numeric(which(apply(data,1,max)==0));
  if(length(indices_blank)!=0){
  data <- as.matrix(data[-indices_blank,]);
  }

  message('Fitting the topic model (due to Matt Taddy)', domain = NULL, appendLF = TRUE)

  Topic_clus_list <- lapply(nclus_vec, function(per_clust) {
    suppressWarnings(out <- maptpx::topics(as.matrix(data), K = per_clust, tol=tol))
    return(out)
  })

  names(Topic_clus_list) <- paste0("clust_",nclus_vec)
  save(Topic_clus_list, file = path_rda);

#  if(plot) {
#  message('Creating the Structure plots', domain = NULL, appendLF = TRUE)
#  for(num in 1:length(nclus_vec))
#  {
#       StructureObj_omega(Topic_clus_list[[num]]$omega,samp_metadata, batch_lab, path_struct,
#                            partition=partition,
#                            control=control)
#  }}

}
