#' @title Grade of Membership (GoM) model fit !
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param K the vector of clusters or topics to be fitted.
#' @param tol the tolerance value for GoM model log posterior increase at successive iterations (set to 0.1 as default)
#' @param path_rda The directory path for saving the GoM model output. If NULL, it will return the output to console.
#' @param control The control parameters as in topics() function of maptpx package.
#'
#' @description This function takes the counts data (no. of samples x no. of features) and the value of K, the number of topics or
#' cluster to fit, fits the GoM model (using topics() function of maptpx package) and saves the model fit putput
#'
#' @return Saves the GoM model fit output for each cluster in vector K at the directory path in path_rda.
#'
#' @references Matt Taddy.On Estimation and Selection for Topic Models. AISTATS 2012, JMLR W\&CP 22.
#'
#'            Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly. Inference of population structure using multilocus genotype data.
#'            Genetics 155.2 (2000): 945-959.
#' @keywords counts data, clustering, Structure plot
#'
#' @export
#'
#' @examples
#'
#' library(GTExV6Brain)
#' gtex.counts <- exprs(GTExV6Brain)
#' gtex.meta_data <- pData(GTExV6Brain)
#' gtex.gene_names <- rownames(gtex.counts)
#' ex.counts <- t(gtex.counts[1:200,])
#' FitGoM(ex.counts,
#'            K=4, tol=1000, control=list(kill=1))



FitGoM <- function(data, K, tol, path_rda=NULL, control=list())
{
  ## dealing with blank rows: we first remove them

  control.default <- list(shape=NULL, initopics=NULL, bf=FALSE,
                          kill=2, ord=TRUE, verb=1)

  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=modifyList(control.default, control)



  indices_blank <- as.numeric(which(apply(data,1,max)==0));
  if(length(indices_blank)!=0){
  data <- as.matrix(data[-indices_blank,]);
  }

  message('Fitting the topic model (due to Matt Taddy)', domain = NULL, appendLF = TRUE)

  Topic_clus_list <- lapply(K, function(per_clust) {
    suppressWarnings(out <- maptpx::topics(as.matrix(data), K = per_clust, control$shape, control$initopics, tol,
                                           control$bf, control$kill, control$ord, control$verb))
    return(out)
  })

  names(Topic_clus_list) <- paste0("clust_",K)
  if(!is.null(path_rda)){
  save(Topic_clus_list, file = path_rda);
  }else{
    return(Topic_clus_list)
  }
#  if(plot) {
#  message('Creating the Structure plots', domain = NULL, appendLF = TRUE)
#  for(num in 1:length(nclus_vec))
#  {
#       StructureObj_omega(Topic_clus_list[[num]]$omega,samp_metadata, batch_lab, path_struct,
#                            partition=partition,
#                            control=control)
#  }}

}
