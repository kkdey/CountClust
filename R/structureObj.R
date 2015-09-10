#' @title Topic model fitting and Structure Plot!
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param nclus the number of clusters or topics to be fitted.
#' @param samp_metadata the sample metadata, samples along the rows and each column representing some metadata information
#'        that will be used to arrange the Structure plot columns (one plot for one arrangement).
#' @param tol the tolerance value for topic model fit (set to 0.001 as default)
#' @param batch_lab the batch labels, the output will have one Structure plot arranged by batch labels too.
#' @param path The directory path where we want to save the data and Structure plots.
#' @param partition A logical vector of same length as metadata. partition[i]=TRUE will imply that for the Structure
#'            plot for i th metadata, no vertical line parititon between classes is used.
#' @param control() A list of control parameters for the Structure plot. The control list has the arguments
#'        struct.width, struct.height, cex.axis, cex.main, las.struct, lwd, las and margin parameters.
#'
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



StructureObj <- function(data, nclus, samp_metadata, tol, batch_lab, path,
                         partition=rep('TRUE',ncol(samp_metadata)),
                         control=list())
{
  control.default <- list(struct.width=800, struct.height=250, cex.axis=0.5, cex.main=1.5, las=2, lwd=2,
                          mar.bottom =14, mar.left=2, mar.top=2, mar.right=2);

   namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=modifyList(control.default, control)

  struct.width <- control$struct.width;
  struct.height <- control$struct.height;
  cex.axis <- control$cex.axis;
  cex.main <- control$cex.main;
  las <- control$las;
  lwd <- control$lwd;
  mar.bottom <- control$mar.bottom;
  mar.left <- control$mar.left;
  mar.top <- control$mar.top;
  mar.right <- control$mar.right;

  ## dealing with blank rows: we first remove them

  indices_blank <- as.numeric(which(apply(data,1,max)==0));
  if(length(indices_blank)!=0){
  data <- as.matrix(data[-indices_blank,]);
  samp_metadata <- as.matrix(samp_metadata[-indices_blank,]);
  batch_lab <- as.vector(batch_lab[-indices_blank]);
  }

  message('Fitting the topic model (due to Matt Taddy)', domain = NULL, appendLF = TRUE)
  Topic_clus <- topics(data, K=nclus, tol=tol);
  docweights <- Topic_clus$omega;
  write.table(Topic_clus$omega,paste0(path,'/omega_mat.txt'));
  write.table(Topic_clus$theta,paste0(path,'/theta_mat.txt'));
  num_metadata <- dim(samp_metadata)[2];

  message('Creating the Structure plots', domain = NULL, appendLF = TRUE)

  StructureObj_omega(docweights,samp_metadata, tol, batch_lab, path,
                            partition=rep('TRUE',ncol(samp_metadata)),
                            control=control)

  ll <- list("omega"=Topic_clus$omega, "theta"=Topic_clus$theta, "bf"=Topic_clus$BF, "blank_indices"=indices_blank)
  return(ll)
}
