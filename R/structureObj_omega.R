#' @title Structure Plot given topic model proportions!
#'
#' @param omega the topic proportion matrix returned by topics() in maptpx,
#'        or the omega component from structureObj, with samples along the rows and topics along columns.
#' @param samp_metadata the sample metadata, samples along the rows and each column representing some metadata information
#'        that will be used to arrange the Structure plot columns (one plot for one arrangement).
#' @param batch_lab the batch labels, the output will have one Structure plot arranged by batch labels too.
#' @param path The directory path where we want to save the data and Structure plots.
#' @param partition A logical vector of same length as metadata. partition[i]=TRUE will imply that for the Structure
#'            plot for i th metadata, no vertical line parititon between classes is used.
#' @param control() A list of control parameters for the Structure plot. The control list has the arguments
#'        struct.width, struct.height, cex.axis, cex.main, lwd, las and color and margin parameters.
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


StructureObj_omega <- function(omega, samp_metadata, tol, batch_lab, path,
                         partition=rep('TRUE',ncol(samp_metadata)),
                         control=list())
{
  control.default <- list(struct.width=600, struct.height=400, cex.axis=0.5, cex.main=1.5, las=2, lwd=2,
                          mar.bottom =14, mar.left=2, mar.top=2, mar.right=2,color=2:(dim(omega)[2]+1));

  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=modifyList(control.default, control)

  struct.width <- control$struct.width;
  struct.height <- control$struct.height;
  cex.axis <- control$cex.axis;
  cex.min <- control$cex.main;
  las <- control$las;
  lwd <- control$lwd;
  mar.bottom <- control$mar.bottom;
  mar.left <- control$mar.left;
  mar.top <- control$mar.top;
  mar.right <- control$mar.right;
  color <- control$color;

  docweights <- as.matrix(omega);
  nclus <- dim(docweights)[2];
  num_metadata <- dim(samp_metadata)[2];

  message('Creating the Structure plots', domain = NULL, appendLF = TRUE)

  for(num in 1:num_metadata)
  {
    metadata_vec <- samp_metadata[,num];
    metadata_ordered <- metadata_vec[order(metadata_vec)];
    docweights_ordered <- docweights[order(metadata_vec),];
    png(filename=paste0(path,'/struct_clus_',nclus,'_',colnames(samp_metadata)[num],'.png'),width=struct.width, height=struct.height);
    par(mar=c(mar.bottom,mar.left, mar.top,mar.right))
    barplot(t(docweights_ordered),col=color,axisnames=F,space=0,border=NA,
            main=paste("Structure arranged by",colnames(samp_metadata)[num],": topics=",nclus),
            las=las,ylim=c(0,1),ylab="admix prop", xlab=paste0(colnames(samp_metadata)[num]),
            cex.axis=cex.axis,cex.main=cex.main);
    labels = match(unique(metadata_ordered), metadata_ordered);
    if(partition[num]=='TRUE') abline(v=labels-1, lty=1, lwd=2)

    labels_low=labels-1;
    labels_up=c(labels_low[2:length(labels_low)],dim(docweights_ordered)[1]);
    mid_point <- labels_low +0.5*(labels_up-labels_low);
    axis(1,at=mid_point, unique(metadata_ordered),las=las,cex.axis=cex.axis,lwd=lwd);
    dev.off()

  }

  if(!is.null(batch_lab)){
    batch_vec <- batch_lab;
    batch_vec_ordered <- batch_vec[order(batch_vec)];
    docweights_ordered <- docweights[order(batch_vec),];
    png(filename=paste0(path,'/struct_clus_',nclus,'_batch.png'),width=struct.width, height=struct.height);
    par(mar=c(mar.bottom,mar.left, mar.top,mar.right))
    barplot(t(docweights_ordered),col=color,axisnames=F,space=0,border=NA,
            main=paste("Structure arranged by batch",": topics=",nclus),
            las=las,ylim=c(0,1),ylab="admix prop", xlab="batch",
            cex.axis=cex.axis,cex.main=cex.main);
    labels = match(unique(batch_vec_ordered), batch_vec_ordered);
    if(partition[num]=='TRUE')  abline(v=labels-1, lty=1, lwd=2)

    labels_low=labels-1;
    labels_up=c(labels_low[2:length(labels_low)],dim(docweights_ordered)[1]);
    mid_point <- labels_low +0.5*(labels_up-labels_low);
    axis(1,at=mid_point, unique(batch_vec_ordered),las=las,cex.axis=cex.axis,lwd=lwd);
    dev.off()
  }
}
