#' @title Topic model fitting and Structure Plot!
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param nclus the number of clusters or topics to be fitted.
#' @param samp_metadata the sample metadata, samples along the rows and each column representing some metadata information 
#'        that will be used to arrange the Structure plot columns (one plot for one arrangement).
#' @param tol the tolerance value for topic model fit (set to 0.001 as default)
#' @param batch_lab the batch labels, the output will have one Structure plot arranged by batch labels too.
#' 
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



StructureObj <- function(data, nclus, samp_metadata, tol, batch_lab)
{
  Topic_clus <- topics(data, K=nclus, tol=tol);
  docweights <- Topic_clus$omega;
  message('Fitting the topic model (due to Matt Taddy)', domain = NULL, appendLF = TRUE)
  write.table(Topic_clus$omega,'Structure/topic_data/omega_mat.txt');
  write.table(Topic_clus$theta,'Structure/topic_data/theta_mat.txt');
  dir.create(paste0('Structure/plots/clus_',nclus));
  num_metadata <- dim(samp_metadata)[2];
  
  message('Creating the Structure plots', domain = NULL, appendLF = TRUE)
  
  for(num in 1:num_metadata)
  {
    metadata_vec <- samp_metadata[,num];
    metadata_ordered <- metadata_vec[order(metadata_vec)];
    docweights_ordered <- docweights[order(metadata_vec),];
    png(filename=paste0('Structure/plots/clus_',nclus,'/struct_',colnames(samp_metadata)[num],'.png'));
    barplot(t(docweights_ordered),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("Structure arranged by",colnames(samp_metadata)[num],
                                                                                       ": topics=",nclus),las=1,ylim=c(0,1),cex.axis=0.3,cex.main=1.4);
    labels = match(unique(metadata_ordered), metadata_ordered);
    abline(v=labels-1)
    
    labels_low=labels-1;
    labels_up=c(labels_low[2:length(labels_low)],dim(docweights_ordered)[1]);
    mid_point <- labels_low +0.5*(labels_up-labels_low);
    axis(1,at=mid_point, unique(metadata_ordered),las=2,cex.axis=0.3);
    dev.off()
    
  }
  
  if(!is.null(batch_labs)){
    batch_vec <- batch_labs;
    batch_vec_ordered <- batch_vec[order(batch_vec)];
    docweights_ordered <- docweights[order(batch_vec),];
    png(filename=paste0('Structure/plots/clus_',nclus,'/struct_batch','.png'));
    barplot(t(docweights_ordered),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("Structure arranged by batch",
                                                                                       ": topics=",nclus),las=1,ylim=c(0,1),cex.axis=0.3,cex.main=1.4);
    labels = match(unique(batch_vec_ordered), batch_vec_ordered);
    abline(v=labels-1)
    
    labels_low=labels-1;
    labels_up=c(labels_low[2:length(labels_low)],dim(docweights_ordered)[1]);
    mid_point <- labels_low +0.5*(labels_up-labels_low);
    axis(1,at=mid_point, unique(batch_vec_ordered),las=2,cex.axis=0.3);
    dev.off()
  }

  ll <- list("omega"=Topic_clus$omega, "theta"=Topic_clus$theta, "bf"=Topic_clus$BF)
  return(ll)
}