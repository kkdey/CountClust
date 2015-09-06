#' @title The CountsClust function- a function for clustering and visualization of counts!
#'
#' @param data counts data, with samples along the rows and features along the columns.
#' @param nclus_vec a vector of topics/clusters that user wants to fit.
#' @param samp_metadata the sample metadata, samples along the rows and each column representing some metadata information 
#' @param tol the tolerance value for topic model fit (set to 0.001 as default)
#' @param batch_lab the batch labels, the output will have one Structure plot arranged by batch labels too.
#' @param use_parallel equal to TRUE if parallel implementation / default is FALSE
#' @param use_tsne equal to TRUE if the user wants to compute t-SNE / default is TRUE
#' @param use_pca equal to TRUE is the user wants to compute the PCA/ default is TRUE
#' @param thresh_prop The threshold proportion used to remove/substitute the columns of data containing NA. Deafult is 0 implying that all columns 
#'        with NA are removed.
#' @param filter_prop The threshold proportion that determines if the feature is too sparse or not. Default is 0.9.
#' @param expand_factor An expansion factor used to obtain batch corrected counts
#' 
#' @description This function takes the counts data (no. of samples x no. of features), the vector of topics/clusters the user wants to fit, along
#' with the sample metadata and batch label information and it produces the Structure plots with and without controlling fro batch effects for 
#' different number of topics, with arrangement determined by the sample metadata. Additionally, it outputs the cluster data output and model Bayes 
#' factor, and also plots t-SNE and the PCA for both with and without batch cluster outputs.
#' 
#' @author Kushal K Dey, Matthew Stephens
#' 
#' @references Matt Taddy.On Estimation and Selection for Topic Models. AISTATS 2012, JMLR W\&CP 22.
#' 
#'             Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly. Inference of population structure using multilocus genotype data. 
#'             Genetics 155.2 (2000): 945-959.
#'             
#'             Broman KW. R/qtlcharts: interactive graphics for quantitative trait locus mapping. 2015. Genetics 199:359-361
#'             
#'             L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research, 2014:3221-3245
#'            
#' @keywords counts data, clustering, topic model, Structure plot
#' 
#' @export
#' 

countsclust <- function(data, 
                        nclus_vec, 
                        samp_metadata,
                        samp_lab_tsne,
                        tol=0.001, 
                        batch_lab, 
                        use_parallel=FALSE, 
                        use_tsne=TRUE, 
                        use_pca=TRUE, 
                        thresh_prop=0, 
                        filter_prop=0.9,
                        expand_factor=100)
{
  data_preprocessed <- handleNA(data,thresh_prop=thresh_prop)$data;
  data_filtered <- RemoveSparseFeatures(data_preprocessed,filter_prop = filter_prop)$data;
  
  counts <- data_filtered;
  rm(data_filtered)
  
  if(!dir.exists("Structure")) dir.create("Structure")
  if(!dir.exists("Structure/topic_data")) dir.create("Structure/topic_data")
  if(!dir.exists("Structure/plots")) dir.create("Structure/plots")
  
  bayesfac_batch_uncorrected <- array(0,length(nclus_vec));
  tsne_out <- vector("list", length(nclus_vec)); 
  
  if(!dir.exists("tSNE")) dir.create("tSNE")
  if(!dir.exists("tSNE/batch_uncorrected")) dir.create("tSNE/batch_uncorrected")
  
  for(num in 1:length(nclus_vec))
  {
    obj <- StructureObj(counts,nclus_vec[num],samp_metadata = samp_metadata, tol=tol, batch_lab = batch_lab);
    bayesfac[num] <- obj$bf;
    if(!dir.exists(paste0("tSNE/batch_uncorrected/clus_",nclus_vec[num]))) dir.create(paste0("tSNE/batch_uncorrected/clus_",nclus_vec[num]))
    tsne_out[[num]] <- tsne_struct(obj$omega,samp_lab_tsne,paste0("tSNE/batch_uncorrected/clus_",nclus_vec[num]));
  }
  
  batch_corrected_counts <- BatchCorrectedCounts(counts,expand_factor = expand_factor,batch_lab)
  if(!dir.exists("tSNE/batch_corrected")) dir.create("tSNE/batch_uncorrected")
  
  bayesfac_batchcorrect <- array(0,length(nclus_vec));
  tsne_out_batchcorrect <- vector("list", length(nclus_vec)); 
  
  for(num in 1:length(nclus_vec))
  {
    obj <- StructureObj(batch_corrected_counts,nclus_vec[num],samp_metadata = samp_metadata, tol=tol, batch_lab = batch_lab);
    bayesfac_batchcorrect[num] <- obj$bf;
    if(!dir.exists(paste0("tSNE/batch_corrected/clus_",nclus_vec[num]))) dir.create(paste0("tSNE/batch_corrected/clus_",nclus_vec[num]))
    tsne_out_batchcorrect[[num]] <- tsne_struct(obj$omega,samp_lab_tsne,paste0("tSNE/batch_corrected/clus_",nclus_vec[num]));
  }
  
  
  
}
  