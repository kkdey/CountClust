#' Apply t-SNE on the topic proportions from structure
#'
#' @param  omega the topic proportion matrix obtained from StructureObj/ topics ()  in maptpx
#' @param  samp_metadata The sample labels used for identifying labels in tSNE plot
#' @param  path  The path where you want to save the t-SNE plot
#' 
#' @description This function transforms the simplex matrix omega to unconstrained matrix and then applies t-SNE (Van der Maaten) and plots the results on 2D
#' 
#' @references L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research, 2014:3221-3245
#' 
#'            L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 2008: 2579-2605.
#'
#'            Broman KW. R/qtlcharts: interactive graphics for quantitative trait locus mapping. 2015. Genetics 199:359-361


tsne_struct <- function(omega,samp_metadata,path)
{
  reverse_transform=function(x) 
  {
    return(log((abs(x[2:length(x)])+1e-7)/(abs(x[1])+1e-7)));
  }
  
  transform <- function(y) 
  {
    temp =c(1,exp(y));
    out=temp/sum(temp);
    return(out)
  }
  
  rev_omega <- apply(omega,1,function(x) reverse_transform(x));
  tsne_out <- tsne(rev_omega,2);
  if(is.null(samp_metadata)){
    img <- iplot(tsne_samples[,1],tsne_samples[,2]);
    htmlwidgets::saveWidget(img, file=paste0(path,"/tSNE_omega.html"),selfcontained=FALSE);
  }
  if(!is.null(samp_metadata)){
    iplot(tsne_samples[,1],tsne_samples[,2],rep(1,length(tsne_samples[,1])),as.factor(samp_metadata),samp_metadata)
    htmlwidgets::saveWidget(img, file=paste0(path,"/tSNE_omega.html"),selfcontained=FALSE);
  }
  return(tsne_out)
}