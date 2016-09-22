#' @title Grade of Membership (GoM) model fit !
#'
#' @description Fits a grade of membership model to count data. Default input
#'                includes a sample-by-feature matrix, the number of clusters
#'                (topics) to fit (K). The function is a wrapper of the topics()
#'                function implemented in Matt Taddy's maptpx pacakge.
#'
#' @param data counts data, with samples along the rows and features
#'              along the columns.
#' @param K the vector of clusters or topics to be fitted.
#' @param tol Tolerance value for GoM model log posterior increase
#'            at successive iterations (set to 0.1 as default).
#' @param path_rda The directory path for saving the GoM model output.
#'                  If NULL, it will return the output to console.
#' @param control Control parameters. Same as topics() function of
#'                 maptpx package.
#'
#' @return Saves the GoM model fit output for each cluster in vector K at the
#'                directory path in path_rda.
#'
#' @references Matt Taddy. On Estimation and Selection for Topic Models.
#'                AISTATS 2012, JMLR W\&CP 22.
#'
#'            Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly.
#'                Inference of population structure using multilocus genotype
#'                data. Genetics 155.2 (2000): 945-959.
#'
#' @keywords counts data, clustering, Structure plot
#'
#'
#' @examples
#'
#' data("ex.counts")
#' out <- FitGoM(ex.counts, K=4, tol=100, control=list(tmax=100))
#'
#' @importFrom maptpx topics
#' @import slam
#' @importFrom utils modifyList
#' @export


FitGoM <- function(data,
                   K,
                   tol=0.1,
                   path_rda = NULL,
                   control=list())
{
    ## dealing with blank rows: we first remove them

    control.default <- list(shape=NULL, initopics=NULL, bf=TRUE,
                            kill=2, ord=TRUE, verb=1, tmax=1000)

    namc=names(control)
    if (!all(namc %in% names(control.default)))
        stop("unknown names in control: ",
             namc[!(namc %in% names(control.default))])
    control=modifyList(control.default, control)



    indices_blank <- as.numeric(which(apply(data,1,max) == 0))
    if(length(indices_blank)!=0){
        data <- as.matrix(data[-indices_blank,]);
    }

    message('Fitting the topic model (due to Matt Taddy)',
            domain = NULL, appendLF = TRUE)

    Topic_clus_list <- lapply(K, function(per_clust) {

        suppressWarnings(out <- maptpx::topics(
            as.matrix(data),
            K = per_clust,
            shape=control$shape,
            initopics=control$initopics,
            tol=tol,
            bf=control$bf,
            kill=control$kill,
            ord=control$ord,
            verb=control$verb,
            tmax=control$tmax))
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
    #       StructureObj_omega(Topic_clus_list[[num]]$omega,samp_metadata,
    #                          batch_lab, path_struct,
    #                          partition=partition,
    #                          control=control)
    #  }}

}
