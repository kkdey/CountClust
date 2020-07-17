#' @title Struture Pie plot using ggplot2
#'
#' @description Make a Pie chart with STRUCTURE grades of membership.
#'
#' @param input_data The input data matrix with samples along the rows
#'   and the columns either represnting either features (for
#'   \code{input_type = apply_tsne/apply_pca}) on which to apply
#'   t-SNE/PCA, or columns representing the 2 co-ordinates to plot in
#'   STRUCTURE pie chart plot.
#' 
#' @param input_type The type of input data provided. One of three
#'   possible options: \code{apply_tsne} : input data on which t-SNE is
#'   applied and the STRUCTURE pie chart co-ordinates defined by first
#'   two dimensions of the t-SNE \code{apply_pca} : input data on which
#'   PCA is applied and the STRUCTURE pie chart co-ordinates are defined
#'   by first two dimensions of the PCA.  \code{coord}: Co-ordinates
#'   data with 2 columns representing the two co-ordinate axes of the
#'   STRUCTURE pie chart.
#' 
#' @param omega Cluster membership probabilities of each
#'   sample. Usually a sample by cluster matrix in the Topic model
#'   output. The cluster weights sum to 1 for each sample.
#' 
#' @param use_voom Use \code{\link[limma]{voom}} to transform the
#'   count data before performing PCA or t-SNE.
#' 
#' @param color_set The set of colors used for the pie charts in
#'   STRUCTURE pie.  Defaults to NULL, which uses a set of 75
#'   qualitatively distinct colors.
#' 
#' @param pie_radius The radius of the pie chart in STRUCTUREpie plot.
#' 
#' @param xlab X-label of the STRUCTURE pie chart.
#' 
#' @param ylab Y-label of the STRUCTURE pie chart.
#' 
#' @param main The title of the STRUCTURE pie chart.
#' 
#' @param control Control paramaters for the STRUCTURE pie chart plot.
#'   User can add padding to the X and Y axes using \code{padding}
#'   option in control, tune legend pie radius and size
#'   (\code{legend_pie_radius} and \code{legend_pie_cex}), the plot
#'   background (\code{bg}), the color intensity of pie charts
#'   (\code{color_intensity}) and legend location (\code{legendx} and
#'   \code{legendy}).
#'
#' @return Plots the Structure Pie chart visualization of the t-SNE or
#'   PCA or user defined co-ordinate decomposition of data.
#' 
#' @examples
#' library(singleCellRNASeqMouseDeng2014)
#' deng.counts <- exprs(Deng2014MouseESC)
#' data("MouseDeng2014.FitGoM")
#' omega <- MouseDeng2014.FitGoM$clust_6$omega
#' set.seed(1000)
#' StructurePie(t(deng.counts), input_type="apply_tsne",
#'             use_voom=FALSE, omega = omega, xlab="TSNE1",
#'             ylab = "TSNE2",
#'             main = "STRUCTURE K=6 pie on tSNE",
#'             control = list(bg = "lightcyan"))
#' StructurePie(t(deng.counts), input_type="apply_pca",
#'             use_voom = TRUE, omega = omega, xlab="PCA1",
#'             ylab = "PCA2",
#'             main = "STRUCTURE K=6 pie on PCA",
#'             control = list(bg = "lightcyan"))
#'
#' @importFrom Rtsne Rtsne
#' @importFrom limma voom
#'
#' @export
#'
StructurePie <- function(input_data,
                        input_type,
                        omega,
                        color_set=NULL,
                        use_voom = TRUE,
                        pie_radius = 0.8,
                        xlab = "Co-ordinate1",
                        ylab = "Co-ordinate2",
                        main = "STRUCTURE pie chart",
                        control = list()){

    control.default <- list(color_intensity = 0.9,
                            bg = "white",
                            padding = c(2,2),
                            legend_pie_radius = 2,
                            legend_pie_cex = 0.8,
                            legendx = NULL,
                            legendy = NULL,
                            edges = 200, clockwise = TRUE,
                            init.angle = 90, density = NULL,
                            angle = 45, border = NA,
                            lty = NULL, label.dist = 1.1)

    namc <- names(control)
    if (!all(namc %in% names(control.default)))
        stop("unknown names in control: ",
             namc[!(namc %in% names(control.default))])
    control <- modifyList(control.default, control)

    pie_control <- list(edges = control$edges,
                       clockwise = control$clockwise,
                       init.angle = control$init.angle,
                       density = control$density,
                       angle = control$angle,
                       border = control$border,
                       lty = control$lty,
                       label.dist = control$label.dist)

    if(any(omega < 0) | any(omega > 1)){
        stop("omega must be values between 0 and 1")
    }
    if(dim(omega)[1] != dim(input_data)[1]){
        stop("Number of rows in omega matrix for grades of membership must
             match the number of rows in the input_data")
    }

    #############   Normalize the rows of the omega matrix ###############
    omega_norm <- t(apply(omega, 1, function(x) return(x/sum(x))))

    ################  Build co-ordinates for the STRUCTURE pie ##############
    if(input_type == "apply_tsne"){
        if(use_voom){
            tsne_out <- Rtsne(t(voom(t(input_data))$E))
        }else{
            tsne_out <- Rtsne(input_data)
        }
        coord1 <- tsne_out$Y[,1]
        coord2 <- tsne_out$Y[,2]
    }else if (input_type == "apply_pca"){
        if(use_voom){
            pca_out <- prcomp(t(voom(t(input_data))$E))
        }else{
            pca_out <- prcomp(input_data)
        }
        coord1 <- pca_out$x[,1]
        coord2 <- pca_out$x[,2]
    }else if (input_type == "coord"){
        if(ncol(input_data) > 2){
            warning("input data has more than 2 columns and the input_type
                    chosen is coord, so the first two columns are only selected
                    for the X and Y axis of the plot. If you want to use t-SNE
                    or PCA to reduce the data to 2 dimensions, please use
                    input_type = apply_tsne or apply_pca. Else provide a 2 column
                    input data for generating the axes of your choice.")
        }
        coord1 <- input_data[,1]
        coord2 <- input_data[,2]
    }else{
         stop("input_type must be apply_tsne, apply_pca or coord")
    }

    if (is.null(control$legendx)){
        control$legendx = max(coord1) + 0.2*control$padding[1]
    }
    if (is.null(control$legendy)){
        control$legendy = max(coord2) + 0.2*control$padding[2]
    }
    range1 <- range(coord1)[2] - range(coord1)[1]
    range2 <- range(coord2)[2] - range(coord2)[1]
    range <- max(range1, range2)
    pie_radius <- pie_radius*10^(max(floor(log(range, 10))-1, 0))
    control$legend_pie_radius <- control$legend_pie_radius*10^(max(0, floor(log(range, 10))-1))

    ##################  Colors to plot in STRUCTURE pie #################
    if(is.null(color_set)){
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual",]
        color_set = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                                  rownames(qual_col_pals)))
        if(length(color_set) < ncol(omega_norm)){
            stop("Enter in color_set a vector of colors larger than the
                 number of clusters")
        }
    }

    color_with_intensity <- c()
    for(dim in 1:ncol(omega_norm)){
        color_with_intensity <- c(color_with_intensity,
                                 alpha(color_set[dim],
                                       control$color_intensity))
    }

    basemap2(xlim=c(min(coord1) - control$padding[1],
                    max(coord1) + control$padding[1]),
             ylim=c(min(coord2) - control$padding[2],
                    max(coord2) + control$padding[2]),
             main = main, xlab = xlab, ylab = ylab, bg=control$bg)
    do.call(draw.pie,
            append(list(
                x = coord1, y = coord2, z = 100*omega_norm,
                radius = pie_radius, col = color_with_intensity),
                pie_control))
    legend.pie(control$legendx, control$legendy,
               labels=paste0("c", 1:ncol(omega_norm)),
               radius=control$legend_pie_radius, bty="n",
               col=color_set[1:ncol(omega_norm)], cex=control$legend_pie_cex,
               label.dist=1.3)
}

basemap2 <- function(xlim,ylim,xlab="Longitude",ylab="Latitude",
                     bg="lightblue",...){
    # more-or-less appropriateaspect ratio
    plot(NA,NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    # background colour (does not work properly if you re-size plot window)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
}
