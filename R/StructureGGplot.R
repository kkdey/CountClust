#' Struture plot with ggplot2 package
#' 
#' Make the traditional Structure histogram plot using the ggplot2 pacakges
#' 
#' @param omega Cluster membership probabilities of each sample. Usually a sample by 
#'              cluster matrix in the Topic model output. The cluster weights sum to 1
#'              for each sample.
#' @param annotation A data.frame of two columns: sample_id and tissue_label. sample_id
#'                  is the unique identifying number of each sample (alphanumeric).
#'                  tissue_lable is a factor of tissue labels, with levels of the factor
#'                  arranged in the order of the tissues in the Structure (left to right).
#' @param palette A vector of colors assigned to the clusters. First color in the vector
#'                is assigned to the cluster labeled as 1, and second color in the vector
#'                is assigned to the cluster labeled as 2, etc. The number of colors must be
#'                the same or greater than the number of clusters. The clusters not assigned
#'                a color are filled with white in the figure. In addition, the recommended
#'                choice of color palette here is RColorBrewer, 
#'                for instance RColorBrewer::brewer.pal(8, "Accent") or 
#'                RColorBrewwer::brewer.pal(9, "Set1").
#' @param figure_title Title of the plot.
#' @param yaxis_label Axis label for the samples.
#' 
#' @author Chiaowen Joyce Hsiao and Kushal K Dey
#' 
#' @export
#' 
#' @examples 
#' ## load the previously analyzed results
#' # load("../../../count-clustering/project/rdas/deng_topic_fit.rda")
#' 
#' ## extract the omega matrix: membership weights of each cell
#' # names(Topic_clus_list)
#' # str(Topic_clus_list$clust_6)
#' # omega <- Topic_clus_list$clust_6$omega
#' 
#' ## make annotation matrix
#' # annotation <- data.frame(
#' #   sample_id = c(1:NROW(omega)),
#' #   tissue_label = factor(rownames(omega),
#' #                         levels = c("zy", "early2cell", "mid2cell", "late2cell",
#' #                                     "4cell", "8cell", "16cell", "earlyblast",
#' #                                     "midblast", "lateblast")) )
#' # StructureGGplot(omega = omega, 
#' #                 tissue_labels = cell_labels)                                    
#' #                 
#' # save_plot("../../../count-clustering/project/plots/deng-figures/deng-ggplot.png", 
#' #         b, base_height = 4, base_width = 2)

StructureGGplot <- function(omega, annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Tissue type",
                            order_sample = TRUE) {
    
    library(ggplot2)
    library(reshape2)
    library(cowplot)
    
    # check if the number of colors is same as or more than the number of clusters
    if (dim(omega)[2] > length(palette)) {
        stop("Color choices is smaller than the number of clusters!")
    }
    
    # check if rownames of omega are unique
    if(length(unique(rownames(omega))) != NROW(omega)) {
        stop("omega rownames are not unique!")
    }
    # check the annotation data.frame
    if (!is.data.frame(annotation)) stop("annotation must be a data.framer object")
    if (!all.equal(colnames(annotation), c("sample_id", "tissue_label")) ) {
        stop("annotation data.frame column names must be sample_id and tissue_label")
    }
    if ( length(unique(annotation$sample_id)) != NROW(omega)) {
        stop("sample_id is not unique")
    }
    
    if (order_sample == TRUE) {
    # make the re-ordered dataframe
    df_ord <- do.call(rbind,
#                      lapply(1:20, function(ii) {
                      lapply(1:nlevels(annotation$tissue_label), function(ii) {
                          temp_label <- levels(annotation$tissue_label)[ii]
                          temp_df <- omega[which(annotation$tissue_label == temp_label), ]
                          
                          # find the dominant cluster in each sample
                          each_sample_order <- apply(temp_df, 1, which.max)
                          
                          # find the dominant cluster across samples
                          sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
                          
                          # reorder the matrix
                          temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
                          
                          temp_df_ord
                      }) )
    } else {
        df_ord <- omega
    }

    df_mlt <- reshape2::melt(t(df_ord))
    df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                               "Var2" = "document"))
    df_mlt$document <- factor(df_mlt$document)
    df_mlt$topic <- factor(df_mlt$topic)

    # set blank background
    theme_set(theme_bw(base_size = 12)) +
        theme_update( panel.grid.minor.x = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_blank() )
    
    # inflat nubmers to avoid rounding errors
    value_ifl <- 10000
    
    # number of ticks for the weight axis, including 0 and 1
    ticks_number <- 6

    # make ggplot
    a <- ggplot(df_mlt, 
                aes(x = document, y = value*10000, fill = factor(topic)) ) + 
        xlab(yaxis_label) + ylab("") +
        scale_fill_manual(values = palette) +
        theme(legend.position = "none",
              axis.text = element_text(size = 4),
              title = element_text(size = 6)) +
        ggtitle(figure_title) + 
        scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                            labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
        coord_flip() 
    
    # width = 1: increase bar width and in turn remove space
    # between bars
    b <- a + geom_bar(stat = "identity", 
                      position = "stack", 
                      width = 1)
    b <- b + panel_border(remove = TRUE)
    b <- ggdraw(switch_axis_position((b), axis = "y"))
    
    # tidy up the figure output using cowplot::plot_grid (optoinal)
    cowplot::plot_grid(b)
}




