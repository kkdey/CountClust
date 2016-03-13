#' Struture plot with ggplot2 package
#'
#' Make the traditional Structure histogram plot of GoM model using ggplot2
#'
#' @param omega Cluster membership probabilities of each sample. Usually a 
#'              sample by cluster matrix in the Topic model output. The 
#'              cluster weights sum to 1 for each sample.
#' @param annotation A data.frame of two columns: sample_id and tissue_label. 
#'                    sample_id is the unique identifying number of each 
#'                    sample (alphanumeric). tissue_lable is a factor of 
#'                    tissue labels, with levels of the factor arranged 
#'                    in the order of the tissues in the Structure 
#'                    (left to right).
#' @param palette A vector of colors assigned to the clusters. First color 
#'                in the vector is assigned to the cluster labeled as 1, 
#'                and second color in the vector is assigned to the cluster 
#'                labeled as 2, etc. The number of colors must be the same 
#'                or greater than the number of clusters. The clusters not 
#'                assigned a color are filled with white in the figure. 
#'                In addition, the recommended
#'                choice of color palette here is RColorBrewer,
#'                for instance RColorBrewer::brewer.pal(8, "Accent") or
#'                RColorBrewwer::brewer.pal(9, "Set1").
#' @param figure_title Title of the plot.
#' @param yaxis_label Axis label for the samples.
#' @param order_sample if TRUE, we order samples in each annotation batch sorted by membership
#'                     of most representative cluster. If FALSE, we keep the order in the data.
#' @param sample_order_decreasing if order_sample TRUE, then this input determines if the ordering
#'                                due to main cluster is in ascending or descending order.
#' @param split_line Control parameters for line splitting the batches in the plot.
#' @param axis_tick Control parameters for x-axis and y-axis tick sizes.
#'
#' @return A ggplot-version of the GoM model visualization
#'
#' @import cowplot
#' @import reshape2
#' @import plyr
#' @import RColorBrewer
#' @export
#'
#' @examples
#' # load the previously analyzed results
#' data("MouseDeng2014.FitGoM")
#'
#' # extract the omega matrix: membership weights of each cell
#' names(MouseDeng2014.FitGoM$clust_6)
#' omega <- MouseDeng2014.FitGoM$clust_6$omega
#'
#' # make annotation matrix
#' annotation <- data.frame(
#' sample_id = paste0("X", c(1:NROW(omega))),
#' tissue_label = factor(rownames(omega),
#'                      levels = rev( c("zy", "early2cell",
#'                                      "mid2cell", "late2cell",
#'                                      "4cell", "8cell", "16cell",
#'                                      "earlyblast","midblast",
#'                                      "lateblast") ) ) )
#' rownames(omega) <- annotation$sample_id;
#' StructureGGplot(omega = omega,
#'                  annotation = annotation,
#'                  palette = RColorBrewer::brewer.pal(8, "Accent"),
#'                  yaxis_label = "development phase",
#'                  order_sample = TRUE,
#'                  axis_tick = list(axis_ticks_length = .1,
#'                                   axis_ticks_lwd_y = .1,
#'                                   axis_ticks_lwd_x = .1,
#'                                   axis_label_size = 7,
#'                                   axis_label_face = "bold"))
#'




StructureGGplot <- function(omega, annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            sample_order_decreasing = TRUE,
                            split_line = list(split_lwd = 1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 3,
                                             axis_label_face = "bold") ) {


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

    df_ord <- do.call(rbind,
                      lapply(1:nlevels(annotation$tissue_label), function(ii) {
                          temp_label <- levels(annotation$tissue_label)[ii]
                          temp_df <- omega[which(annotation$tissue_label == temp_label), ]

                          # find the dominant cluster in each sample
                          each_sample_order <- apply(temp_df, 1, which.max)

                          # find the dominant cluster across samples
                          sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])

                          if (order_sample == TRUE) {
                          # reorder the matrix
                              temp_df_ord <- temp_df[order(temp_df[ , sample_order],
                                                           decreasing = sample_order_decreasing), ]
                          } else {
                              temp_df_ord <- temp_df
                          }
                          temp_df_ord
                      }) )

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

    # set axis tick positions
    tissue_count <- c(0, cumsum(table(droplevels(annotation$tissue_label)) ) )
    tissue_names <- levels(droplevels(annotation$tissue_label))
    tissue_breaks <- sapply(1:(length(tissue_count)-1), function(i) {
        round((tissue_count[i] + tissue_count[i+1])/2)
    })
    names(tissue_breaks) <- tissue_names
    #cbind(tissue_breaks, cumsum(table(annotation$tissue_label)))

    # make ggplot
    a <- ggplot(df_mlt,
                aes(x = df_mlt$document, y = df_mlt$value*10000, fill = factor(df_mlt$topic)) ) +
        xlab(yaxis_label) + ylab("") +
        scale_fill_manual(values = palette) +
        theme(legend.position = "right",
              legend.key.size = unit(.2, "cm"),
              legend.text = element_text(size = 5),
##<-- TBD: center legend title
#              legend.title = element_text(hjust = 1),
              axis.text = element_text(size = axis_tick$axis_label_size,
                                       face = axis_tick$axis_label_face),
              axis.ticks.y = element_line(size = axis_tick$axis_ticks_lwd_y),
              axis.ticks.x = element_line(size = axis_tick$axis_ticks_lwd_x),
              axis.ticks.length = unit(axis_tick$axis_ticks_length, "cm"),
              title = element_text(size = 6) ) +
        ggtitle(figure_title) +
        scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                            labels = seq(0, 1, 1/(ticks_number -1 ) ) ) +
        # Add tissue axis labels
        scale_x_discrete(breaks = levels(df_mlt$document)[tissue_breaks],
                         labels = names(tissue_breaks)) +
        # Add legend title
        labs(fill = "Clusters") +
        coord_flip()


    # width = 1: increase bar width and in turn remove space
    # between bars
    b <- a + geom_bar(stat = "identity",
                      position = "stack",
                      width = 1)
    b <- b + cowplot::panel_border(remove = TRUE)
    # Add demarcation (TBI)
    b <- b + geom_vline(
        xintercept = cumsum(table(droplevels(annotation$tissue_label)))[
            -length(table(droplevels(annotation$tissue_label)))],
        col = split_line$split_col,
        size = split_line$split_lwd)
    b <- cowplot::ggdraw(cowplot::switch_axis_position((b), axis = "y"))
    b
    # tidy up the figure output using cowplot::plot_grid (optoinal)
#    cowplot::plot_grid(b)
}




