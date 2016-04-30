#' Struture plot with ggplot2 package
#'
#' Make the traditional Structure histogram plot of GoM model using ggplot2
#'
#' @param omega Cluster membership probabilities of each sample. Usually a
#' sample by cluster matrix in the Topic model output. The cluster weights
#' sum to 1 for each sample.
#' @param annotation A data.frame of two columns: sample_id and tissue_label.
#' sample_id is the unique identifying number of each sample (alphanumeric).
#' tissue_lable is a factor of tissue labels, with levels of the factor
#' arranged in the order of the tissues in the Structure (left to right).
#' @param palette A vector of colors assigned to the clusters. First color in
#' the vector is assigned to the cluster labeled as 1, and second color in the
#' vector is assigned to the cluster labeled as 2, etc. The number of colors
#' must be the same or greater than the number of clusters. The clusters not
#' assigned a color are filled with white in the figure. In addition, the
#' recommended choice of color palette here is RColorBrewer, for instance
#' RColorBrewer::brewer.pal(8, "Accent") or RColorBrewwer::brewer.pal(9, "Set1").
#' @param figure_title Title of the plot.
#' @param yaxis_label Axis label for the samples.
#' @param order_sample if TRUE, we order samples in each annotation batch
#' sorted by membership of most representative cluster. If FALSE, we keep
#' the order in the data.
#' @param sample_order_decreasing if order_sample TRUE, then this input
#' determines if the ordering due to main cluster is in ascending or descending
#' order.
#' @param split_line Control parameters for line splitting the batches in the
#' plot.
#' @param axis_tick Control parameters for x-axis and y-axis tick sizes.
#' @param plot_labels A logical parameter, if TRUE the function plots the axis labels.
#'
#' @return Plots the Structure plot visualization of the GoM model
#'
#' @examples
#' # Example 1
#' data("MouseDeng2014.FitGoM")
#'
#' # extract the omega matrix: membership weights of each cell
#' names(MouseDeng2014.FitGoM$clust_6)
#' omega <- MouseDeng2014.FitGoM$clust_6$omega
#'
#' # make annotation matrix
#' annotation <- data.frame(
#'   sample_id = paste0("X", c(1:NROW(omega))),
#'   tissue_label = factor(rownames(omega),
#'                      levels = rev( c("zy", "early2cell",
#'                                      "mid2cell", "late2cell",
#'                                      "4cell", "8cell", "16cell",
#'                                      "earlyblast","midblast",
#'                                      "lateblast") ) ) )
#' head(annotation)
#' rownames(omega) <- annotation$sample_id
#'
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
#' # Example 2
#' # Import Deng et al data
#'
#' # function to read Deng data from GitHub
#' read.data <- function() {
#'   x <- tempfile()
#'   download.file(paste0("https://cdn.rawgit.com/kkdey/",
#'                          "singleCellRNASeqMouseDeng2014",
#'                          "/master/data/Deng2014MouseEsc.rda"),
#'                 destfile = x, quiet = TRUE)
#'   z <- get(load((x)))
#'   return(z)
#'   }
#' Deng2014MouseESC <-read.data()
#'
#' deng.counts <- Biobase::exprs(Deng2014MouseESC)
#' deng.meta_data <- Biobase::pData(Deng2014MouseESC)
#' deng.gene_names <- rownames(deng.counts)
#'
#' samples_subvector <- which(!duplicated(deng.meta_data$cell_type))[1:3]
#'
#' # Fit GoM on 3 samples with K = 3
#' fit_k3 <- FitGoM( t(deng.counts[,samples_subvector]),
#'                   K = 3, tol=0.1)
#'
#' names(fit_k3$clust_3)
#' omega <- fit_k3$clust_3$omega
#'
#' # make annotation matrix
#' annotation <- data.frame(
#'      sample_id = paste0("X", c(1:NROW(omega))),
#'      tissue_label = factor( as.character(deng.meta_data$cell_type[samples_subvector]),
#'              levels = rev( as.character(deng.meta_data$cell_type[samples_subvector])  ) )
#'      )
#' rownames(omega) <- annotation$sample_id
#' StructureGGplot(omega = omega,
#'                  annotation = annotation,
#'                  palette = RColorBrewer::brewer.pal(3, "Accent"),
#'                  yaxis_label = "development phase",
#'                  order_sample = TRUE,
#'                  axis_tick = list(axis_ticks_length = .1,
#'                                   axis_ticks_lwd_y = .1,
#'                                   axis_ticks_lwd_x = .1,
#'                                   axis_label_size = 7,
#'                                   axis_label_face = "bold"))
#'
#' @import ggplot2
#' @importFrom cowplot ggdraw panel_border switch_axis_position plot_grid
#' @import plyr
#' @import reshape2
#' @export

StructureGGplot <- function(omega, annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            sample_order_decreasing = TRUE,
                            split_line = list(split_lwd = 1,
                                              split_col = "white"),
                            plot_labels = TRUE,
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
    if (!is.data.frame(annotation))
        stop("annotation must be a data.frame")
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

                          is_single_sample <-
                                  ( length(temp_df) == nlevels(annotation$tissue_label)|
                                           is.null(dim(temp_df)) )
                          # find the dominant cluster in each sample
                          if ( is_single_sample ) {
                              each_sample_order <- which.max(temp_df)
                          } else {
                              each_sample_order <- apply(temp_df, 1, which.max)
                          }

                          # find the dominant cluster across samples
                          sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])

                          if (order_sample == TRUE & !is_single_sample) {
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
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 12)) +
        ggplot2::theme_update( panel.grid.minor.x = ggplot2::element_blank(),
                               panel.grid.minor.y = ggplot2::element_blank(),
                               panel.grid.major.x = ggplot2::element_blank(),
                               panel.grid.major.y = ggplot2::element_blank() )

    # inflat nubmers to avoid rounding errors
    value_ifl <- 10000

    # number of ticks for the weight axis, including 0 and 1
    ticks_number <- 6

    # set axis tick positions
    tissue_count <- table(droplevels(annotation$tissue_label))
    tissue_count_cumsum <- cumsum(table(droplevels(annotation$tissue_label)))

    tissue_names <- levels(droplevels(annotation$tissue_label))
    tissue_breaks <- sapply(1:length(tissue_count), function(i) {
        if (i == 1) {
            if (tissue_count[i] == 1) bk <- 1
            if (tissue_count[i] > 1)  bk <- (tissue_count_cumsum[i] - 0)/2
            return(bk)
        }
        if (i > 1) {
            if (tissue_count[i] == 1) bk_interval <- 1
            if (tissue_count[i] > 1 ) {
                bk_interval <- (tissue_count_cumsum[i] - tissue_count_cumsum[i-1])/2 }
            bk <- tissue_count_cumsum[i-1] + bk_interval
            return(bk)
        }
    })
    names(tissue_breaks) <- tissue_names

    # make ggplot
    a <- ggplot2::ggplot(df_mlt,
                         ggplot2::aes(x = df_mlt$document,
                                      y = df_mlt$value*10000,
                                      fill = factor(df_mlt$topic)) ) +
        ggplot2::xlab(yaxis_label) + ggplot2::ylab("") +
        ggplot2::scale_fill_manual(values = palette) +
        ggplot2::theme(legend.position = "right",
                       legend.key.size = ggplot2::unit(.2, "cm"),
                       legend.text = ggplot2::element_text(size = 5),
                       ##<-- TBD: center legend title
                       #              legend.title = element_text(hjust = 1),
                       axis.text = ggplot2::element_text(size = axis_tick$axis_label_size,
                                                         face = axis_tick$axis_label_face),
                       axis.ticks.y = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_y),
                       axis.ticks.x = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_x),
                       axis.ticks.length = ggplot2::unit(axis_tick$axis_ticks_length, "cm"),
                       title = ggplot2::element_text(size = 6) ) +
        ggplot2::ggtitle(figure_title) +
        ggplot2::scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                                     labels = seq(0, 1, 1/(ticks_number -1 ) ) ) +
        # Add tissue axis labels
        # ggplot2::scale_x_discrete(breaks = as.character(as.numeric(levels(df_mlt$document)[round(tissue_breaks)])),
        #                           labels = names(tissue_breaks)) +
        ggplot2::scale_x_discrete(breaks = as.character((levels(df_mlt$document)[round(tissue_breaks)])),
                                  labels = names(tissue_breaks)) +
        # Add legend title
        ggplot2::labs(fill = "Clusters") +
        ggplot2::coord_flip()


    # width = 1: increase bar width and in turn remove space
    # between bars
    b <- a + ggplot2::geom_bar(stat = "identity",
                               position = "stack",
                               width = 1)
    # sample labels option
    if (plot_labels == TRUE) {
        b
    } else {
        b <- b + theme(axis.text.y = element_blank())
    }
    
    # remove plot border
    b <- b + cowplot::panel_border(remove = TRUE)
    
    # Add demarcation
    b <- b + ggplot2::geom_vline(
        xintercept = cumsum(table(droplevels(annotation$tissue_label)))[
            -length(table(droplevels(annotation$tissue_label)))] + .5,
        col = split_line$split_col,
        size = split_line$split_lwd)
    b
    
    # if (!plot_labels) {
    #     b
    # } else {
    #     b <- cowplot::ggdraw(cowplot::switch_axis_position((b), axis = "y"))
    #     b
    # }
}
