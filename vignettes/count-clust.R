## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,tidy.opts=list(width.cutoff=30),dev="png",fig.show="hide",
               fig.width=4,fig.height=7,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------------------
BiocStyle::latex()

## ----options, results="hide", echo=FALSE--------------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----install_countclust, eval=FALSE-------------------------------------------
#  BiocInstaller::biocLite("CountClust")

## ----install_countclust_github, eval=TRUE-------------------------------------
library(devtools)
install_github('kkdey/CountClust')

## ----install_github_maptpx, eval=TRUE-----------------------------------------
library(devtools)
install_github("TaddyLab/maptpx")

## ----load_countclust, cache=FALSE, eval=TRUE, warning=FALSE-------------------
library(CountClust)

## ----data_install, eval=TRUE--------------------------------------------------
devtools::install_github('kkdey/singleCellRNASeqMouseDeng2014')
devtools::install_github('kkdey/GTExV6Brain')

## ----data_load_deng, eval=FALSE-----------------------------------------------
#  library(singleCellRNASeqMouseDeng2014)
#  deng.counts <- exprs(Deng2014MouseESC)
#  deng.meta_data <- pData(Deng2014MouseESC)
#  deng.gene_names <- rownames(deng.counts)

## ----data_load_gtex, eval=FALSE-----------------------------------------------
#  library(GTExV6Brain)
#  gtex.counts <- exprs(GTExV6Brain)
#  gtex.meta_data <- pData(GTExV6Brain)
#  gtex.gene_names <- rownames(gtex.counts)

## ----topic_fit_gtex, eval=FALSE-----------------------------------------------
#  FitGoM(t(gtex.counts),
#              K=4, tol=0.1,
#              path_rda="../data/GTExV6Brain.FitGoM.rda")

## ----topic_fit_deng, eval=FALSE-----------------------------------------------
#  FitGoM(t(deng.counts),
#              K=2:7, tol=0.1,
#              path_rda="../data/MouseDeng2014.FitGoM.rda")

## ----plot_topic_deng_annot,eval=TRUE, warning=FALSE---------------------------
data("MouseDeng2014.FitGoM")
names(MouseDeng2014.FitGoM$clust_6)
omega <- MouseDeng2014.FitGoM$clust_6$omega


annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(rownames(omega),
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                         "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

## ----plot_topic_deng,eval=TRUE, warning=FALSE, fig.show="asis", dpi=144, fig.width=3, fig.height=7, out.width="3in", out.height="7in"----
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))



## ----gtex_annot, eval=FALSE---------------------------------------------------
#  data("GTExV6Brain.FitGoM")
#  omega <- GTExV6Brain.FitGoM$omega;
#  dim(omega)
#  colnames(omega) <- c(1:NCOL(omega))
#  
#  tissue_labels <- gtex.meta_data[,3];
#  
#  
#  annotation <- data.frame(
#      sample_id = paste0("X", 1:length(tissue_labels)),
#      tissue_label = factor(tissue_labels,
#                            levels = rev(unique(tissue_labels) ) ) );
#  
#  cols <- c("blue", "darkgoldenrod1", "cyan", "red")

## ----plot_topic_gtex,eval=FALSE, warning=FALSE, fig.show="asis", dpi=144, fig.width=5, fig.height=7, out.width="5in", out.height="7in"----
#  StructureGGplot(omega = omega,
#                  annotation= annotation,
#                  palette = cols,
#                  yaxis_label = "",
#                  order_sample = TRUE,
#                  split_line = list(split_lwd = .4,
#                                    split_col = "white"),
#                  axis_tick = list(axis_ticks_length = .1,
#                                   axis_ticks_lwd_y = .1,
#                                   axis_ticks_lwd_x = .1,
#                                   axis_label_size = 7,
#                                   axis_label_face = "bold"))

## ----extract_features_deng, eval=FALSE, warning=FALSE-------------------------
#  theta_mat <- MouseDeng2014.FitGoM$clust_6$theta;
#  top_features <- ExtractTopFeatures(theta_mat, top_features=100,
#                                     method="poisson", options="min");
#  gene_list <- do.call(rbind, lapply(1:dim(top_features)[1],
#                          function(x) deng.gene_names[top_features[x,]]))

## ----top_genes_clusters_deng, eval=FALSE--------------------------------------
#  library(xtable)
#  xtable(gene_list[,1:5])

## ----extract_features_gtex, eval=FALSE, warning=FALSE-------------------------
#  theta_mat <- GTExV6Brain.FitGoM$theta;
#  top_features <- ExtractTopFeatures(theta_mat, top_features=100,
#                                     method="poisson", options="min");
#  gene_list <- do.call(rbind, lapply(1:dim(top_features)[1],
#                          function(x) gtex.gene_names[top_features[x,]]))

## ----top_genes_clusters_gtex, eval=FALSE--------------------------------------
#  library(xtable)
#  xtable(gene_list[,1:3])

## ----data_install_jaitin, echo=TRUE, eval=FALSE-------------------------------
#  devtools::install_github('jhsiao999/singleCellRNASeqMouseJaitinSpleen')

## ----data_load_jaitin, echo=TRUE, eval=FALSE----------------------------------
#  library(singleCellRNASeqMouseJaitinSpleen)
#  jaitin.counts <- exprs(MouseJaitinSpleen)
#  jaitin.meta_data <- pData(MouseJaitinSpleen)
#  jaitin.gene_names <- rownames(jaitin.counts)

## ----non_ercc, eval=FALSE, echo=TRUE------------------------------------------
#  ENSG_genes_index <- grep("ERCC", jaitin.gene_names, invert = TRUE)
#  jaitin.counts_ensg <- jaitin.counts[ENSG_genes_index, ]
#  filter_genes <- c("M34473","abParts","M13680","Tmsb4x",
#                    "S100a4","B2m","Atpase6","Rpl23","Rps18",
#                    "Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4",
#                    "Rps26","EF437368")
#  fcounts <- jaitin.counts_ensg[ -match(filter_genes, rownames(jaitin.counts_ensg)), ]
#  sample_counts <- colSums(fcounts)
#  
#  filter_sample_index <- which(jaitin.meta_data$number_of_cells == 1 &
#                                jaitin.meta_data$group_name == "CD11c+" &
#                                   sample_counts > 600)
#  fcounts.filtered <- fcounts[,filter_sample_index];
#  

## ----metadata, eval=FALSE, echo=TRUE------------------------------------------
#  jaitin.meta_data_filtered <- jaitin.meta_data[filter_sample_index, ]

## ----topic_fit_jaitin, eval=FALSE, echo=TRUE----------------------------------
#  StructureObj(t(fcounts),
#              nclus_vec=7, tol=0.1,
#               path_rda="../data/MouseJaitinSpleen.FitGoM.rda")

## ----plot_topic_annot, eval=FALSE, echo=TRUE----------------------------------
#  data("MouseJaitinSpleen.FitGoM")
#  names(MouseJaitinSpleen.FitGoM$clust_7)
#  omega <- MouseJaitinSpleen.FitGoM$clust_7$omega
#  
#  amp_batch <- as.numeric(jaitin.meta_data_filtered[ , "amplification_batch"])
#  annotation <- data.frame(
#      sample_id = paste0("X", c(1:NROW(omega)) ),
#      tissue_label = factor(amp_batch,
#                            levels = rev(sort(unique(amp_batch))) ) )

## ----plot_topic, eval=FALSE, echo=TRUE, warning=FALSE, fig.show="asis", dpi=144, fig.width=3, fig.height=7, out.width="3in", out.height="7in"----
#  StructureGGplot(omega = omega,
#                  annotation = annotation,
#                  palette = RColorBrewer::brewer.pal(9, "Set1"),
#                  yaxis_label = "Amplification batch",
#                  order_sample = FALSE,
#                  axis_tick = list(axis_ticks_length = .1,
#                                   axis_ticks_lwd_y = .1,
#                                   axis_ticks_lwd_x = .1,
#                                   axis_label_size = 7,
#                                   axis_label_face = "bold"))
#  

## ----batch_correct, eval=FALSE, echo=TRUE-------------------------------------
#  batchcorrect.fcounts <- BatchCorrectedCounts(t(fcounts.filtered),
#                                               amp_batch, use_parallel = TRUE);
#  dim(batchcorrect.fcounts)

## ----session_info, eval=TRUE--------------------------------------------------
sessionInfo()

