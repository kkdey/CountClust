###--- GTEx ---###

# prepare data

omega <- read.table("../project/rdas/omega_cis_genes_0_1_2.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
tissue_labels <- vector("numeric", NROW(sample_labels))
tissue_labels <- sample_labels[ ,3]

# clean labels
tissue_labels[grep("Nucleus", tissue_labels)] <- "Brain -N. accumbens"
tissue_labels[grep("Putamen", tissue_labels)] <- "Brain -Putamen"
tissue_labels[grep("Caudate", tissue_labels)] <- "Brain -Caudate"
tissue_labels[grep("Gastroe", tissue_labels)] <- "Esophagus -Gastroesophageal Jn."
tissue_labels[grep("cingulate", tissue_labels)] <- "Brain - Anterior cortex (BA24)."
tissue_labels[grep("EBV", tissue_labels)] <- "Cells -EBV-lymphocytes"
tissue_labels[grep("Suprapubic", tissue_labels)] <- "Skin - Unexposed (Suprapubic)"
tissue_labels[grep("Lower Leg", tissue_labels)] <- "Skin - Sun Exposed (Lower Leg)"


# find sample orders in hierarchical clustering
docweights_per_tissue_mean <- apply(omega, 2, 
                                    function(x) { tapply(x, tissue_labels, mean) })
ordering <- heatmap(docweights_per_tissue_mean)$rowInd

# order tissue by hierarhical clustering results
tissue_levels_reordered <- unique(tissue_labels)[ordering]

# move Artery - Coronary next to the other Artery tissues
which(tissue_levels_reordered == "Artery - Coronary")
tissue_levels_reordered[45:50] <- c("Artery - Tibial", 
                                    "Artery - Coronary",
                                    "Esophagus - Muscularis",
                                    "Colon - Sigmoid",
                                    "Esophagus - Mucosa", 
                                    "Bladder")
# rearrange brain tissue order
tissue_levels_reordered[1:13] <- c(
    "Brain - Cerebellar Hemisphere", "Brain - Cerebellum",
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Anterior cortex (BA24).", "Brain - Frontal Cortex (BA9)",
    "Brain - Cortex", 
    "Brain - Hippocampus", "Brain - Substantia nigra", 
    "Brain - Amygdala", "Brain -Putamen", 
    "Brain -Caudate", 
    "Brain - Hypothalamus", "Brain -N. accumbens")


annotation <- data.frame(
    sample_id = paste0("X", 1:length(tissue_labels)),
    tissue_label = factor(tissue_labels,
                          levels = rev(tissue_levels_reordered ) ) )




##<----- Structure ggplot of all tissues

source("R/StructureGGplot.R")

# define colors of the clusers
# Joyce's color scheme
cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3"))

# Kushal's color scheme
cols2 <- c("red", "blue", "cornflowerblue", "black", "cyan", "darkblue",
           "brown4", "burlywood", "darkgoldenrod1", "darkgray", "deepskyblue",
           "darkkhaki", "firebrick", "darkorchid", "hotpink", "green",
           "magenta", "yellow", "azure1", "azure4")


source("R/StructureGGplot.R")
pdf("plots/gtex-figures/all-tissues-2.pdf", 
    height = 9, width=4)
StructureGGplot(omega = omega,
                annotation= annotation, 
                palette = cols1, 
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .1,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 5))
dev.off()





# # ##<----- Polar histogram of all tissues
# 
# source("R/polarHistogramStructure.R")
# 
# # rename columns to satisfy function input
# df_mlt <- reshape2::melt(t(omega_ordered))
# df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
#                                            "Var2" = "document"))
# df_mlt$document <- factor(df_mlt$document)
# df_mlt$topic <- factor(df_mlt$topic)
# head(df_mlt)
# df_plot <- df_mlt
# colnames(df_plot) <- c("score", "item", "value")
# df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
# df_plot$family <- factor(df_plot$family)
# df_plot$score <- factor(df_plot$score)
# df_plot$item <- factor(df_plot$item)
# levels(df_plot$family) <- levels(annotation$tissue_label)
# 
# 
# # being preparing for computing on midway
# save(polarHistogramStructure,
#      cols1, cols2,
#      df_plot, file = "rdas/rda-for-midway-all.rda")
# 
# # the following code was run on midway to make the plot
# # run large interactive job
# # gsinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000
# pdf("main-figure-polar-histogram-all.pdf",
#     height = 12, width = 12)
# polarHistogramStructure(df_plot, palette = cols1,
#                        outerRadius = 1.8,
#                        innerRadius = 0.3,
#                        familyLabelDistance = 2, # same metric as outerRadius
#                        binSize = 1,
#                        spaceFamily = 4,
#                        circleProportion = 0.90,
#                        familyLabels = TRUE)
# dev.off()
# 
# 


##<- thinned version


# prepare data

omega <- read.table("../project/rdas/omega_cis_genes_0_0001.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
tissue_labels <- vector("numeric", NROW(sample_labels))
tissue_labels <- sample_labels[ ,3]



# # move Artery - Coronary next to the other Artery tissues
# which(tissue_levels_reordered == "Artery - Coronary")
# tissue_levels_reordered[45:50] <- c("Artery - Tibial", 
#                                     "Artery - Coronary",
#                                     "Esophagus - Muscularis",
#                                     "Colon - Sigmoid",
#                                     "Esophagus - Mucosa", 
#                                     "Bladder")
# 
annotation_thinned <- data.frame(
    sample_id = paste0("X", 1:length(tissue_labels)),
    tissue_label = annotation$tissue_label ) 


source("R/StructureGGplot.R")

# define colors of the clusers
# Joyce's color scheme
cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3"))

# Match color with the unthinned version
cols_thinned <- cols1
cols_thinned[15] <- cols1[14];
cols_thinned[14] <- cols1[13];
cols_thinned[10] <- cols1[10];
cols_thinned[3] <- cols1[4];
cols_thinned[7] <- cols1[7];
cols_thinned[8] <- cols1[9];
cols_thinned[11] <- cols1[12];
cols_thinned[12] <- cols1[11];
cols_thinned[1] <- cols1[1];
cols_thinned[2] <- cols1[2];
cols_thinned[4] <- cols1[3];
cols_thinned[5] <- cols1[8];
cols_thinned[6] <- cols1[6];
cols_thinned[9] <- cols1[5];
cols_thinned[13] <- cols1[15];



pdf("plots/gtex-figures/all-tissues-thinned.pdf", 
    height = 9, width=4)
StructureGGplot(omega = omega,
                annotation= annotation, 
                palette = cols_thinned, 
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .1,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 5))
dev.off()







##<-- brain tissue plots


omega <- read.table("../project/rdas/omega_cis_genes_brain.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))
head(omega)


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
brain_labels <- sample_labels[grep("Brain", sample_labels[,3]), 3]

# assign tissue labels
rownames(omega) <- paste0("X", 1:length(brain_labels))
annotation <- data.frame(
    sample_id = paste0("X", 1:length(brain_labels)),
    tissue_label = factor(brain_labels,
                          levels = rev(c("Brain - Cerebellar Hemisphere",
                                     "Brain - Cerebellum",
                                     "Brain - Spinal cord (cervical c-1)",
                                     "Brain - Anterior cingulate cortex (BA24)",
                                     "Brain - Frontal Cortex (BA9)",
                                     "Brain - Cortex",
                                     "Brain - Hippocampus",
                                     "Brain - Substantia nigra",
                                     "Brain - Amygdala",
                                     "Brain - Putamen (basal ganglia)",
                                     "Brain - Caudate (basal ganglia)",
                                     "Brain - Hypothalamus",
                                     "Brain - Nucleus accumbens (basal ganglia)") ) ) )
                                     
# define colors of the clusers
cols <- c("blue", "darkgoldenrod1", "cyan", "red")

##<-- make barplot
source("R/StructureGGplot.R")

# subset_example <- which(annotation$tissue_label %in% 
#     levels(annotation$tissue_label)[1:2] )
pdf("plots/gtex-figures/brain-barplot.pdf",
    height = 4, width = 4)
StructureGGplot(omega = omega,
                annotation= annotation, 
                palette = cols, 
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .4,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 3,
                                 axis_label_face = "bold"))
StructureGGplot(omega = omega,
                annotation= annotation, 
                palette = cols, 
                yaxis_label = "",
                order_sample = TRUE,
                split_line = list(split_lwd = .4,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 3))
dev.off()



##<-- make polar histogram
source("R/polarHistogramStructure.R")

# rename columns to satisfy function input
df_mlt <- reshape2::melt(t(omega))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))

df_mlt$document <- factor(df_mlt$document)
df_mlt$topic <- factor(df_mlt$topic)
head(df_mlt)
df_plot <- df_mlt
colnames(df_plot) <- c("score", "item", "value")
df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
df_plot$family <- factor(df_plot$family)
df_plot$score <- factor(df_plot$score)
df_plot$item <- factor(df_plot$item)
levels(df_plot$family) <- levels(annotation$tissue_label)


# preparing for computing on midway
save(polarHistogramStructure,
     cols1, cols2, cols,
     df_plot, file = "rdas/rda-for-midway.rda")


# the following code was run on midway to make the plot
# run large interactive job
# 
# sinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000
pdf("main-figure-polar-histogram.pdf",
    height = 8, width = 8)
    polarHistogramStructure(
        df = df_plot, 
        palette = cols,
        outerRadius = 1.8,
        innerRadius = 0.3,
        familyLabelDistance = 2, # same metric as outerRadius
        binSize = 1,
        spaceFamily = 4,
        circleProportion = 0.90,
        familyLabels = TRUE)
dev.off()


 





#####<------ Deng et al. 6 clusters ------>######

source("R/StructureGGplot.R")

# load the previously analyzed results
load("rdas/deng_topic_fit.rda")

# extract the omega matrix: membership weights of each cell
names(Topic_clus_list)
str(Topic_clus_list$clust_6)
omega <- Topic_clus_list$clust_6$omega

# import embryonl labels
embryo_label <- read.table("rdas/cell_labels_phase_embryo.txt",
                           quote = "\"", 
                           header = TRUE,
                           stringsAsFactors = FALSE)$x
head(embryo_label, 20)
table(embryo_label)
stopifnot(length(embryo_label) == NROW(omega))


# make annotation matrix for the plot of all tissues
# sample_id has to be unique
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(rownames(omega),
                        levels = rev( c("zy", "early2cell", "mid2cell", "late2cell",
                                       "4cell", "8cell", "16cell", "earlyblast",
                                       "midblast", "lateblast") ) ) )

# make annotation for early stage plot
# sample_id has to be unique
annotation_embryo <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(embryo_label,
      levels = rev( c("zy_.",
                      paste("early2cell",c("0r", c(1:3)), sep = "_"), 
                      paste("mid2cell",c("0r", c(3:7)), sep = "_"),
                      paste("late2cell",c("0r", c(5:9)), sep = "_"),
                      paste("4cell",c("0r", c(1:4)), sep = "_"),
                      paste("8cell",c("0r", c(1,2,5,8)), sep = "_"),
                      paste("16cell",c("0r", c(1,4,5,6)), sep = "_"),
                      paste("earlyblast",c("0r", c(2:4)), sep = "_"),
                      paste("midblast",c("0r", c(1:3)), sep = "_"),
                      paste("lateblast",c("0r", c(1:3)), sep = "_") ) ) ) ) 

      
# after extracting tissue type of each sample
# recode each sample to have unique rownames
rownames(omega) <- paste0("X", annotation$sample_id)

pdf(file = "plots/deng-figures/deng-ggplot-all.pdf", 
    height = 4, width = 3)
StructureGGplot(omega = omega, 
               annotation = annotation,
               palette = RColorBrewer::brewer.pal(8, "Accent"),
               figure_title = "",
               yaxis_label = "Cell type",
               sample_order_decreasing = FALSE)
dev.off()


pdf(file = "plots/deng-figures/deng-ggplot-embryo.pdf", 
    height = 4, width = 3)
embryo_plot <- which(
    annotation_embryo$tissue_label %in% 
        levels(annotation_embryo$tissue_label)[c(13:16, 18:21, 23:26)])
StructureGGplot(omega = omega[embryo_plot, ], 
                annotation = annotation_embryo[embryo_plot, ],
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                figure_title = "",
                yaxis_label = "Cell type",
                sample_order_decreasing = FALSE)
dev.off()


# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
cell_labels_count <- table(droplevels(annotation_embryo$tissue[embryo_plot]) )
png(file = "plots/deng-figures/labels-early.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()

# 
# 
# ##<-- Display for some early stages
# 
# # make the re-ordered dataframe
# df_ord <- do.call(rbind,
#       lapply(5:7, function(ii) {
#           temp_label <- levels(annotation$tissue_label)[ii]
#           temp_df <- omega[which(annotation$tissue_label == temp_label), ]
# 
#           # find the dominant cluster in each sample
#           each_sample_order <- apply(temp_df, 1, which.max)
#           
#           # find the dominant cluster across samples
#           sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
#           
#           # reorder the matrix
#           temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
#           
#           temp_df_ord
#       }) )
# dim(df_ord)
# 
# #df_mlt <- reshape2::melt(t(df_ord))
# df_mlt <- reshape2::melt(t(df_ord))
# df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
#                                            "Var2" = "document"))
# head(df_mlt)
# 
# # set blank background
# theme_set(theme_bw(base_size = 12)) +
#     theme_update( panel.grid.minor.x = element_blank(),
#                   panel.grid.minor.y = element_blank(),
#                   panel.grid.major.x = element_blank(),
#                   panel.grid.major.y = element_blank() )
# 
# value_ifl <- 10000
# ticks_number <- 6
# # Use RColorBrewer color
# # http://bxhorn.com/rcolorbrewer-palettes/
# a <- ggplot(df_mlt, 
#             aes(x = document, y = value*10000, fill = factor(topic)) ) + 
#     xlab("Cell types") + ylab("") +
#     scale_fill_brewer(palette = "Accent") +
#     theme(legend.position = "none",
#           axis.text = element_text(size = 4),
#           title = element_text(size = 6)) +
#     ggtitle("STRUCTURE plot by developmental phase") + 
#     scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
#                         labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
#     coord_flip() 
# 
# # width = 1: increase bar width to 1 so that there's no 
# # space between bars
# b <- a + geom_bar(stat = "identity", 
#                   position = "stack", 
#                   width = 1)
# b <- b + panel_border(remove = TRUE)
# b <- ggdraw(switch_axis_position((b), axis = "y"))
# cowplot::plot_grid(b)
# save_plot("../../../count-clustering/project/plots/deng-figures/deng-ggplot-early.png", 
#           b, base_height = 4, base_width = 2)
# #ggplotly(b)










######<--- Jaitin et al. 7 clusters --->######
# see here for analysis steps
# http://stephenslab.github.io/count-clustering/project/src/jaitin_structure_genes.html

# load the previously analyzed results
topic_fit <- readRDS("rdas/MouseJaitinSpleen-topicFit.rds")

# extract the omega matrix: membership weights of each cell
names(topic_fit$clust_7)
omega <- topic_fit$clust_7$omega

# load phenotype information
library(singleCellRNASeqMouseJaitinSpleen)
counts <- exprs(MouseJaitinSpleen)
meta_data <- pData(MouseJaitinSpleen)
gene_names <- rownames(counts)

# follow previous anaysis steps to extract
# amplification batch information
# 
# exclude ERCC genes 
ENSG_genes_index <- grep("ERCC", gene_names, invert = TRUE)

# expression matrix without ENSG genes
counts_ensg <- counts[ENSG_genes_index, ]
filter_genes <- c("M34473","abParts","M13680","Tmsb4x",
                  "S100a4","B2m","Atpase6","Rpl23","Rps18",
                  "Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4",
                  "Rps26","EF437368") 
fcounts <- counts_ensg[ -match(filter_genes, rownames(counts_ensg)), ]
sample_counts <- colSums(fcounts)

filter_sample_index <- which(meta_data$number_of_cells == 1 & 
                                 meta_data$group_name == "CD11c+" & 
                                 sample_counts > 600)

# make filterd phenotype data
meta_data_filtered <- meta_data[filter_sample_index, ]
stopifnot(dim(meta_data_filtered)[1] == dim(omega)[1])

# load packages
library(ggplot2)
library(plotly)
library(reshape2)
library(cowplot)
source("R/StructureGGplot.R")


###<--- Order by amplification batch

# making annotation matrix
annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega)) ),
    tissue_label = factor(amp_batch,
                          levels = rev(sort(unique(amp_batch))) ) )
pdf("plots/jaitin-figures/amplification.pdf", 
     height = 4, width = 3)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(9, "Set1"),
                yaxis_label = "Amplification batch",
                order_sample = TRUE)
dev.off()

# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
png(file = "plots/jaitin-figures/amplfication-labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()




###<--- Order by sequencing batch

# making annotation matrix
seq_batch <- as.numeric(meta_data_filtered[ , "sequencing_batch"])
annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega)) ),
    tissue_label = factor(seq_batch,
                          levels = rev(sort(unique(seq_batch))) ) )

pdf("plots/jaitin-figures/sequencing.pdf", 
    height = 4, width = 3)
StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(9, "Set1"),
                yaxis_label = "Amplification batch",
                order_sample = TRUE)
dev.off()

# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
png(file = "plots/jaitin-figures/sequencing-labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()

