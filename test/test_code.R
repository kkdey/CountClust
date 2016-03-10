
###
library(devtools)
install_github('kkdey/CountClust')
library(CountClust)

###  Reading the GTEX V6 data

brain_data <- data.frame(data.table::fread("cis_gene_expression_brain.txt"))[,-(1:2)];
dim(brain_data)

### Reading the sample metadata that you wish toorganize the Structure plot

brain_labels <- as.vector((read.table("metadata_brain.txt")[,3]));


## Using the StructureObj to fit the Structure model (using the maptpx package due
## to Matt Taddy) for clusters 2 and 3

StructureObj(brain_data, nclus_vec=2:3, tol=0.1, path_rda="brain_structure.rda")

#The rda file contains the W and T files,

## Load the rda file
Topic_clus <- get(load("topics_data.rda"));
omega <- Topic_clus$clust_3$omega;

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


cols <- c("blue", "darkgoldenrod1", "cyan")


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
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


