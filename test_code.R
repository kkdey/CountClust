
###
library(devtools)
install_github('kkdey/CountClust')
library(CountClust)

###  Reading the GTEX V6 data

brain_data <- t(data.frame(data.table::fread("test/cis_gene_expression_brain.txt"))[,-(1:2)]);
brain_data <- brain_data[,1:1000];
dim(brain_data)

### Reading the sample metadata that you wish toorganize the Structure plot

brain_metadata <- cbind.data.frame(read.table("test/metadata_brain.txt")[,3]);
colnames(brain_metadata) <- "brain_labs"
dim(brain_metadata)

## Using the StructureObj to fit the Structure model (using the maptpx package due to Matt Taddy) 

if(!dir.exists("test/Structure")) dir.create("test/Structure")
StructureObj(brain_data, nclus_vec=3, samp_metadata = brain_metadata, tol=0.1, batch_lab = NULL,
             plot=TRUE, path_rda="test/brain_structure.rda", path_struct="test/Structure")

#The rda file contains the W and T files, 
brain_structure <- get(load("test/brain_structure.rda"));
StructureObj_omega(brain_structure[[1]]$omega, samp_metadata = brain_metadata, 
                          batch_lab = NULL, path_struct="test/Structure",
                          control=list(cex.axis=1, struct.width=500, struct.height=500))