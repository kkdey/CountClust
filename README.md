# CountClust
A R package for counts clustering

To download and install this package, 

```
install.packages("devtools")
library(devtools)
install_github("kkdey/CountClust")
```

Load the package in R

```
library(CountClust)
```

Load the counts data (samples along the rows, variables or features along the columns) in R. Here we use a test data
comprising of the brain samples in GTEX V6 data

```
data <- t(data.frame(data.table::fread("test/cis_gene_expression_brain.txt"))[,-(1:2)]);
data <- data[,1:1000];
dim(data)
```
Next we load the brain metadata 

```
metadata <- cbind.data.frame(read.table("test/metadata_brain.txt")[,3]);
colnames(metadata) <- "labs"
dim(metadata)
```
We  apply the StructureObj function to fit the topic model (due to the **maptpx** package of Matt Taddy) 
and plot the Structure plot. 

```
if(!dir.exists("test/Structure")) dir.create("test/Structure")
StructureObj(data, nclus_vec=3, samp_metadata = metadata, tol=0.1, batch_lab = NULL,
             plot=TRUE, path_rda="test/topics_data.rda", path_struct="test/Structure")
```

This function will output two things:

-  a rda file containing all topic model information named "test/topics_data.rda" 
-  The Structure plot for K=3 arranged by the metadata provided in test/Structure. If samp_metadata is a matrix 
   with several columns, then Structure plots arranged by each column in samp_metadata will be created.

The rda file contains the W (the topic proportion matrix) and T (topic distribution) files, so if you need to change
the Structure plot, you do not need to rerun the StructureObj fit again and you can just load the W matrix in rda 
file into out StrutureObj_omega() function.

```
brain_structure <- get(load("test/topics_data.rda"));
StructureObj_omega(brain_structure[[1]]$omega, samp_metadata = brain_metadata, 
                          batch_lab = NULL, path_struct="test/Structure",
                          control=list(cex.axis=1, struct.width=500, struct.height=500))
```

Then go to `test/Structure` and you will be able to see the Structure plot for K=3 with *cex.axis* modified


