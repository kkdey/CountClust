
###
library(devtools)
install_github('kkdey/CountClust')
library(CountClust)

###  Reading the GTEX V6 data

library(data.table)
brain_data <- data.frame(data.table::fread("test/cis_gene_expression_brain.txt"))[,-(1:2)];
dim(brain_data)

### Applying FLASH

sqr_brain_data <- sqrt(brain_data);
out <- flash_VEM(as.matrix(sqr_brain_data[1:10,1:20]));

samples_id_mat <- read.table("test/samples_id.txt");
samples_id <- samples_id_mat[which(samples_id_mat[,2]=="Brain"),3];

docweights <- out$l;
K=dim(docweights)[2];
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");
brain_ids <- samples_id;
ordering <- order(brain_ids);
samples_id_ordered <- brain_ids[ordering];
docweights_ordering <- docweights[ordering,];
#png(filename=paste0('../plots/GTEX_V6_brain_thin_',0,'.png'),width=700,height=300)
par(mar=c(14,2,2,1))
barplot(t(docweights_ordering),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,cex.axis=1,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=0.8);