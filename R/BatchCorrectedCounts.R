#' @title Obtain Batch effect Corrected counts
#'
#' @description This function first converts counts data to log CPM data ,
#'    then apply a linear model with the batch effect as a factor. We take
#'    the sum of intercept, residuals and mean batch effect across all the
#'    batches and then inverse transform it back to counts to get rid of
#'    batch effects.
#'
#' @param data count matrix, with samples along the rows and features
#'             along the columns.
#' @param batch_lab batch label vector.
#' @param use_parallel if TRUE, we do a parallel analysis over features,
#'        else serial application.
#'
#' @return Returns a counts data. with same dimension as the input data,
#'         but which is corrected for batch_lab.
#'
#' @keywords counts data, batch effect
#'
#' @examples
#'
#' # Simulation example
#' N=500;
#' K=4;
#' G=100;
#' Label.Batch=c(rep(1,N/4),rep(2,N/4),rep(3,N/4),rep(4,N/4));
#' alpha_true=matrix(rnorm((K)*G,0.5,1),nrow=(K));
#' library(gtools)
#' tt <- 10;
#' omega_true = matrix(rbind(rdirichlet(tt*10,c(3,4,2,6)),
#'                          rdirichlet(tt*10,c(1,4,6,3)),
#'                          rdirichlet(tt*10,c(4,1,2,2)),
#'                          rdirichlet(tt*10,c(2,6,3,2)),
#'                          rdirichlet(tt*10,c(3,3,5,4))), nrow=N);
#' B=max(Label.Batch);
#' sigmab_true=2;
#' beta_true=matrix(0,B,G);
#' for(g in 1:G)
#' {
#'     beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
#' }
#' read_counts=matrix(0,N,G);
#' for(n in 1:N){
#'     for(g in 1:G)
#'     {
#'         read_counts[n,g]=rpois(1, omega_true[n,]%*%exp(alpha_true[,g]
#'                                                       + beta_true[Label.Batch[n],g]));
#'    }
#'}

#' K <- 4
#' barplot(t(omega_true),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",3),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
#' topic_clus <- maptpx::topics(read_counts, K=4, tol=0.1)
#' barplot(t(topic_clus$omega),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",3),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
#'
#' batchcorrect_counts <- BatchCorrectedCounts(read_counts, Label.Batch, use_parallel=TRUE)
#' topic_clus <- maptpx::topics(batchcorrect_counts, K=4, tol=0.1)
#' library(permute);
#' library("BioPhysConnectoR");
#' perm_set=rbind(1:K,allPerms(1:K));
#' diff=array(0,dim(perm_set)[1]);
#' for (p in 1:dim(perm_set)[1])
#' {
#'     temp=topic_clus$omega[,perm_set[p,]];
#'     diff[p]=fnorm(temp,omega_true);
#' }
#' p_star=which(diff==min(diff));
#' docweights=topic_clus$omega[,perm_set[p_star,]];
#' barplot(t(docweights),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",3),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
#'
#' @importFrom gtools rdirichlet
#' @importFrom  limma voom
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @export
#'
BatchCorrectedCounts <- function(data, batch_lab, use_parallel = TRUE)
{
    trans_data <- voom2(data);
    lib_size <- rowSums(data);
    batch_lab <- as.factor(batch_lab)
    if(use_parallel){
        batch_removed_counts_mean <-
            do.call(cbind,
                    parallel::mclapply(1:dim(trans_data)[2], function(g)
                    {
                        out <- lm(trans_data[,g] ~  batch_lab,
                                  contrasts =  list(batch_lab="contr.sum"))
                        return(round((2^{out$coefficients[1] + out$residuals - 6*log(10, base=2)})*(lib_size+1) - 0.5))
                    }, mc.cores=parallel::detectCores()));
    }

    if(!use_parallel){
        batch_removed_counts_mean <-
            do.call(cbind, lapply(1:dim(trans_data)[2], function(g)
            {
                out <- lm(trans_data[,g] ~  batch_lab,
                            contrasts =  list(batch_lab="contr.sum"))
                return(round((2^{out$coefficients[1] + out$residuals - 6*log(10, base=2)})*(lib_size+1) - 0.5))
                #return(round(exp((out$coefficients[1] + out$residuals)/6)*(lib_size+1)-0.4));
            }));
    }

    if (dim(batch_removed_counts_mean)[2]!=dim(data)[2])
        stop("The batch corrected data is not of same dimension as the counts data : try changing use_parallel")
    batch_corrected_counts <- round(batch_removed_counts_mean);
    rownames(batch_corrected_counts) = rownames(data);
    colnames(batch_corrected_counts) = colnames(data);
    return(batch_corrected_counts)
}

voom2 <- function(counts){
    libsize.mat <- rep.col(rowSums(counts), dim(counts)[2]);
    voom.out <- log((counts+0.5), base=2) - log((libsize.mat+1), base=2)+ 6* log(10, base=2);
    return(voom.out)
}
rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
